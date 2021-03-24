/* -*- c++ -*- */
/*
 * Copyright 2018, 2019, 2020 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 * Copyright 2021 Jacob Gilbert
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_IMPL_H
#define INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_IMPL_H

#include <gnuradio/fft/fft.h>
#include <fhss_utils/fft_burst_tagger.h>
#include <mutex>
//#define __USE_MKL__
#ifdef __USE_MKL__
#include "mkl_dfti.h"
#endif

/*
 * DO_TIMER is used for profiling timer
 *
 * WARNING: Use of profiling timers has performance implications, they are only complied
 * if this option is enabled.
 */
//#define DO_TIMER
#include <chrono>
using namespace std::chrono;
using std::chrono::high_resolution_clock;

namespace gr {
namespace fhss_utils {

struct peak {
    uint32_t bin;
    float relative_magnitude;
};

class timer
{
public:
    timer() : total(0) {}

#ifndef DO_TIMER
    void start() {}
    void end() {}
#else
    void start() { start_time = high_resolution_clock::now(); }
    void end()
    {
        high_resolution_clock::time_point end_time = high_resolution_clock::now();
        duration<double> span = duration_cast<duration<double>>(end_time - start_time);
        total += span.count();
    }
#endif
    double elapsed() { return total; }
    void reset() { total = 0; }

private:
    double total;
    high_resolution_clock::time_point start_time;
}; // end class timer

/**
 * Calculates a moving average utilizing circular buffer, used in this block for the noise
 * estimate in each bin.
 */
class moving_average
{
public:
    /**
     * Constructor
     *
     * @param size - size in samples of the moving average circular buffer
     */
    moving_average(size_t size) : N(size)
    {
        sum = current_index = 0.0;
        pp = 0.0;
        hist.resize(N);
        memset(&hist[0], 0, sizeof(float) * N);
    }

    /**
     * Add a value to circular buffer
     *
     * @param p - value to add
     * @return float - current accumulator value
     */
    float add(float p)
    {
        sum += pp - hist[current_index];
        hist[current_index++] = pp;
        pp = p;
        if (current_index == N)
            current_index = 0;
        return sum / N;
    }

    /**
     * Print the contents of the a given moving average buffer in a python-friendly
     * way. This is a debug function only that is not normally used
     */
    std::stringstream print(void)
    {
        std::stringstream sout;
        sout << "[" << hist[0];
        for (auto ii = 1; ii < N; ii++)
            sout << ", " << hist[ii];
        sout << "]";
        return sout;
    }

private:
    std::vector<float> hist;
    size_t current_index;
    size_t N;
    float sum;
    float pp;   // delay noise floor sum by one FFT
}; // end class moving_average

/*
 * A pre burst is energy that has been detected as being above the configured threshold,
 * but is not yet sustained long enough to be considered a burst.
 */
struct pre_burst {
    uint64_t start;
    uint64_t peak_count;
    uint64_t id;
    uint64_t thresh_count;
    int center_bin;
    int start_bin;
    int stop_bin;
    float max_relative_magnitude;
};

/*
 * A burst is sustained energy above the configured threshold and will result in a tag
 * being emitted.
 */
struct burst {
    uint64_t start;
    uint64_t stop;
    uint64_t last_active;
    uint64_t id;
    int center_bin;
    int start_bin;
    int stop_bin;
    float magnitude;
    float bandwidth;
    float center_freq;
    float noise;
    bool valid;
};

struct owners {
    owners() { clear(); }
    uint64_t ids[4];
    uint64_t _size;
    uint64_t uid;

    size_t size() { return _size; }

    void clear()
    {
        ids[0] = ids[1] = ids[2] = ids[3] = (uint64_t)-1;
        _size = 0;
    }

    /* return of zero is success, otherwise error */
    int push_back(uint64_t id)
    {
        if (_size > 3) {
            return 1;
            // printf("Owners::push_back - Trying to add too many points for  bin id=%zu,
            // size=%zu, burst=%zu\n", uid, _size, id); printf("Owner bins are %zu, %zu,
            // %zu %zu\n", ids[0], ids[1], ids[2], ids[3]);
        }
        ids[_size] = id;
        _size++;
        return 0;
    }

    int update(uint64_t oldID, uint64_t newID)
    {
        if (ids[0] == oldID)
            ids[0] = newID;
        else if (ids[1] == oldID)
            ids[1] = newID;
        else if (ids[2] == oldID)
            ids[2] = newID;
        else if (ids[3] == oldID)
            ids[3] = newID;
        else
            return 1;
        // printf("Owners::Update - Couldn't find id to update. This should never
        // happen\n");
        return 0;
    }

    int erase(uint64_t id)
    {
        _size--;
        if (ids[0] == id) {
            if (_size) {
                ids[0] = ids[_size];
                ids[_size] = (uint64_t)-1;
            } else
                ids[0] = (uint64_t)-1;
        } else if (ids[1] == id) {
            if (_size > 1) {
                ids[1] = ids[_size];
                ids[_size] = (uint64_t)-1;
            } else
                ids[1] = (uint64_t)-1;
        } else if (ids[2] == id) {
            if (_size > 2) {
                ids[2] = ids[_size];
                ids[_size] = (uint64_t)-1;
            } else
                ids[2] = (uint64_t)-1;
        } else if (ids[3] == id) {
            ids[3] = (uint64_t)-1;
        } else {
            _size++; // Increment size because we didn't remove anything.
            return 1;
            // printf("Owners::Remove - Couldn't find id to erase. This should never
            // happen - bin id = %zu, burst=%zu\n", uid, id);
        }
        return 0;
    }
};

/**
 * Implementing class of fft_burst_tagger
 */
class fft_burst_tagger_impl : public fft_burst_tagger
{
private:
    bool d_history_primed;
    bool d_debug;

    int d_fft_size;
    int d_fine_fft_size;
    int d_lookahead;
    int d_burst_pre_len;
    int d_history_size;
    int d_burst_width;
    int d_history_index;
    int d_burst_post_len;
    int d_max_bursts;
    int d_sample_rate;
    int d_max_burst_len;
    uint64_t d_abs_fft_index;
    uint64_t d_burst_id;
    uint64_t d_pre_burst_id;
    uint64_t d_n_tagged_bursts;
    uint64_t d_current_peaks;
    uint64_t extra;
    uint64_t d_rel_mag_hist;
    uint64_t d_rel_hist_index;
    uint32_t d_work_history_nffts;
    uint32_t d_work_sample_offset;

    float d_bin_width_db;

    float* d_window_f;
    float* d_magnitude_shifted_f;
    float* d_fine_window_f;
    float* d_fine_magnitude_shifted_f;
    float* d_baseline_sum_f;
    float* d_baseline_history_f;
    float* d_relative_magnitude_f;
    float* d_relative_history_f;
    uint32_t* d_burst_mask_i;
    uint32_t* d_burst_mask_j;  // 1 FFT buffer to prevent burst energy from impacting noise
    float* d_ones_f;
    float d_threshold;
    float d_threshold_low;
    float d_center_freq;
    float d_filter_bandwidth;

    FILE* d_burst_debug_file;

    gr::fft::fft_complex* d_fft;
    gr::fft::fft_complex* d_fine_fft;
#ifdef __USE_MKL__
    DFTI_DESCRIPTOR_HANDLE m_fft;
    DFTI_DESCRIPTOR_HANDLE m_fine_fft;
#endif
    std::vector<peak> d_peaks;
    std::list<pre_burst> d_pre_bursts;
    std::list<burst> d_bursts;
    std::list<burst> d_new_bursts;
    std::list<burst> d_gone_bursts;
    std::vector<owners> d_mask_owners;
    std::vector<moving_average> d_bin_averages;

    bool compute_relative_magnitude(void);
    void update_circular_buffer(void);
    void extract_peaks(void);
    void save_peaks_to_debug_file(char* filename);
    void remove_currently_tracked_bursts(void);
    void update_active_bursts(void);
    void update_potential_bursts(void);
    void delete_gone_bursts(void);
    void create_new_bursts(const gr_complex* input, int fft);
    void create_new_potential_bursts(void);
    void tag_new_bursts(void);
    void tag_gone_bursts(int noutput_items);

    void add_ownership(const pre_burst& b);
    void update_ownership(const pre_burst& pb, const burst& b);
    void remove_ownership(const pre_burst& b);
    void remove_ownership(const burst& b);

    bool check_prev_magnitude(size_t bin);
    void publish_debug(void);

    std::mutex d_work_mutex;

    void _reset();

    timer d_fft_timer;
    timer d_update_pb_timer;
    timer d_update_ab_timer;
    timer d_remove_tb_timer;
    timer d_extract_timer;
    timer d_delete_timer;
    timer d_new_pb_timer;
    timer d_new_b_timer;
    timer d_update_cb_timer;
    timer d_total_timer;
    timer d_rel_mag_timer;
    timer d_other;

public:
    /**
     * Constructor
     *
     * @param center_freq - center frequency of data stream, unit Hz
     * @param fft_size - number of bins in the primary FFT
     * @param sample_rate - sample rate of incoming stream
     * @param burst_pre_len - number of FFTs before the burst to place the tag
     * @param burst_post_len - number of FFTs after the burst to place the tag
     * @param burst_width - estimated bandwidth of the burst in Hz
     * @param max_bursts - maximum number of bursts simultaneously detected
     * @param max_burst_len - bursts exceeding this length wwill be tagged immediately
     * @param threshold - threshold above dynamic noise average (dB)
     * @param history_size - number of FFTs to compute noise estimate over
     * @param lookahead - number of FFTs a preburst must be present to convert to a burst
     * @param debug - true enables debug functionality
     *
     */
    fft_burst_tagger_impl(float center_freq,
                          int fft_size,
                          int sample_rate,
                          int burst_pre_len,
                          int burst_post_len,
                          int burst_width,
                          int max_bursts,
                          int max_burst_len,
                          float threshold,
                          int history_size,
                          int lookahead,
                          bool debug);

    /**
     * Deconstructor
     */
    ~fft_burst_tagger_impl();


    /**
     * Returns total number of bursts seen
     *
     * @return uint64_t - number of bursts
     */
    uint64_t get_n_tagged_bursts();

    /**
     * Resets burst tagger
     */
    void reset();

    bool stop();

    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items);

    /**
     * Sets max burst bandwidth
     * Unit: Hz
     *
     * @param bw - bandwidth
     */
    void set_max_burst_bandwidth(double bw) { d_filter_bandwidth = bw; }
};

} // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_IMPL_H */
