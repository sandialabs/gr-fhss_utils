/* -*- c++ -*- */
/*
 * Copyright 2018, 2019, 2020 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_FHSS_UTILS_TAGGED_BURST_TO_PDU_IMPL_H
#define INCLUDED_FHSS_UTILS_TAGGED_BURST_TO_PDU_IMPL_H

#include <gnuradio/blocks/rotator.h>
#include <gnuradio/filter/fir_filter.h>
#include <fhss_utils/tagged_burst_to_pdu.h>
#include <fhss_utils/constants.h>
#include <boost/lockfree/queue.hpp>
#include <queue>

namespace gr {
namespace fhss_utils {

struct burst_data {
    uint64_t id;
    uint64_t offset;
    uint64_t data_skip;
    float magnitude;
    float relative_frequency;
    float center_frequency;
    float sample_rate;
    size_t len;
    size_t rot_skip;
    pmt::pmt_t dict;
    gr_complex* data;
    gr_complex* rot_tmp;
    blocks::rotator rotate;
};

struct two_gr_complex {
    two_gr_complex() {}
    two_gr_complex(gr_complex* a, gr_complex* b) : one(a), two(b) {}
    gr_complex* one;
    gr_complex* two;
};

struct buffer {
    buffer(size_t capacity = 0)
    {
        reset();
        if (capacity > 0)
            data.resize(int(capacity * 1.1));
        data.resize(0);
        desired_max = capacity;
    }
    std::vector<tag_t> new_burst_tags;
    std::vector<tag_t> rx_time_tags;
    std::vector<tag_t> gone_burst_tags;
    std::vector<gr_complex> data;
    size_t desired_max;
    bool end_flag;
    size_t start;

    bool is_full() { return data.size() >= desired_max; }

    void reset()
    {
        new_burst_tags.resize(0);
        rx_time_tags.resize(0);
        gone_burst_tags.resize(0);
        data.resize(0);
        end_flag = false;
    }

    void add_data(const gr_complex* d, size_t length)
    {
        // size_t clength = std::min(length, data.capacity() - data.size());
        // printf("Going to add %zu data of %zu, cap = %zu, size = %zu, ptr = %p\n",
        // clength, length, data.capacity(), data.size(), this);
        data.insert(data.end(), d, d + length);
    }

    void add_tags(const std::vector<tag_t> &new_b,
                  const std::vector<tag_t> &rx_time,
                  const std::vector<tag_t> &gone_b)
    {
        new_burst_tags.insert(new_burst_tags.end(), new_b.begin(), new_b.end());
        rx_time_tags.insert(rx_time_tags.end(), rx_time.begin(), rx_time.end());
        gone_burst_tags.insert(gone_burst_tags.end(), gone_b.begin(), gone_b.end());
    }
};


/**
 * Implementing class of tagged_burst_to_pdu
 */
class tagged_burst_to_pdu_impl : public tagged_burst_to_pdu
{
private:
    bool d_debug;
    float d_relative_center_frequency;
    float d_relative_span;
    float d_relative_sample_rate;
    float d_sample_rate;
    int d_min_burst_size;
    int d_max_burst_size;
    int d_outstanding;
    int d_max_outstanding;
    uint64_t d_n_dropped_bursts;
    size_t d_block_increment;
    bool d_blocked;
    size_t d_max_id;
    int d_num_threads;

    buffer* d_current_buffer;

    tag_t d_current_rx_time_tag;
    std::queue<tag_t> d_rx_time_tags;
    std::queue<two_gr_complex> d_alloced_arrays;
    std::vector<float> d_taps;
    size_t d_decimation;
    std::vector<filter::kernel::fir_filter_ccf*> d_input_fir_filters;

    boost::lockfree::queue<buffer*> d_write_queue;
    boost::lockfree::queue<buffer*> d_work_queue;
    boost::thread* process_thread;

    float d_lower_border;
    float d_upper_border;

    std::map<uint64_t, burst_data> d_bursts;
    std::map<uint64_t, burst_data> d_new_bursts;

    void append_to_burst(burst_data& burst, const gr_complex* data, size_t n);
    void publish_burst(burst_data& burst);

    void create_new_bursts(const buffer& work_buffer);
    void publish_and_remove_old_bursts(const buffer& work_buffer);
    void update_current_bursts(int noutput_items, const gr_complex* in);
    void process_data();

    int get_output_queue_size();
    int get_output_max_queue_size();
    void burst_handled(pmt::pmt_t msg);

    double convert_rx_time(const tag_t& rx_time_tag);

    // The gnuradio filter function reads memory that it shouldn't.  We offset our array
    // by this amount to prevent bad things from happening.  If we fix the filter
    // function, then we can get rid of this.
    const size_t ROT_TMP_OFFSET = 4;

public:
    // d_block_size is used as the buffer size and should always be set to
    // far more than half of the tap length, since that is the overlap region
    static const size_t d_block_size = 32 * 1024;
    // d_num_buffers was empirically chosen to be 5. If we need more than 5,
    // the block is probably backed up significantly and should back-pressure
    // the gnuradio scheduler
    static const size_t d_num_buffers = 5;

    /**
     * Constructor
     *
     * @param decimation -
     * @param taps -
     * @param min_burst_time -
     * @param max_burst_time -
     * @param relative_center_frequency -
     * @param relative_span -
     * @param relative_sample_rate -
     * @param sample_rate -
     * @param num_threads -
     */
    tagged_burst_to_pdu_impl(size_t decimation,
                             const std::vector<float>& taps,
                             float min_burst_time,
                             float max_burst_time,
                             float relative_center_frequency,
                             float relative_span,
                             float relative_sample_rate,
                             float sample_rate,
                             int num_threads);

    /**
     * Deconstructor
     */
    ~tagged_burst_to_pdu_impl();

    bool stop();

    /*uint64_t get_n_dropped_bursts();*/

    /**
     * Processes input streadm
     *
     * @param noutput_items -
     * @param input_items -
     * @param output_items -
     * @return int -
     */
    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items);
};

} // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_TAGGED_BURST_TO_PDU_IMPL_H */
