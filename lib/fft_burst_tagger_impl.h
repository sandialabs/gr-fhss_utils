/* -*- c++ -*- */
/*
 * Copyright 2018 <+YOU OR YOUR COMPANY+>.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_IMPL_H
#define INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_IMPL_H

#include <fhss_utils/fft_burst_tagger.h>
#include <gnuradio/fft/fft.h>
#include <mutex>
//#define __USE_MKL__
#ifdef __USE_MKL__
#include "mkl_dfti.h"
#endif

// Used for profiling timer
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

    class timer {
      public:
        timer() : total(0) {}
        // The timers can take a lot of extra cycles.  Only compile them in if wanted.
        #ifndef DO_TIMER
        void start() {}
        void end() {}
        #else
        void start() {start_time = high_resolution_clock::now();}
        void end() {
          high_resolution_clock::time_point end_time = high_resolution_clock::now();
          duration<double> span = duration_cast<duration<double> >(end_time - start_time);
          total += span.count();
        }
        #endif
        double elapsed() { return total;}
        void reset() { total = 0; }
      private:
        double total;
        high_resolution_clock::time_point start_time;
    };

    class movingAverage {
      // Calculate a moving average of floats using a circular buffer.
      public:
        movingAverage(size_t window) : N(window) {
          sum = curIndex = 0;
          hist.resize(N);
          memset(&hist[0], 0, sizeof(float)*N);
        }
        float add(float p) {
          sum += p - hist[curIndex];
          hist[curIndex++] = p;
          if (curIndex == N) curIndex = 0;
          return sum;
        }
      private:
        std::vector<float> hist;
        size_t curIndex;
        size_t N;
        float sum;
    };

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
        bool valid;
    };

    struct owners {
        owners() { clear(); }
        uint64_t ids[4];
        uint64_t _size;
        uint64_t uid;

        size_t size() { return _size;}

        void clear() {
          ids[0] = ids[1] = ids[2] = ids[3] = (uint64_t)-1;
          _size = 0;
        }

        void push_back(uint64_t id) {
          if (_size > 3) {
            printf("Owners::push_back - Trying to add too many points for  bin id=%zu, size=%zu, burst=%zu\n", uid, _size, id);
            printf("Owner bins are %zu, %zu, %zu %zu\n", ids[0], ids[1], ids[2], ids[3]);
            return;
          }
          ids[_size] = id;
          _size++;
        }

        void update(uint64_t oldID, uint64_t newID) {
          if (ids[0] == oldID) ids[0] = newID;
          else if (ids[1] == oldID) ids[1] = newID;
          else if (ids[2] == oldID) ids[2] = newID;
          else if (ids[3] == oldID) ids[3] = newID;
          else {
            printf("Owners::Update - Couldn't find id to update. This should never happen\n");
          }
        }

        void erase(uint64_t id) {
          _size--;
          if (ids[0] == id) {
            if (_size) {
              ids[0] = ids[_size];
              ids[_size] = (uint64_t)-1;
            } else ids[0] = (uint64_t)-1;
          } else if (ids[1] == id) {
            if (_size > 1) {
              ids[1] = ids[_size];
              ids[_size] = (uint64_t)-1;
            } else ids[1] = (uint64_t)-1;
          } else if (ids[2] == id) {
            if (_size > 2) {
              ids[2] = ids[_size];
              ids[_size] = (uint64_t)-1;
            } else ids[2] = (uint64_t)-1;
          } else if (ids[3] == id) {
            ids[3] = (uint64_t)-1;
          } else {
            // Increment size because we didn't remove anything.
            _size++;
            printf("Owners::Remove - Couldn't find id to erase. This should never happen - bin id = %zu, burst=%zu\n", uid, id);
          }
        }
    };

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

      float * d_window_f;
      float * d_magnitude_shifted_f;
      float * d_fine_window_f;
      float * d_fine_magnitude_shifted_f;
      float * d_baseline_sum_f;
      float * d_baseline_history_f;
      float * d_relative_magnitude_f;
      float * d_relative_history_f;
      uint32_t * d_burst_mask_i;
      float * d_ones_f;
      float d_threshold;
      float d_threshold_low;
      float d_center_frequency;
      float d_filter_bandwidth;

      FILE * d_burst_debug_file;

      gr::fft::fft_complex *d_fft;
      gr::fft::fft_complex *d_fine_fft;
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
      std::vector<movingAverage> d_bin_averages;

      bool compute_relative_magnitude(void);
      void update_circular_buffer(void);
      void extract_peaks(void);
      void save_peaks_to_debug_file(char * filename);
      void remove_currently_tracked_bursts(void);
      void update_active_bursts(void);
      void update_potential_bursts(void);
      void delete_gone_bursts(void);
      void create_new_bursts(const gr_complex* input);
      void create_new_potential_bursts(void);
      void tag_new_bursts(void);
      void tag_gone_bursts(int noutput_items);

      void add_ownership(const pre_burst& b);
      void update_ownership(const pre_burst& pb, const burst& b);
      void remove_ownership(const pre_burst& b);
      void remove_ownership(const burst& b);

      bool check_prev_magnitude(size_t bin);

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

      // PMT definitions
      const pmt::pmt_t PMT_RX_RATE = pmt::mp("rx_rate");
      const pmt::pmt_t PMT_RX_FREQ = pmt::mp("rx_freq");
      const pmt::pmt_t PMT_BURST_ID = pmt::mp("burst_id");
      const pmt::pmt_t PMT_REL_FREQ = pmt::mp("relative_frequency");
      const pmt::pmt_t PMT_CENTER_FREQ = pmt::mp("center_frequency");
      const pmt::pmt_t PMT_MAG = pmt::mp("magnitude");
      const pmt::pmt_t PMT_SAMPLE_RATE = pmt::mp("sample_rate");
      const pmt::pmt_t PMT_NOISE_POWER = pmt::mp("noise_power");
      const pmt::pmt_t PMT_BANDWIDTH = pmt::mp("bandwidth");

     public:
      fft_burst_tagger_impl(float center_frequency, int fft_size, int sample_rate, int burst_pre_len, int burst_post_len, int burst_width, int max_bursts, int max_burst_len, float threshold, int history_size, int lookahead, bool debug);
      ~fft_burst_tagger_impl();

      uint64_t get_n_tagged_bursts();

      void reset();
      bool stop();

      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);

       const pmt::pmt_t NEW_BURST_TAG = pmt::mp("new_burst");
       const pmt::pmt_t GONE_BURST_TAG = pmt::mp("gone_burst");
       const pmt::pmt_t RX_TIME_TAG = pmt::mp("rx_time");
       const pmt::pmt_t META_RELATIVE_FREQ = pmt::mp("relative_frequency");
       const pmt::pmt_t META_ID = pmt::mp("burst_id");
       const pmt::pmt_t META_MAGNITUDE = pmt::mp("magnitude");
       const pmt::pmt_t META_CENTER_FREQ = pmt::mp("center_frequency");
       const pmt::pmt_t META_SAMP_RATE = pmt::mp("sample_rate");
       const pmt::pmt_t META_START_TIME = pmt::mp("start_time");
       const pmt::pmt_t META_DURATION = pmt::mp("duration");
       const pmt::pmt_t META_NOISE_POWER = pmt::mp("noise_power");
       const pmt::pmt_t META_BANDWIDTH = pmt::mp("bandwidth");

       const pmt::pmt_t PDU_OUT_PORT = pmt::mp("cpdus");

       void set_max_burst_bandwidth(double bw) { d_filter_bandwidth = bw; }
    };

  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_IMPL_H */
