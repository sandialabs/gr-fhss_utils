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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>

#include "fft_burst_tagger_impl.h"

#include <volk/volk.h>

#include <stdio.h>
#include <inttypes.h>

#ifdef __AVX__
#include <immintrin.h>
#endif

namespace gr {
  namespace fhss_utils {

    fft_burst_tagger::sptr
    fft_burst_tagger::make(float center_frequency, int fft_size, int sample_rate, int burst_pre_len, int burst_post_len, int burst_width, int max_bursts, int max_burst_len, float threshold, int history_size, int lookahead, bool debug)
    {
      return gnuradio::get_initial_sptr
        (new fft_burst_tagger_impl(center_frequency, fft_size, sample_rate, burst_pre_len, burst_post_len, burst_width, max_bursts, max_burst_len, threshold, history_size, lookahead, debug));
    }

    /*
     * The private constructor
     */
    fft_burst_tagger_impl::fft_burst_tagger_impl(float center_frequency, int fft_size, int sample_rate, int burst_pre_len, int burst_post_len, int burst_width, int max_bursts, int max_burst_len, float threshold, int history_size, int lookahead, bool debug)
      : gr::sync_block("fft_burst_tagger",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_center_frequency(center_frequency), d_sample_rate(sample_rate),
        d_fft_size(fft_size), d_burst_pre_len(burst_pre_len),
        d_lookahead(lookahead),
        d_burst_id(0),
        d_pre_burst_id(0),
        d_n_tagged_bursts(0),
        d_abs_fft_index(burst_pre_len-1),
        d_max_burst_len(max_burst_len),
        d_fft(NULL), d_history_size(history_size),
        d_history_primed(false), d_history_index(0),
        d_burst_post_len(burst_post_len), d_debug(debug), d_burst_debug_file(NULL),
        d_mask_owners(d_fft_size)

    {
        const int nthreads = 1;
        d_fine_fft_size = d_fft_size * std::min(16, (int)pow(2, int(log(d_lookahead) / log(2))));
        d_fft = new fft::fft_complex(d_fft_size, true, nthreads);
        d_fine_fft = new fft::fft_complex(d_fine_fft_size, true, nthreads);

        #ifdef __USE_MKL__
        MKL_LONG status = DftiCreateDescriptor( &m_fft, DFTI_SINGLE, DFTI_COMPLEX, 1, d_fft_size);
        status = DftiSetValue( m_fft, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
        status = DftiCommitDescriptor(m_fft);
        status = DftiCreateDescriptor( &m_fine_fft, DFTI_SINGLE, DFTI_COMPLEX, 1, d_fine_fft_size);
        status = DftiSetValue( m_fine_fft, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
        status = DftiCommitDescriptor(m_fine_fft);
        #endif

        set_output_multiple(d_fft_size);

        // We need to keep d_burst_pre_len samples
        // in the buffer to be able to tag a burst at it's start.
        // Set the history to this + 1, so we always have
        // this amount of samples available at the start of
        // our input buffer.
        set_history(d_burst_pre_len * d_fft_size + 1);
        // This makes sure we have at least d_lookahead FFTs in the
        // buffer, not including history. There will always be enough
        // data to convert a pre_burst to a burst, and still tag it
        // d_burst_pre_len FFTs in the past.
        // TODO: this call doesn't do anything? https://github.com/gnuradio/gnuradio/issues/1483
        set_min_noutput_items((d_lookahead + 1) * d_fft_size);

        d_window_f = (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
        std::vector<float> window = fft::window::build(fft::window::WIN_BLACKMAN, d_fft_size, 0);
        memcpy(d_window_f, &window[0], sizeof(float) * d_fft_size);

        d_fine_window_f = (float *)volk_malloc(sizeof(float) * d_fine_fft_size, volk_get_alignment());
        std::vector<float> fine_window = fft::window::build(fft::window::WIN_BLACKMAN, d_fine_fft_size, 0);
        memcpy(d_fine_window_f, &fine_window[0], sizeof(float) * d_fine_fft_size);

        d_baseline_history_f = (float *)volk_malloc(sizeof(float) * d_fft_size * d_history_size, volk_get_alignment());
        d_baseline_sum_f = (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
        d_magnitude_shifted_f = (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
        d_fine_magnitude_shifted_f = (float *)volk_malloc(sizeof(float) * d_fine_fft_size, volk_get_alignment());
        d_relative_magnitude_f = (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
        d_burst_mask_i = (uint32_t *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());

        d_rel_mag_hist = 1;
        d_rel_hist_index = 0;
        d_relative_history_f = (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());

        memset(d_baseline_history_f, 0, sizeof(float) * d_fft_size * d_history_size);
        memset(d_baseline_sum_f, 0, sizeof(float) * d_fft_size);
        memset(d_magnitude_shifted_f, 0, sizeof(float) * d_fft_size);
        memset(d_fine_magnitude_shifted_f, 0, sizeof(float) * d_fft_size);
        memset(d_relative_magnitude_f, 0, sizeof(float) * d_fft_size);
        extra = 0;

        memset(d_burst_mask_i, ~0, sizeof(uint32_t)*d_fft_size);

        d_threshold = pow(10, threshold/10) / d_history_size;
        //d_threshold_low = d_threshold / 2.0;  // hysteresis
        d_threshold_low = d_threshold;
        if(d_debug) {
          fprintf(stderr, "threshold=%f, d_threshold=%f (%f/%d)\n",
              threshold, d_threshold, d_threshold * d_history_size, d_history_size);
        }

        d_peaks.resize(d_fft_size);
        d_current_peaks = 0;
        d_bin_averages.resize(d_fft_size, movingAverage(d_history_size));

        if(max_bursts){
          d_max_bursts = max_bursts;
        } else {
          // Consider the signal to be invalid if more than 80%
          // of all channels are in use.
          d_max_bursts = (sample_rate / burst_width) * 0.8;
        }

        for (size_t i = 0; i < d_fft_size; i++) d_mask_owners[i].uid = i;

        // Area to ignore around an already found signal in FFT bins
        // Internal representation is in FFT bins
        d_burst_width = burst_width / (sample_rate / fft_size);

        d_filter_bandwidth = 0;

        if(d_debug) {
          fprintf(stderr, "d_max_bursts=%d\n", d_max_bursts);
        }

        if(d_debug) {
          d_burst_debug_file = fopen("/tmp/fft_burst_tagger-bursts.log", "w");
        }
    }

    /*
     * Our virtual destructor.
     */
    fft_burst_tagger_impl::~fft_burst_tagger_impl()
    {
      fprintf(stderr, "Tagged %" PRIu64 " bursts\n", d_n_tagged_bursts);
      delete d_fft;
      delete d_fine_fft;
      volk_free(d_window_f);
      volk_free(d_baseline_history_f);
      volk_free(d_baseline_sum_f);
      volk_free(d_relative_magnitude_f);
      volk_free(d_relative_history_f);
      volk_free(d_magnitude_shifted_f);
      volk_free(d_fine_window_f);
      volk_free(d_fine_magnitude_shifted_f);
      volk_free(d_burst_mask_i);
      if(d_burst_debug_file) {
        fclose(d_burst_debug_file);
      }
      #ifdef __USE_MKL__
      MKL_LONG status = DftiFreeDescriptor(&m_fft);
      status = DftiFreeDescriptor(&m_fine_fft);
      #endif
    }

    bool fft_burst_tagger_impl::stop() {
      #ifdef DO_TIMER
      printf("total time: %f\n", d_total_timer.elapsed());
      printf("rel mag: %f\n", d_rel_mag_timer.elapsed());
      printf("fft time: %f\n", d_fft_timer.elapsed());
      printf("update pb time: %f\n",d_update_pb_timer.elapsed());
      printf("update ab time: %f\n",d_update_ab_timer.elapsed());
      printf("remove tb time: %f\n",d_remove_tb_timer.elapsed());
      printf("extract time: %f\n",d_extract_timer.elapsed());
      printf("delete time: %f\n",d_delete_timer.elapsed());
      printf("new pb time: %f\n",d_new_pb_timer.elapsed());
      printf("new burst time: %f\n",d_new_b_timer.elapsed());
      printf("update cb time: %f\n",d_update_cb_timer.elapsed());
      printf("other timer: %f\n", d_other.elapsed());
      #endif
      printf("extra = %zu\n", extra);
    }

    /*
     *  Compute the magnitude of current data normalized by nearby magnitudes
     */
    bool
    fft_burst_tagger_impl::compute_relative_magnitude(void)
    {
      // Only compute if we have sufficient data in our circular buffer
      d_rel_mag_timer.start();
      if(!d_history_primed) {
        return false;
      }

      volk_32f_x2_divide_32f(d_relative_magnitude_f, d_magnitude_shifted_f, d_baseline_sum_f, d_fft_size);
      d_rel_mag_timer.end();
      return true;
    }

#define HIST(i) (d_baseline_history_f + (i % d_history_size) * d_fft_size)
    void
    fft_burst_tagger_impl::update_circular_buffer(void)
    {
      // We only update the average if there is no burst going on at the moment
      // TODO: also block if there are any pre_bursts. Is that what we want?
      if (!d_history_primed) {
        for (size_t i = 0; i < d_fft_size; i++) {
          d_baseline_sum_f[i] = d_bin_averages[i].add(d_magnitude_shifted_f[i]);
        }
        d_history_index++;
        //printf("updated baseline %u\n", d_history_index);

        if(d_history_index == d_history_size) {
          d_history_primed = true;
        }
      } else if ((d_abs_fft_index & 0x3) == 0) {
        for (size_t i = 0; i < d_fft_size; i++) {
          if (d_burst_mask_i[i] != 0) {
            d_baseline_sum_f[i] = d_bin_averages[i].add(d_magnitude_shifted_f[i]);
          }
        }
      }
      // Copy the relative magnitude history into circular buffer
      memcpy(d_relative_history_f + d_rel_hist_index * d_fft_size, d_relative_magnitude_f, sizeof(float) *d_fft_size);
      if (d_rel_mag_hist > 0)
        d_rel_hist_index = (d_rel_hist_index + 1) % d_rel_mag_hist;
    }

    void
    fft_burst_tagger_impl::update_active_bursts(void)
    {
      auto b = std::begin(d_bursts);
      while (b != std::end(d_bursts)) {
        if(d_relative_magnitude_f[b->center_bin-1] > d_threshold_low ||
           d_relative_magnitude_f[b->center_bin] > d_threshold_low ||
           d_relative_magnitude_f[b->center_bin+1] > d_threshold_low) {
          b->last_active = d_abs_fft_index;
        }
        ++b;
      }
    }

    void
    fft_burst_tagger_impl::reset() {
      std::lock_guard<std::mutex> lock(d_work_mutex);
      _reset();
    }

    void
    fft_burst_tagger_impl::_reset() {

      // close and tag all current bursts
      //d_peaks.clear();
      d_current_peaks = 0;
      d_pre_bursts.clear();
      d_new_bursts.clear();
      for(burst b : d_bursts) {
        b.stop = d_abs_fft_index;
        if (b.valid) d_gone_bursts.push_back(b);
      }
      d_bursts.clear();
      for (size_t i = 0; i < d_fft_size; i++)
        d_mask_owners[i].clear();

      tag_gone_bursts(d_abs_fft_index*d_fft_size);

      // reset baseline
      d_history_index = 0;
      d_history_primed = false;
      memset(d_baseline_sum_f, 0, sizeof(float)*d_fft_size);
      memset(d_baseline_history_f, 0, sizeof(float)*d_fft_size*d_history_size);
      memset(d_burst_mask_i, ~0, sizeof(uint32_t)*d_fft_size);
    }

    void
    fft_burst_tagger_impl::delete_gone_bursts(void)
    {
      auto b = std::begin(d_bursts);

      while (b != std::end(d_bursts)) {
        if((b->last_active + d_burst_post_len) < d_abs_fft_index
            || (d_max_burst_len && d_abs_fft_index - b->start > d_max_burst_len)) {
          //printf("Deleting gone burst %" PRIu64 " (start=%" PRIu64 ", d_abs_fft_index=%" PRIu64 ")\n", b->id, b->start, d_abs_fft_index);
          b->stop = d_abs_fft_index;
          if (b->valid) d_gone_bursts.push_back(*b);
          remove_ownership(*b);
          b = d_bursts.erase(b);
        } else {
          ++b;
        }
      }
    }

    bool fft_burst_tagger_impl::check_prev_magnitude(size_t bin) {
      // This is a cicular buffer, order doesn't really matter though.
      size_t cIndex = 0;
      size_t start_bin = std::max(bin - d_burst_width / 2, (size_t)0);
      size_t stop_bin = std::min(bin + d_burst_width / 2, (size_t)(d_fft_size-1));
      for (size_t i = 0; i < d_rel_mag_hist; i++) {
        bool found = false;
        for (size_t b = start_bin; b <= stop_bin; b++) {
          if (d_relative_history_f[cIndex+b] > d_threshold) {
            found = true;
          }
        }
        cIndex += d_fft_size;
        if (found == false) return false;
      }
      return true;
    }

    void
    fft_burst_tagger_impl::create_new_potential_bursts(void)
    {
      // A potential burst is one that we are not yet sure about.  The burst needs to exist for a sufficent amount of time
      // in order for us to keep it.
      bool allow[d_fft_size];
      memset(allow, 1, sizeof(bool)*d_fft_size);
      for (size_t i = 0; i < d_current_peaks; i++) {
        peak& p = d_peaks[i];
        if(d_burst_mask_i[p.bin] && allow[p.bin]) {
          if (check_prev_magnitude(p.bin) ) {
            pre_burst b;
            b.center_bin = p.bin;
            b.start_bin = std::max(b.center_bin - d_burst_width / 2, 0);
            b.stop_bin = std::min(b.center_bin + d_burst_width / 2, d_fft_size-1);

            b.peak_count = 1 + d_rel_mag_hist;
            b.max_relative_magnitude = p.relative_magnitude;
            size_t start_bin = b.start_bin;
            size_t stop_bin = b.stop_bin;

            // find the largest peak within d_burst_width/2 bins on either side of the last detection
            for (size_t j = 0; j < d_rel_mag_hist; j++) {
              for(int k = start_bin; k < stop_bin; k++) {
                if (d_relative_history_f[k] > b.max_relative_magnitude) {
                  b.max_relative_magnitude = d_relative_history_f[k];
                  b.center_bin = k;
                }
              }
            }


            // The burst might have started one FFT earlier
            b.start = d_abs_fft_index - d_burst_pre_len - d_rel_mag_hist;
            b.id = d_pre_burst_id++;
            b.thresh_count = 1 + d_rel_mag_hist;
            //printf("creating new pre_burst (%lu) with offset %lu, at bin %u\n", b.id, b.start*d_fft_size, b.center_bin);

            d_pre_bursts.push_back(b);
            add_ownership(b);

            if(d_burst_debug_file) {
              fprintf(d_burst_debug_file, "%" PRIu64 ",%d,x\n", b.start, b.center_bin);
              //float f_rel = (b.center_bin - d_fft_size / 2) / float(d_fft_size);
              //fprintf(d_burst_debug_file, "%f,%f,x\n", b.start/4e6, f_rel * 4e6 + 1624800000);
            }
          } else {
            size_t start_bin = std::max((int)p.bin - d_burst_width / 2, 0);
            size_t stop_bin = std::min((int)p.bin + d_burst_width / 2, (d_fft_size-1));
            for (size_t j = start_bin; j <= stop_bin; j++) {
              allow[j] = 0;
            }
            //printf("Rejecting pre_burst with offset %lu, at bin %u\n",  (d_abs_fft_index - d_burst_pre_len - d_rel_mag_hist)*d_fft_size, p.bin);
          }
        }
      }
      // TODO: move this
      if(d_max_bursts > 0 && d_bursts.size() > d_max_bursts) {
        printf("=========== RESETTING BURST DETECTOR ===========\n");
        _reset();
        /*fprintf(stderr, "Detector in burst squelch at %f\n", d_abs_fft_index * d_fft_size / float(d_sample_rate));
        d_new_bursts.clear();
        for(burst b : d_bursts) {
          if(b.start != d_abs_fft_index - d_burst_pre_len) {
            b.stop = d_abs_fft_index;
            d_gone_bursts.push_back(b);
          }
        }
        d_bursts.clear();*/
      }
    }

    void
    fft_burst_tagger_impl::add_ownership(const pre_burst& b) {
      //printf("%zu: Adding pre burst %zu: start = %d, stop = %d \n", d_abs_fft_index, b.id, b.start_bin, b.stop_bin);

      for (size_t i = b.start_bin; i <= b.stop_bin; i++) {
        d_burst_mask_i[i] = 0;
        d_mask_owners[i].push_back(b.id);
      }
      //printf("Yo: %zu, %zu\n", d_mask_owners[151].ids[0], d_mask_owners[151].ids[1]);
    }

    void
    fft_burst_tagger_impl::update_ownership(const pre_burst& pb, const burst& b) {
      //printf("%zu: Updating pre burst %zu: start = %d, stop = %d, burst %zu: start = %d, stop = %d\n",  d_abs_fft_index, pb.id, pb.start_bin, pb.stop_bin, b.id, b.start_bin, b.stop_bin);
      if (pb.start_bin != b.start_bin || pb.stop_bin != b.stop_bin) {
        // There is a chance that we will have two pre-bursts that
        for (size_t i = pb.start_bin; i <= pb.stop_bin; i++) {
          if (d_mask_owners[i].size() == 1) {
            d_burst_mask_i[i] = ~0;
          }
          d_mask_owners[i].erase(pb.id);
        }
        for (size_t i = b.start_bin; i <= b.stop_bin; i++) {
          d_burst_mask_i[i] = 0;
          d_mask_owners[i].push_back(b.id);
        }
      } else {
        for (size_t i = b.start_bin; i <= b.stop_bin; i++) {
          d_mask_owners[i].update(pb.id, b.id);
        }
      }
      //printf("Yo: %zu, %zu %d %d\n", d_mask_owners[151].ids[0], d_mask_owners[151].ids[1], pb.start_bin, b.start_bin);
    }

    void
    fft_burst_tagger_impl::remove_ownership(const pre_burst& b) {
      //printf("%zu: Removing pre burst %zu: start = %d, stop = %d\n", d_abs_fft_index, b.id, b.start_bin, b.stop_bin);
      //printf("Yo: %zu, %zu\n", d_mask_owners[151].ids[0], d_mask_owners[151].ids[1]);
      for (size_t i = b.start_bin; i <= b.stop_bin; i++) {
        if (d_mask_owners[i].size() == 1) {
          d_burst_mask_i[i] = ~0;
        }
        d_mask_owners[i].erase(b.id);
      }

    }

    void
    fft_burst_tagger_impl::remove_ownership(const burst& b) {
      //printf("%zu: Removing burst %zu: start = %d, stop = %d\n", d_abs_fft_index, b.id, b.start_bin, b.stop_bin);
      for (size_t i = b.start_bin; i <= b.stop_bin; i++) {
        if (d_mask_owners[i].size() == 1) {
          d_burst_mask_i[i] = ~0;
        }
        d_mask_owners[i].erase(b.id);
      }

    }

    void
    fft_burst_tagger_impl::create_new_bursts(const gr_complex* input)
    {
      // Convert potential bursts into burst once they have existed for a sufficient amount of time.
      auto b = d_pre_bursts.begin();
      while (b != d_pre_bursts.end()) {

        // if d_lookahead peaks have been measured, convert the pre_burst into a burst
        if (b->peak_count >= d_lookahead + 1) {
          burst new_b;
          new_b.start = b->start;
          new_b.last_active = d_abs_fft_index;
          // use the max peak's bin and magnitude for the burst
          // take the last N samples up to this point.
          size_t offset = d_abs_fft_index * d_fft_size - d_fine_fft_size - nitems_read(0);
          volk_32fc_32f_multiply_32fc(d_fine_fft->get_inbuf(), &input[offset], d_fine_window_f, d_fine_fft_size);
          #ifdef __USE_MKL2__
          DftiComputeForward(m_fine_fft, d_fine_fft->get_inbuf(), d_fine_fft->get_outbuf());
          #else
          d_fine_fft->execute();
          #endif

          // Get the mag squared of the fft and fft shift
          size_t fftd2 = d_fine_fft_size/2;
          volk_32fc_magnitude_squared_32f(d_fine_magnitude_shifted_f + fftd2, d_fine_fft->get_outbuf(), fftd2);
          volk_32fc_magnitude_squared_32f(d_fine_magnitude_shifted_f, d_fine_fft->get_outbuf() + fftd2, fftd2);

          // Search around the peak to try to estimate the bandwidth
          size_t factor = d_fine_fft_size / d_fft_size;
          // Ensure that a center bin - burst_width < 0 becomes a 1, not a very large positive number.
          size_t search_start = (size_t)std::max((int)factor*(b->center_bin - d_burst_width / 2), 1);
          size_t search_stop = std::min(factor*(b->center_bin + d_burst_width / 2), (size_t)(d_fine_fft_size-2));

          // Normalize the magnitude by the noise floor
          for (size_t i = search_start; i < search_stop; i++) {
            d_fine_magnitude_shifted_f[i] /= d_baseline_sum_f[i / factor] / d_lookahead;
          }

          float* max_entry = std::max_element(d_fine_magnitude_shifted_f + search_start, d_fine_magnitude_shifted_f + search_stop);
          // Smooth to get max value.  This is safe because of the bounds on search start and stop.
          float max_value = .5*max_entry[0] + .25*(max_entry[-1] + max_entry[1]);
          // How high is this above the noise floor?
          float snr_est = max_value;
          float *n = d_baseline_sum_f + b->center_bin;

          size_t minIndex = search_start;
          size_t maxIndex = search_stop;
          // 6 db from peak
          float thresh = max_value * .0625;
          //float thresh = (max_value - noise_floor) * 0.0625;
          while(minIndex <= search_stop && d_fine_magnitude_shifted_f[minIndex] < thresh) {
            minIndex++;
          }
          while(maxIndex >= search_start && d_fine_magnitude_shifted_f[maxIndex] < thresh) {
            maxIndex--;
          }
          new_b.bandwidth = (maxIndex - minIndex + 1) * d_sample_rate / d_fine_fft_size;
          new_b.center_freq = ((ssize_t)(maxIndex + minIndex) - (ssize_t)d_fine_fft_size) / 2 * d_sample_rate / d_fine_fft_size;

          new_b.center_bin = (maxIndex + minIndex) / 2 / factor;
          new_b.start_bin = std::max(new_b.center_bin - d_burst_width / 2, 0);
          new_b.stop_bin = std::min(new_b.center_bin + d_burst_width / 2, d_fft_size-1);
          float max_relative_magnitude = b->max_relative_magnitude;

          // now that we know the center bin, make sure this burst doesn't already exist
          // TODO: improve this, maybe use the mask
          bool already_exists = false;
          for (burst existing_b : d_bursts) {
              if (abs(new_b.center_bin - existing_b.center_bin) <= (d_burst_width+1)/2) {
                already_exists = true;
                break;
              }
          }
          if (already_exists) {
            remove_ownership(*b);
            b = d_pre_bursts.erase(b);
            continue;
          }
          new_b.id = d_burst_id++;
          //printf("%zu: center = %d\n", new_b.id, new_b.center_bin);
          //printf("bandwidth info 1: snr = %f, min = %f, max = %f\n", snr_est, (float)minIndex * d_sample_rate / d_fine_fft_size - d_sample_rate / 2,
          //  (float)maxIndex * d_sample_rate / d_fine_fft_size - d_sample_rate /2);
          // normalize the magnitude
          new_b.magnitude = 10 * log10(max_relative_magnitude * d_history_size);
          // Allow for us to filter out the data if needed
          if (d_filter_bandwidth > 0 && new_b.bandwidth > d_filter_bandwidth) {
            // keep tracking the burst but don't write out the tags so we don't tune/filter/decimate
            new_b.valid = false;
            //printf("*************** Dropping burst\n");
          } else {
            d_new_bursts.push_back(new_b);
            new_b.valid = true;
          }
          d_bursts.push_back(new_b);
          //printf("created new burst (%lu) at offset %lu, bin = %zu\n", new_b.id, new_b.start*d_fft_size, new_b.center_bin);
          //printf("center bin is %d (%d), magnitude is %f (%f)\n", new_b.center_bin, b->peaks.front().bin, max_peak->relative_magnitude, b->peaks.front().relative_magnitude);
          update_ownership(*b, new_b);
          b = d_pre_bursts.erase(b);
          continue;
        }
        ++b;
      }
    }

    void
    fft_burst_tagger_impl::update_potential_bursts(void)
    {
      float* c_magnitude = d_relative_magnitude_f;
      float threshold_low = d_threshold_low*.5;
      auto b = d_pre_bursts.begin();
      while (b != d_pre_bursts.end()) {

        // initialize a peak to store the largest new peak in
        peak p;
        size_t start_bin = b->start_bin;
        size_t stop_bin = b->stop_bin;
        p.bin = stop_bin;
        p.relative_magnitude = c_magnitude[stop_bin];

        // find the largest peak within d_burst_width/2 bins on either side of the first detection
        for(int i = start_bin; i < stop_bin; i++) {
          if (c_magnitude[i] > p.relative_magnitude) {
            p.relative_magnitude = c_magnitude[i];
            p.bin = i;
          }
        }

        // if the largest peak is below the threshold, drop the pre_burst
        if (p.relative_magnitude < threshold_low) {
          // erase
          //printf("dropped pre_burst %lu****************************************\n", b->id);
          // Think about not creating a potential burst until we see it a few times.  In my dataset, 60% of pb vanish after 1 frame, and 90% after 4 frames.
          // I could also change this to be a hash table, with potentially cheaper add/delete.
          if (b->peak_count == 4) extra++;
          //printf("Deleted short pre_burst id=%zu, bin=%d, sample=%zu\n", b->id, b->center_bin, (b->start + b->peak_count)*d_fft_size);
          remove_ownership(*b);
          b = d_pre_bursts.erase(b);
          continue;
        } else {
          if (p.relative_magnitude > d_threshold_low) {
            b->thresh_count++;
          } else if (b->thresh_count < b->peak_count * .75) {
            remove_ownership(*b);
            b = d_pre_bursts.erase(b);
            continue;
          }
        }

        // otherwise, add the peak to the list of measurements
        b->peak_count++;
        if (p.relative_magnitude > b->max_relative_magnitude) {
          b->max_relative_magnitude = p.relative_magnitude;
          b->center_bin = p.bin;
        }
        ++b;
      }
    }

    void
    fft_burst_tagger_impl::remove_currently_tracked_bursts(void)
    {
      // Remove peaks from consideration that we are already tracking.  This function is called after we have
      // updated bursts, so we only want to see new ones present.
      // Volk function requires ints as inputs, but we are really just zeroing out and bin that we are already tracking
      volk_32i_x2_and_32i((int32_t*)d_relative_magnitude_f, (int32_t *)d_relative_magnitude_f, (int32_t*)d_burst_mask_i, d_fft_size);
    }

    inline void swap_peaks(std::vector<peak>& p, size_t i, size_t j) {
      peak temp = p[i];
      p[i] = p[j];
      p[j] = temp;
    }

    void
    fft_burst_tagger_impl::extract_peaks(void)
    {
      const uint32_t start_bin = d_burst_width / 2;
      const uint32_t end_bin = d_fft_size - d_burst_width / 2;
      size_t length = end_bin - start_bin;
      size_t length8 = 0;
      size_t index = 0;

      #ifdef __AVX__
      length8 = length >> 3;
      float* c_magnitude = d_relative_magnitude_f + start_bin;
      // Copy d_threshold to every position in 256 bit array
      __m256 threshold = _mm256_set1_ps(d_threshold);

      for (size_t j = 0; j < length8; j++, c_magnitude+=8) {
        // Unalgined load.  About ~1% slower than aligned load, not worth changing
        __m256 magnitude = _mm256_loadu_ps(c_magnitude);
        // Sub data from threshold
        __m256 diff = _mm256_sub_ps(threshold, magnitude);
        // negative numbers have a sign bit of 1.  Pull the sign bit off all 8 values and store in an int.
        int which = _mm256_movemask_ps(diff);

        if (which) {
          size_t loc = start_bin + (j << 3);
          while(which) {
            if (which & 0x1) {
              d_peaks[index].bin = loc;
              d_peaks[index++].relative_magnitude = d_relative_magnitude_f[loc];
            }
            loc++;
            which >>= 1;
          }
        }
      }
      #endif
      // Do the trailer and get the non multiple of 8 tail.
      for (size_t bin = (length8 << 3) + start_bin; bin < end_bin; bin++) {
        float rel_mag = d_relative_magnitude_f[bin];
        if(rel_mag > d_threshold) {
          d_peaks[index].bin = bin;
          d_peaks[index++].relative_magnitude = rel_mag;
          //printf("ts %" PRIu64 " bin %d val %f\n", d_abs_fft_index, p.bin, p.relative_magnitude);
        }
      }
      // We only need to sort if the last time we had peaks.  Otherwise there can't be any that will pass this time.
      bool doSort = d_current_peaks > 0 && d_rel_mag_hist > 0;
      d_current_peaks = index;

      // sort by descending relative magnitude
      // Special case a few things here. To save time
      // Case 0,1: Just return
      // Case 2: swap if needed
      // Case 3 is a little more complicated, but not bad
      if (!doSort) return;
      if (index == 2) {
        if (d_peaks[0].relative_magnitude < d_peaks[1].relative_magnitude) {
          swap_peaks(d_peaks, 0, 1);
        }
      } else if (index == 3) {
        if (d_peaks[1].relative_magnitude > d_peaks[0].relative_magnitude) {
          swap_peaks(d_peaks, 0, 1);
        }
        if (d_peaks[2].relative_magnitude > d_peaks[0].relative_magnitude) {
          swap_peaks(d_peaks, 0, 2);
        }
        if (d_peaks[2].relative_magnitude > d_peaks[1].relative_magnitude) {
          swap_peaks(d_peaks, 2, 1);
        }
      } else if (index > 3) {
        std::sort(d_peaks.begin(), d_peaks.begin() + index,
                [] (const peak &a, const peak&b) {
                  return a.relative_magnitude > b.relative_magnitude;
                });
      }



    }

    void
    fft_burst_tagger_impl::save_peaks_to_debug_file(char * filename)
    {
      FILE * file = fopen(filename, "a");
      for (size_t index = 0; index < d_current_peaks; index++) {
        const peak& p = d_peaks[index];
        fprintf(file, "%" PRIu64 ",%d,x\n", d_abs_fft_index, p.bin);
        //float f_rel = (p.bin - d_fft_size / 2) / float(d_fft_size);
        //fprintf(file, "%f,%f,x\n", d_abs_fft_index/4e6, f_rel * 4e6 + 1624800000);
      }
      fclose(file);
    }

    void
    fft_burst_tagger_impl::tag_new_bursts(void)
    {
      for(burst b : d_new_bursts) {
        //printf("new burst %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", nitems_read(0), b.start, nitems_read(0) - b.start);
        pmt::pmt_t key = pmt::string_to_symbol("new_burst");
        //float relative_frequency = b.center_freq / d_sample_rate;
        //float relative_frequency = (b.center_bin - d_fft_size / 2) / float(d_fft_size);
        float relative_frequency = (b.center_freq) / d_sample_rate;

        // get a noise floor estimate
        //int start_bin = std::max(b.center_bin - d_burst_width / 2, 0);
        //int stop_bin = std::min(b.center_bin + d_burst_width / 2, d_fft_size-1);
        float noise_power = 0;
        for(int i = b.start_bin; i <= b.stop_bin; i++) {
            noise_power += d_baseline_sum_f[i];
        }
        noise_power /= b.stop_bin - b.start_bin + 1;
        noise_power = 10 * log10(noise_power / (d_history_size * d_fft_size));

        pmt::pmt_t value = pmt::make_dict();
        value = pmt::dict_add(value, PMT_BURST_ID, pmt::from_uint64(b.id));
        value = pmt::dict_add(value, PMT_REL_FREQ, pmt::from_float(relative_frequency));
        value = pmt::dict_add(value, PMT_CENTER_FREQ, pmt::from_float(d_center_frequency));
        value = pmt::dict_add(value, PMT_MAG, pmt::from_float(b.magnitude));
        value = pmt::dict_add(value, PMT_SAMPLE_RATE, pmt::from_float(d_sample_rate));
        value = pmt::dict_add(value, PMT_NOISE_POWER, pmt::from_float(noise_power));
        value = pmt::dict_add(value, PMT_BANDWIDTH, pmt::from_float(b.bandwidth));

        //printf("Tagging new burst %" PRIu64 " on sample %" PRIu64 " (nitems_read(0)=%" PRIu64 ")\n", b.id, b.start + d_burst_pre_len, nitems_read(0));
        add_item_tag(0, b.start * d_fft_size, key, value);
        //printf("%zu: %10f: %zu\n", b.id, relative_frequency, b.start * d_fft_size);
      }
      d_new_bursts.clear();
    }

    void
    fft_burst_tagger_impl::tag_gone_bursts(int noutput_items)
    {
      auto b = std::begin(d_gone_bursts);

      while (b != std::end(d_gone_bursts)) {
        uint64_t output_index = b->stop * d_fft_size;

        if(nitems_read(0) <= output_index && output_index < nitems_read(0) + noutput_items) {
          pmt::pmt_t key = pmt::string_to_symbol("gone_burst");
          pmt::pmt_t value = pmt::make_dict();
          value = pmt::dict_add(value, pmt::mp("burst_id"), pmt::from_uint64(b->id));
          //printf("Tagging gone burst %" PRIu64 " on sample %" PRIu64 " (nitems_read(0)=%" PRIu64 ", noutput_items=%u)\n", b->id, output_index, nitems_read(0), noutput_items);
          add_item_tag(0, output_index, key, value);
          d_n_tagged_bursts++;

          b = d_gone_bursts.erase(b);
        } else {
          ++b;
        }
      }
    }

    uint64_t
    fft_burst_tagger_impl::get_n_tagged_bursts()
    {
      return d_n_tagged_bursts;
    }

    int
    fft_burst_tagger_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      std::lock_guard<std::mutex> lock(d_work_mutex);
      d_total_timer.start();

      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      assert(noutput_items % d_fft_size == 0);

      // check for center frequency updates
      std::vector<tag_t> freq_tags;
    	get_tags_in_window(freq_tags, 0, 0, noutput_items, PMT_RX_FREQ);
      if (freq_tags.size() > 0) {
        d_center_frequency = pmt::to_double(freq_tags[freq_tags.size()-1].value);
      }

      // start at next unprocessed FFT
      int start_offset = (d_abs_fft_index + 1) * d_fft_size - nitems_read(0);
      // we have noutput_items + d_burst_pre_len*d_fft_size items in the buffer (due to history)
      for(int i = start_offset; i < noutput_items + d_burst_pre_len*d_fft_size; i += d_fft_size) {

        // keep track of current absolute FFT number
        d_abs_fft_index = (nitems_read(0) + i) / d_fft_size;
        //printf("fft_burst_tagger: nitems read: %lu, noutput_items: %d, i: %d, fft: %lu\n", nitems_read(0), noutput_items, i, d_abs_fft_index);

        // Apply a blackman window and take fft
        d_fft_timer.start();
        volk_32fc_32f_multiply_32fc(d_fft->get_inbuf(), &in[i], d_window_f, d_fft_size);
        #ifdef __USE_MKL__
        DftiComputeForward(m_fft, d_fft->get_inbuf(), d_fft->get_outbuf());
        #else
        d_fft->execute();
        #endif

        // Get the mag squared of the fft and fft shift
        size_t fftd2 = d_fft_size/2;
        volk_32fc_magnitude_squared_32f(d_magnitude_shifted_f + fftd2, d_fft->get_outbuf(), fftd2);
        volk_32fc_magnitude_squared_32f(d_magnitude_shifted_f, d_fft->get_outbuf() + fftd2, fftd2);
        d_fft_timer.end();

        if(compute_relative_magnitude()) {
          d_update_pb_timer.start();
          update_potential_bursts();
          d_update_pb_timer.end();

          d_update_ab_timer.start();
          update_active_bursts();
          d_update_ab_timer.end();

          if(d_debug) {
            extract_peaks();
            save_peaks_to_debug_file((char *)"/tmp/fft_burst_tagger-peaks.log");
          }
          d_remove_tb_timer.start();
          remove_currently_tracked_bursts();
          d_remove_tb_timer.end();

          d_extract_timer.start();
          extract_peaks();
          d_extract_timer.end();

          if(d_debug) {
            save_peaks_to_debug_file((char *)"/tmp/fft_burst_tagger-peaks-filtered.log");
          }
          d_delete_timer.start();
          delete_gone_bursts();
          d_delete_timer.end();

          d_new_pb_timer.start();
          create_new_potential_bursts();
          d_new_pb_timer.end();

          d_new_b_timer.start();
          create_new_bursts(in);
          d_new_b_timer.end();
        }

        d_update_cb_timer.start();
        update_circular_buffer();
        d_update_cb_timer.end();
      }

      d_other.start();
      int ffts_to_output;
      if (d_pre_bursts.size() != 0) {
        ffts_to_output = d_pre_bursts.begin()->start - nitems_read(0)/d_fft_size - d_rel_mag_hist;
      } else {
        ffts_to_output = noutput_items / d_fft_size - d_rel_mag_hist;
      }
      if (ffts_to_output < 0) ffts_to_output = 0;
      memcpy(out, in, sizeof(gr_complex) * ffts_to_output * d_fft_size);

      // TODO: these may write tags into the future, is that a problem?
      tag_new_bursts();
      tag_gone_bursts(noutput_items);
      d_other.end();
      d_total_timer.end();

      if (ffts_to_output == 0) usleep(100);

      // Tell runtime system how many output items we produced.
      return ffts_to_output * d_fft_size;
    }

  } /* namespace fhss_utils */
} /* namespace gr */
