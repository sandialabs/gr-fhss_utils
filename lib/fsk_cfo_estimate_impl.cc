/* -*- c++ -*- */
/* 
 * Copyright 2018 gr-fhss_utils author.
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
#include <gnuradio/math.h>
#include <volk/volk.h>
#include "fsk_cfo_estimate_impl.h"

namespace gr {
  namespace fhss_utils {

    fsk_cfo_estimate::sptr
      fsk_cfo_estimate::make(int search_depth, int start_offset, int samples_to_average, float threshold, const std::vector<float> &start_finder_taps)
    {
      return gnuradio::get_initial_sptr
        (new fsk_cfo_estimate_impl(search_depth, start_offset, samples_to_average, threshold, start_finder_taps));
    }

    /*
     * The private constructor
     */
    fsk_cfo_estimate_impl::fsk_cfo_estimate_impl(int search_depth, int start_offset, int samples_to_average, float threshold, const std::vector<float> &start_finder_taps)
      : gr::block("fsk_cfo_estimate",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0)),
      d_offset(start_offset),
      d_n_to_average(samples_to_average),
      d_search_depth(search_depth),
      d_threshold(threshold),
      d_max_burst_size(0),
      d_magnitude_f(nullptr),
      d_magnitude_filtered_f(nullptr),
      d_tmp_c(nullptr),
      d_start_finder_fir(1, start_finder_taps)
    {
      message_port_register_in(pmt::mp("pdu_in"));
      message_port_register_out(pmt::mp("pdu_out"));
      set_msg_handler(pmt::mp("pdu_in"), boost::bind(&fsk_cfo_estimate_impl::pdu_handler, this, _1));
    }

    /*
     * Our virtual destructor.
     */
    fsk_cfo_estimate_impl::~fsk_cfo_estimate_impl()
    {
      if(d_magnitude_f) {
        volk_free(d_magnitude_f);
      }
      if(d_magnitude_filtered_f) {
        volk_free(d_magnitude_filtered_f);
      }
      if(d_tmp_c) {
        volk_free(d_tmp_c);
      }
    }

    void fsk_cfo_estimate_impl::update_buffer_sizes(size_t burst_size)
    {
      if(burst_size > d_max_burst_size) {
        d_max_burst_size = burst_size;

        if(d_tmp_c) {
          volk_free(d_tmp_c);
        }
        d_tmp_c = (gr_complex *)volk_malloc(d_max_burst_size * sizeof(gr_complex), volk_get_alignment());

        if(d_magnitude_f) {
          volk_free(d_magnitude_f);
        }
        d_magnitude_f = (float *)volk_malloc(d_max_burst_size * sizeof(float), volk_get_alignment());

        if(d_magnitude_filtered_f) {
          volk_free(d_magnitude_filtered_f);
        }
        d_magnitude_filtered_f = (float *)volk_malloc(d_max_burst_size * sizeof(float), volk_get_alignment());
      }
    }

    void fsk_cfo_estimate_impl::pdu_handler(pmt::pmt_t pdu) {
      if (!pmt::is_pair(pdu)) {
        GR_LOG_WARN(d_logger, "PDU not a pair, dropping");
        return;
      }

      pmt::pmt_t metadata = pmt::car(pdu);
      pmt::pmt_t data = pmt::cdr(pdu);

      if (!pmt::is_dict(metadata)) {
        GR_LOG_WARN(d_logger, "PDU metadata not a dictionary, dropping");
        return;
      }

      float center_frequency = pmt::to_float(pmt::dict_ref(metadata, pmt::mp("center_frequency"), pmt::PMT_NIL));
      float relative_frequency = pmt::to_float(pmt::dict_ref(metadata, pmt::mp("relative_frequency"), pmt::PMT_NIL));
      float sample_rate = pmt::to_float(pmt::dict_ref(metadata, pmt::mp("sample_rate"), pmt::PMT_NIL));

      if (!pmt::is_c32vector(data)) {
        GR_LOG_WARN(d_logger, "PDU data not complex, dropping");
        return;
      }

      size_t burst_size = pmt::length(data);
      const gr_complex* burst = (const gr_complex*) pmt::c32vector_elements(data, burst_size);

      update_buffer_sizes(burst_size + 10000);

      /*
       * Estimate fine CFO by averaging the FM demodulated burst, looking
       * up to d_search_depth samples into the burst
       */
      int half_fir_size = (d_start_finder_fir.ntaps() - 1) / 2;

      // The burst might be shorter than d_search_depth.
      int n_to_search = std::min(d_search_depth, (int)burst_size);

      // calc magnitudes
      volk_32fc_magnitude_32f(d_magnitude_f + half_fir_size, burst, n_to_search);

      // zero padding
      memset(d_magnitude_f, 0, sizeof(float) * half_fir_size);
      memset(d_magnitude_f + half_fir_size + n_to_search, 0, sizeof(float) * half_fir_size);

      d_start_finder_fir.filterN(d_magnitude_filtered_f, d_magnitude_f, n_to_search);

      float* max = std::max_element(d_magnitude_filtered_f, d_magnitude_filtered_f + n_to_search);
      float threshold = *max * d_threshold;

      int start;
      for(start = 0; start < n_to_search; start++) {
        if(d_magnitude_filtered_f[start] >= threshold) {
          break;
        }
      }

      start += d_offset;
      if (start >= burst_size) {
        GR_LOG_WARN(d_logger, "could not find start, dropping burst");
        return;
      }

      // Need the minus one here, because we are doing the offset conjugate mult below and don't want
      // to overrun our array.
      int n_samps_to_avg(std::min(int(burst_size - start - 1), d_n_to_average));
      volk_32fc_x2_multiply_conjugate_32fc(&d_tmp_c[0], &burst[start+1], &burst[start], n_samps_to_avg);
      float sum = 0.0;
      for(int i = 0; i < n_samps_to_avg; i++) {
        sum += gr::fast_atan2f(std::imag(d_tmp_c[i]), std::real(d_tmp_c[i]));
      }
      sum /= n_samps_to_avg;

      /*
       * Apply fine frequency correction
       */
      d_r.set_phase_incr(exp(gr_complex(0, -sum)));
      d_r.set_phase(gr_complex(1, 0));
      d_r.rotateN(d_tmp_c, burst, burst_size);
      center_frequency += (sum / M_PI / 2.0) * sample_rate;
      relative_frequency += (sum / M_PI / 2.0) * sample_rate;

      // build and publish the new PDU
      pmt::pmt_t pdu_vector = pmt::init_c32vector(burst_size, d_tmp_c);

      metadata = pmt::dict_delete(metadata, pmt::mp("center_frequency"));
      metadata = pmt::dict_add(metadata, pmt::mp("center_frequency"), pmt::mp(center_frequency));
      metadata = pmt::dict_add(metadata, pmt::mp("relative_frequency"), pmt::mp(relative_frequency));

      pmt::pmt_t out_msg = pmt::cons(metadata, pdu_vector);
      message_port_pub(pmt::mp("pdu_out"), out_msg);

      return;
    }

  } /* namespace fhss_utils */
} /* namespace gr */

