/* -*- c++ -*- */
/*
 * Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
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
#include <gnuradio/math.h>

#include "burst_downmix_impl.h"

#include <gnuradio/filter/firdes.h>
#include <volk/volk.h>

#include <inttypes.h>
#include <chrono>

namespace gr {
  namespace fhss_utils {

    burst_downmix::sptr
    burst_downmix::make(const std::vector<float> &taps, int decimation)
    {
      return gnuradio::get_initial_sptr
        (new burst_downmix_impl(taps, decimation));
    }

    /*
     * The private constructor
     */
    burst_downmix_impl::burst_downmix_impl(const std::vector<float> &taps, int decimation)
      : gr::block("burst_downmix",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0)),
          d_decimation(decimation),
          d_max_burst_size(0),

          d_frame(NULL),
          d_tmp_a(NULL),

          d_input_fir(0, taps)
    {
      message_port_register_in(pmt::mp("cpdus"));
      message_port_register_out(pmt::mp("cpdus"));

      set_msg_handler(pmt::mp("cpdus"), boost::bind(&burst_downmix_impl::handler, this, _1));
    }

    /*
     * Our virtual destructor.
     */
    burst_downmix_impl::~burst_downmix_impl()
    {
      if(d_frame) {
        volk_free(d_frame);
      }
      if(d_tmp_a) {
        volk_free(d_tmp_a);
      }
    }

    void burst_downmix_impl::update_buffer_sizes(size_t burst_size)
    {
      if(burst_size > d_max_burst_size) {
        d_max_burst_size = burst_size;
        if(d_frame) {
          volk_free(d_frame);
        }
        d_frame = (gr_complex *)volk_malloc(d_max_burst_size * sizeof(gr_complex), volk_get_alignment());

        if(d_tmp_a) {
          volk_free(d_tmp_a);
        }
        d_tmp_a = (gr_complex *)volk_malloc(d_max_burst_size * sizeof(gr_complex), volk_get_alignment());
      }
    }

    void burst_downmix_impl::handler(pmt::pmt_t msg)
    {
      /*
       * Extract the burst and meta data from the cpdu
       */
      pmt::pmt_t samples = pmt::cdr(msg);
      size_t burst_size = pmt::length(samples);
      const gr_complex * burst = (const gr_complex*)pmt::c32vector_elements(samples, burst_size);

      pmt::pmt_t meta = pmt::car(msg);
      float relative_frequency = pmt::to_float(pmt::dict_ref(meta, pmt::mp("relative_frequency"), pmt::PMT_NIL));
      float center_frequency = pmt::to_float(pmt::dict_ref(meta, pmt::mp("center_frequency"), pmt::PMT_NIL));
      float sample_rate = pmt::to_float(pmt::dict_ref(meta, pmt::mp("sample_rate"), pmt::PMT_NIL));

      // This burst might be larger than the one before.
      // Update he buffer sizes if needed.
      update_buffer_sizes(burst_size + 10000);

      /*
       * Shift the center frequency of the burst to the provided rough CFO estimate.
       */
      float phase_inc = 2 * M_PI * -relative_frequency;
      d_r.set_phase_incr(exp(gr_complex(0, phase_inc)));
      d_r.set_phase(gr_complex(1, 0));
      d_r.rotateN(d_tmp_a, burst, burst_size);
      center_frequency += relative_frequency * sample_rate;

      /*
       * Apply the low pass filter and decimate the burst.
       */
      burst_size = (burst_size - d_input_fir.ntaps() + 1) / d_decimation;
      d_input_fir.filterNdec(d_frame, d_tmp_a, burst_size, d_decimation);
      sample_rate /= d_decimation;

      pmt::pmt_t pdu_vector = pmt::init_c32vector(burst_size, d_frame);

      meta = pmt::dict_delete(meta, pmt::mp("sample_rate"));
      meta = pmt::dict_add(meta, pmt::mp("sample_rate"), pmt::mp(sample_rate));
      meta = pmt::dict_delete(meta, pmt::mp("center_frequency"));
      meta = pmt::dict_add(meta, pmt::mp("center_frequency"), pmt::mp(center_frequency));
      meta = pmt::dict_delete(meta, pmt::mp("relative_frequency")); // no longer needed, since center frequency reflects it

      pmt::pmt_t out_msg = pmt::cons(meta, pdu_vector);
      message_port_pub(pmt::mp("cpdus"), out_msg);

      return;
    }

  } /* namespace fhss_utils */
} /* namespace gr */

