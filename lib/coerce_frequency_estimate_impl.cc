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
#include "coerce_frequency_estimate_impl.h"

namespace gr {
  namespace fhss_utils {

    coerce_frequency_estimate::sptr
    coerce_frequency_estimate::make(std::vector<float> valid_freqs)
    {
      return gnuradio::get_initial_sptr
        (new coerce_frequency_estimate_impl(valid_freqs));
    }

    /*
     * The private constructor
     */
    coerce_frequency_estimate_impl::coerce_frequency_estimate_impl(std::vector<float> valid_freqs)
      : gr::block("coerce_frequency_estimate",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0))
    {
      set_freqs(valid_freqs);

      if (d_n_freqs == 0) {
        GR_LOG_WARN(d_logger, "No frequencies provided! Block effectively disabled!");
      }

      message_port_register_in(pmt::mp("cpdus"));
      message_port_register_out(pmt::mp("cpdus"));
      set_msg_handler(pmt::mp("cpdus"), boost::bind(&coerce_frequency_estimate_impl::handler, this, _1));
    }

    /*
     * Our virtual destructor.
     */
    coerce_frequency_estimate_impl::~coerce_frequency_estimate_impl()
    {
    }


    void coerce_frequency_estimate_impl::set_freqs(const std::vector<float> freqs) {
      gr::thread::scoped_lock l(d_setlock);

      d_valid_freqs = freqs;
      d_n_freqs = d_valid_freqs.size();
    }


    void coerce_frequency_estimate_impl::handler(pmt::pmt_t pdu) {
      // if there aren't any valid frequencies loaded, just pass the PDU on and exit
      gr::thread::scoped_lock l(d_setlock);

      if (d_valid_freqs.size() == 0) {
        message_port_pub(pmt::mp("cpdus"), pdu);
        return;
      }

      // make sure PDU data is formed properly
      if (!(pmt::is_pair(pdu))) {
        GR_LOG_WARN(d_logger, "received unexpected PMT (non-pair)");
        return;
      }

      pmt::pmt_t meta  = pmt::car(pdu);

      if (pmt::is_dict(meta)) {

        double sample_rate = pmt::to_float(pmt::dict_ref(meta, pmt::mp("sample_rate"), pmt::PMT_NIL));
        double relative_freq = pmt::to_float(pmt::dict_ref(meta, pmt::mp("relative_frequency"), pmt::PMT_NIL));
        double center_freq = pmt::to_float(pmt::dict_ref(meta, pmt::mp("center_frequency"), pmt::PMT_NIL));
        double actual_freq_est = center_freq + sample_rate * relative_freq;

        double min_diff = std::abs(actual_freq_est - d_valid_freqs[0]);
        double min_freq = d_valid_freqs[0];
        for (int ii=1; ii<d_n_freqs; ii++) {
          double diff = std::abs(d_valid_freqs[ii]-actual_freq_est);
          if (diff < min_diff) {
            min_diff = diff;
            min_freq = d_valid_freqs[ii];
          }
        }

        GR_LOG_DEBUG(d_logger, boost::format("actual frequency estimate is %d kHz  --  closest valid freq is %d kHz  --  error is %d kHz")
                                          % int(actual_freq_est/1000.0)
                                          % int(min_freq/1000.0)
                                          % int(std::abs(actual_freq_est-min_freq)/1000));

        double relative_error = (actual_freq_est-min_freq) / sample_rate;

        meta = pmt::dict_add(meta, pmt::mp("relative_frequency"), pmt::from_float(relative_freq - relative_error));
        message_port_pub(pmt::mp("cpdus"), pmt::cons(meta, pmt::cdr(pdu)));

      } else {
        GR_LOG_WARN(d_logger, "received malformed PDU");
      }

    }

  } /* namespace fhss_utils */
} /* namespace gr */
