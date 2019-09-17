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
#include <gnuradio/math.h>
#include "pdu_quadrature_demod_cf_impl.h"

#include <volk/volk.h>
#include <inttypes.h>


namespace gr {
  namespace fhss_utils {

    pdu_quadrature_demod_cf::sptr
    pdu_quadrature_demod_cf::make(const std::vector<float> &taps, float sensitivity)
    {
      return gnuradio::get_initial_sptr
        (new pdu_quadrature_demod_cf_impl(taps, sensitivity));
    }

    /*
     * The private constructor
     */
    pdu_quadrature_demod_cf_impl::pdu_quadrature_demod_cf_impl(const std::vector<float> &taps, float sensitivity)
      : gr::block("pdu_quadrature_demod_cf",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0)),
        d_input_fir(0, taps),
        d_sensitivity(sensitivity)
    {
      message_port_register_in(pmt::mp("cpdus"));
      message_port_register_out(pmt::mp("fpdus"));
      set_msg_handler(pmt::mp("cpdus"), boost::bind(&pdu_quadrature_demod_cf_impl::handle_pdu, this, _1));
    }

    /*
     * Our virtual destructor.
     */
    pdu_quadrature_demod_cf_impl::~pdu_quadrature_demod_cf_impl()
    {
    }

    void
    pdu_quadrature_demod_cf_impl::handle_pdu(pmt::pmt_t pdu)
    {
      pmt::pmt_t samples = pmt::cdr(pdu);
      size_t burst_size;
      const gr_complex * burst = (const gr_complex*)pmt::c32vector_elements(samples, burst_size);
      // Subtract off 1 from burst size because we do an offset for the conjugate multiply.
      burst_size--;

      std::vector<gr_complex> tmp(burst_size);
      std::vector<float> out(burst_size);

      volk_32fc_x2_multiply_conjugate_32fc(&tmp[0], &burst[1], &burst[0], burst_size);
      for(int i = 0; i < burst_size; i++) {
        out[i] = d_sensitivity * gr::fast_atan2f(imag(tmp[i]), real(tmp[i]));
      }


      pmt::pmt_t pdu_vector = pmt::init_f32vector(burst_size, out);
      pmt::pmt_t out_msg = pmt::cons(pmt::car(pdu),
          pdu_vector);
      message_port_pub(pmt::mp("fpdus"), out_msg);


      //d_input_fir.filterNdec(d_frame, d_tmp_a, burst_size, decimation);
    }


  } /* namespace fhss_utils */
} /* namespace gr */
