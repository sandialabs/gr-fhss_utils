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

#ifndef INCLUDED_FHSS_UTILS_PDU_QUADRATURE_DEMOD_CF_IMPL_H
#define INCLUDED_FHSS_UTILS_PDU_QUADRATURE_DEMOD_CF_IMPL_H

#include <fhss_utils/pdu_quadrature_demod_cf.h>
#include <gnuradio/filter/fir_filter.h>

namespace gr {
  namespace fhss_utils {

    class pdu_quadrature_demod_cf_impl : public pdu_quadrature_demod_cf
    {
     private:
       filter::kernel::fir_filter_ccf d_input_fir;
       float d_sensitivity;

     public:
      pdu_quadrature_demod_cf_impl(const std::vector<float> &taps, float sensitivity);
      ~pdu_quadrature_demod_cf_impl();

      void handle_pdu(pmt::pmt_t pdu);

    };

  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_PDU_QUADRATURE_DEMOD_CF_IMPL_H */
