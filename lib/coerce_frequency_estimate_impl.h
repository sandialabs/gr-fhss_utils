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

#ifndef INCLUDED_FHSS_UTILS_COERCE_FREQUENCY_ESTIMATE_IMPL_H
#define INCLUDED_FHSS_UTILS_COERCE_FREQUENCY_ESTIMATE_IMPL_H

#include <fhss_utils/coerce_frequency_estimate.h>

namespace gr {
  namespace fhss_utils {

    class coerce_frequency_estimate_impl : public coerce_frequency_estimate
    {
     private:
      std::vector<float> d_valid_freqs;
      uint32_t d_n_freqs;

     public:
      coerce_frequency_estimate_impl(std::vector<float> valid_freqs);
      ~coerce_frequency_estimate_impl();

      void handler(pmt::pmt_t msg);
      void set_freqs(const std::vector<float>);

    };

  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_COERCE_FREQUENCY_ESTIMATE_IMPL_H */
