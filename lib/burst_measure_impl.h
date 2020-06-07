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

#ifndef INCLUDED_FHSS_UTILS_BURST_MEASURE_IMPL_H
#define INCLUDED_FHSS_UTILS_BURST_MEASURE_IMPL_H

#include <fhss_utils/burst_measure.h>
#include <gnuradio/fft/fft.h>

namespace gr {
  namespace fhss_utils {

    class burst_measure_impl : public burst_measure
    {
     private:
      void handler(pmt::pmt_t msg);

      float* d_mags;

      std::vector<gr::fft::fft_complex*> d_ffts;
      std::vector<float*> d_windows;

      void fft_setup(int power);
      void fft_cleanup();

     public:
      burst_measure_impl();
      ~burst_measure_impl();
    };

  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_BURST_MEASURE_IMPL_H */

// vim: ts=2:sw=2:sts=2

