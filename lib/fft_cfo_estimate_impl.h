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

#ifndef INCLUDED_FHSS_UTILS_FFT_CFO_ESTIMATE_IMPL_H
#define INCLUDED_FHSS_UTILS_FFT_CFO_ESTIMATE_IMPL_H

#include <fhss_utils/fft_cfo_estimate.h>
#include <gnuradio/blocks/rotator.h>
#include <gnuradio/fft/fft.h>

namespace gr {
  namespace fhss_utils {

    class fft_cfo_estimate_impl : public fft_cfo_estimate
    {
     private:
       int d_fftsize;
       int d_offset;
       int d_nffts;
       float d_threshold;

       std::vector<float> d_valid_freqs;
       std::vector<gr_complex> d_output;

       float* d_mag;
       float* d_mag_sum;

       fft::fft_complex d_fft;
       blocks::rotator d_r;

       double get_closest_freq(double);
       void pdu_handler(pmt::pmt_t pdu);
       float calc_center_freq(float sample_rate);
       float calc_bandwidth(float& power, float sample_rate);

     public:
      fft_cfo_estimate_impl(int fftsize, int offset, int num_ffts, float threshold, std::vector<float> valid_freqs);
      ~fft_cfo_estimate_impl();

      void set_freqs(std::vector<float> valid_freqs);
    };

  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_FFT_CFO_ESTIMATE_IMPL_H */
