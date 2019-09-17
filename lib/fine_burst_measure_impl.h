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

#ifndef INCLUDED_FHSS_UTILS_FINE_BURST_MEASURE_IMPL_H
#define INCLUDED_FHSS_UTILS_FINE_BURST_MEASURE_IMPL_H

#include <fhss_utils/fine_burst_measure.h>
#include <gnuradio/fft/fft.h>
#include <gnuradio/blocks/rotator.h>

namespace gr {
  namespace fhss_utils {

    struct peak {
      int64_t center_bin;
      uint64_t value;
    };

    struct fsk_info {
      std::vector<float> freqs;
      float shift;
    };

    class fine_burst_measure_impl : public fine_burst_measure
    {
     private:
      uint64_t d_fft_size;
      float d_freq_resolution;
      uint32_t d_analyze_samples;
      float d_threshold;
      float d_min_freq;
      float d_max_freq;

      float* d_phase_diff_f;
      gr_complex* d_tmp_c;
      std::vector<gr_complex> d_tmp_burst;

      blocks::rotator d_r;

      std::vector<uint32_t> d_histogram;
      std::vector<uint32_t> d_histogram_smooth;
      std::vector<size_t> d_fsk_peaks;
      std::vector<peak> d_peak_list;

      gr::fft::fft_complex *d_fft;
      float * d_magnitude_f;
      float * d_magnitude_shifted_f;
      float * d_magnitude_smooth_f;
      float * d_window_f;

      void pdu_handler(pmt::pmt_t pdu);
      void resize_bins(size_t numBins);

      float bin2Freq(size_t bin, size_t numBins) {return ((float)bin - numBins / 2) / numBins * 2 * M_PI;} 
      fsk_info findFreqOffset(float* input, size_t length, const std::vector<float>& freqs);
      float findPower(const gr_complex* in, const std::vector<float>& freqs, size_t inSize);
      float findOccupiedBW(float power, size_t inSize);

      FILE* fid;
      FILE* fid2;

     public:
      fine_burst_measure_impl(float freq_resolution, uint32_t analyze_samples, float threshold);
      ~fine_burst_measure_impl();
    };

  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_FINE_BURST_MEASURE_IMPL_H */


