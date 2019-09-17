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

#ifndef INCLUDED_FHSS_UTILS_FSK_CFO_ESTIMATE_IMPL_H
#define INCLUDED_FHSS_UTILS_FSK_CFO_ESTIMATE_IMPL_H

#include <fhss_utils/fsk_cfo_estimate.h>
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/blocks/rotator.h>

namespace gr {
  namespace fhss_utils {

    class fsk_cfo_estimate_impl : public fsk_cfo_estimate
    {
     private:
      int d_offset;
      int d_n_to_average;
      int d_search_depth;
      float d_threshold;
      size_t d_max_burst_size;
      filter::kernel::fir_filter_fff d_start_finder_fir;

      float* d_magnitude_f;
      float* d_magnitude_filtered_f;
      gr_complex* d_tmp_c;

      blocks::rotator d_r;

      void pdu_handler(pmt::pmt_t pdu);
      void update_buffer_sizes(size_t burst_size);

     public:
      fsk_cfo_estimate_impl(int search_depth, int start_offset, int samples_to_average, float threshold, const std::vector<float> &start_finder_taps);
      ~fsk_cfo_estimate_impl();
    };

  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_FSK_CFO_ESTIMATE_IMPL_H */

