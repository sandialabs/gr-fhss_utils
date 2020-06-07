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

#ifndef INCLUDED_FHSS_UTILS_BURST_DOWNMIX_IMPL_H
#define INCLUDED_FHSS_UTILS_BURST_DOWNMIX_IMPL_H

#include <gnuradio/blocks/rotator.h>
#include <gnuradio/filter/fir_filter.h>

#include <fhss_utils/burst_downmix.h>

namespace gr {
  namespace fhss_utils {

    class burst_downmix_impl : public burst_downmix
    {
     private:

      int d_decimation;
      size_t d_max_burst_size;

      gr_complex * d_frame;
      gr_complex * d_tmp_a;

      filter::kernel::fir_filter_ccf d_input_fir;

      blocks::rotator d_r;

      void handler(pmt::pmt_t msg);

      void update_buffer_sizes(size_t burst_size);

     public:
      burst_downmix_impl(const std::vector<float> &taps, int decimation);
      ~burst_downmix_impl();

      size_t get_input_queue_size();
      uint64_t get_n_dropped_bursts();

    };

  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_BURST_DOWNMIX_IMPL_H */
