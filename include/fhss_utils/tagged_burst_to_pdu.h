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


#ifndef INCLUDED_FHSS_UTILS_TAGGED_BURST_TO_PDU_H
#define INCLUDED_FHSS_UTILS_TAGGED_BURST_TO_PDU_H

#include <fhss_utils/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace fhss_utils {

    /*!
     * \brief <+description of block+>
     * \ingroup fhss_utils
     *
     */
    class FHSS_UTILS_API tagged_burst_to_pdu : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<tagged_burst_to_pdu> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of fhss_utils::tagged_burst_to_pdu.
       *
       * To avoid accidental use of raw pointers, fhss_utils::tagged_burst_to_pdu's
       * constructor is in a private implementation
       * class. fhss_utils::tagged_burst_to_pdu::make is the public interface for
       * creating new instances.
       */
      static sptr make(size_t decimation, const std::vector<float>& taps, float min_burst_time, float max_burst_time, 
          float relative_center_frequency, float relative_span, float relative_sample_rate, float sample_rate, int num_threads);

      //virtual uint64_t get_n_dropped_bursts() = 0;
      //virtual int get_output_queue_size() = 0;
      //virtual int get_output_max_queue_size() = 0;

    };

  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_TAGGED_BURST_TO_PDU_H */
