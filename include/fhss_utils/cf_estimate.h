/* -*- c++ -*- */
/* 
 * Copyright 2020 gr-fhss_utils author.
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


#ifndef INCLUDED_FHSS_UTILS_CF_ESTIMATE_H
#define INCLUDED_FHSS_UTILS_CF_ESTIMATE_H

#include <fhss_utils/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace fhss_utils {

    /*!
     * \brief <+description of block+>
     * \ingroup fhss_utils
     *
     */
    class FHSS_UTILS_API cf_estimate : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<cf_estimate> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of fhss_utils::cf_estimate.
       *
       * To avoid accidental use of raw pointers, fhss_utils::cf_estimate's
       * constructor is in a private implementation
       * class. fhss_utils::cf_estimate::make is the public interface for
       * creating new instances.
       */
      static sptr make(int method=0, std::vector<float> channel_freqs = std::vector<float>());

      /*!
       * \brief Set channel center frequencies
       *
       * If the selected method is COERCE, the center frequency for each
       * burst will be coerced to the nearest channel frequency.
       *
       * \param channel_freqs Channel frequencies to set.
       */
      virtual void set_freqs(std::vector<float> channel_freqs) = 0;

      /*!
       * \brief Set center frequency estimation algorithm
       *
       * COERCE, RMS, HALF_POWER are the valid options, and are enumerated
       * types in the fhss_utils namespace
       *
       * \param int method Method enum to select
       */
      virtual void set_method(int method) = 0;
    };


    /*!
     * \brief Center frequency estimation method options
     *
     */
		enum cf_method
		{
			RMS = 0,
			HALF_POWER,
			COERCE
		};


  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_CF_ESTIMATE_H */

