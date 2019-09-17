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


#ifndef INCLUDED_FHSS_UTILS_PDU_QUADRATURE_DEMOD_CF_H
#define INCLUDED_FHSS_UTILS_PDU_QUADRATURE_DEMOD_CF_H

#include <fhss_utils/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace fhss_utils {

    /*!
     * \brief <+description of block+>
     * \ingroup fhss_utils
     *
     */
    class FHSS_UTILS_API pdu_quadrature_demod_cf : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<pdu_quadrature_demod_cf> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of fhss_utils::pdu_quadrature_demod_cf.
       *
       * To avoid accidental use of raw pointers, fhss_utils::pdu_quadrature_demod_cf's
       * constructor is in a private implementation
       * class. fhss_utils::pdu_quadrature_demod_cf::make is the public interface for
       * creating new instances.
       */
      static sptr make(const std::vector<float> &taps, float sensitivity);
    };

  } // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_PDU_QUADRATURE_DEMOD_CF_H */

