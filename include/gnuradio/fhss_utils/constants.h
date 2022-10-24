/* -*- c++ -*- */
/*
 * Copyright 2018-2021 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */


#ifndef INCLUDED_FHSS_UTILS_CONSTANTS_H
#define INCLUDED_FHSS_UTILS_CONSTANTS_H

#include <gnuradio/fhss_utils/api.h>
#include <pmt/pmt.h>

namespace gr {
namespace fhss_utils {

FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__in();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__out();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__center_frequency();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__relative_frequency();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__sample_rate();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__bandwidth();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__pwr_db();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__snr_db();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__debug();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__rx_freq();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__burst_id();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__magnitude();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__noise_density();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__new_burst();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__gone_burst();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__rx_time();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__start_time();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__duration();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__cpdus();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__start_offset();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__end_offset();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__input_rate();
FHSS_UTILS_API const pmt::pmt_t PMTCONSTSTR__cut_short();

} // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_CONSTANTS_H */
