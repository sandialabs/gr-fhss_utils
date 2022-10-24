/* -*- c++ -*- */
/*
 * Copyright 2018-2021 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/fhss_utils/constants.h>
#include <gnuradio/io_signature.h>

namespace gr {
namespace fhss_utils {

const pmt::pmt_t PMTCONSTSTR__in()
{
    static const pmt::pmt_t val = pmt::mp("in");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__out()
{
    static const pmt::pmt_t val = pmt::mp("out");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__center_frequency()
{
    static const pmt::pmt_t val = pmt::mp("center_frequency");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__relative_frequency()
{
    static const pmt::pmt_t val = pmt::mp("relative_frequency");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__sample_rate()
{
    static const pmt::pmt_t val = pmt::mp("sample_rate");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__bandwidth()
{
    static const pmt::pmt_t val = pmt::mp("bandwidth");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__pwr_db()
{
    static const pmt::pmt_t val = pmt::mp("pwr_db");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__snr_db()
{
    static const pmt::pmt_t val = pmt::mp("snr_db");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__debug()
{
    static const pmt::pmt_t val = pmt::mp("debug");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__rx_freq()
{
    static const pmt::pmt_t val = pmt::mp("rx_freq");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__burst_id()
{
    static const pmt::pmt_t val = pmt::mp("burst_id");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__magnitude()
{
    static const pmt::pmt_t val = pmt::mp("magnitude");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__noise_density()
{
    static const pmt::pmt_t val = pmt::mp("noise_density");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__new_burst()
{
    static const pmt::pmt_t val = pmt::mp("new_burst");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__gone_burst()
{
    static const pmt::pmt_t val = pmt::mp("gone_burst");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__rx_time()
{
    static const pmt::pmt_t val = pmt::mp("rx_time");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__start_time()
{
    static const pmt::pmt_t val = pmt::mp("start_time");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__duration()
{
    static const pmt::pmt_t val = pmt::mp("duration");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__cpdus()
{
    static const pmt::pmt_t val = pmt::mp("cpdus");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__start_offset()
{
    static const pmt::pmt_t val = pmt::mp("start_offset");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__end_offset()
{
    static const pmt::pmt_t val = pmt::mp("end_offset");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__input_rate()
{
    static const pmt::pmt_t val = pmt::mp("input_rate");
    return val;
}
const pmt::pmt_t PMTCONSTSTR__cut_short()
{
    static const pmt::pmt_t val = pmt::mp("cut_short");
    return val;
}

} /* namespace fhss_utils */
} /* namespace gr */
