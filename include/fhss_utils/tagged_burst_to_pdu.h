/* -*- c++ -*- */
/*
 * Copyright 2018-2021 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_FHSS_UTILS_TAGGED_BURST_TO_PDU_H
#define INCLUDED_FHSS_UTILS_TAGGED_BURST_TO_PDU_H

#include <gnuradio/sync_block.h>
#include <fhss_utils/api.h>

namespace gr {
namespace fhss_utils {

/*!
 * \brief Tagged Burst to PDU
 * \ingroup fhss_utils
 *
 * Converts Tagged burst stream to PDU messages
 */
class FHSS_UTILS_API tagged_burst_to_pdu : virtual public gr::sync_block
{
public:
    typedef std::shared_ptr<tagged_burst_to_pdu> sptr;

    /**
     * \brief Creates a new instance of fhss_utils::tagged_burst_to_pdu.
     *
     * @param decimation -
     * @param taps -
     * @param min_burst_time -
     * @param max_burst_time -
     * @param relative_center_frequency -
     * @param relative_span -
     * @param relative_sample_rate -
     * @param sample_rate -
     * @param num_threads -
     */
    static sptr make(size_t decimation,
                     const std::vector<float>& taps,
                     float min_burst_time,
                     float max_burst_time,
                     float relative_center_frequency,
                     float relative_span,
                     float relative_sample_rate,
                     float sample_rate,
                     int num_threads);

    // virtual uint64_t get_n_dropped_bursts() = 0;
    // virtual int get_output_queue_size() = 0;
    // virtual int get_output_max_queue_size() = 0;
};

} // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_TAGGED_BURST_TO_PDU_H */
