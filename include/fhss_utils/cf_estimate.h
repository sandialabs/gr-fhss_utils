/* -*- c++ -*- */
/*
 * Copyright 2018, 2019, 2020 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_FHSS_UTILS_CF_ESTIMATE_H
#define INCLUDED_FHSS_UTILS_CF_ESTIMATE_H

#include <gnuradio/block.h>
#include <fhss_utils/api.h>

namespace gr {
namespace fhss_utils {

/*!
 * \brief Estimates center frequency of a burst and re-centers it.
 * \ingroup fhss_utils
 *\n
 * Three Methods are provided for Center Frequency Estimating\n
 * -Coerce quantizes to the nearest frequency in the channel_freqs list\n
 * -Root-Mean_Squared method\n
 * -Half Power is a half-power point in the cumulative sum of the PSD\n
 *\n
 * Parameters:\n
 * -method selects method of center frequency estimation.\n
 *\n
 * Required metadata: 'sample_rate' and 'center_frequency'\n
 *\n
 * Optional metadata: 'relative_frequency'\n
 *\n
 * Produced metdata: 'bandwidth', 'snr_db', corrections to 'center_frequency' and
 *'relative_frequency'\n
 *
 */
class FHSS_UTILS_API cf_estimate : virtual public gr::block
{
public:
    typedef boost::shared_ptr<cf_estimate> sptr;

    /*!
     * \brief Creates a new instance of fhss_utils::cf_estimate.
     *
     *
     * \param method Center Frequency Estimation method #cf_method
     * \param channel_freqs Channel frequencies to set for coerce method.
     */
    static sptr make(int method = 0,
                     std::vector<float> channel_freqs = std::vector<float>());

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
     * \param method Method enum to select
     */
    virtual void set_method(int method) = 0;
};


/*!
 * \brief Center frequency estimation method options
 *
 */
enum cf_method { RMS = 0, HALF_POWER, COERCE };


} // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_CF_ESTIMATE_H */
