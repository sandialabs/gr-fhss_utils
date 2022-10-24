/* -*- c++ -*- */
/*
 * Copyright 2018-2021 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_FHSS_UTILS_CF_ESTIMATE_H
#define INCLUDED_FHSS_UTILS_CF_ESTIMATE_H

#include <gnuradio/block.h>
#include <gnuradio/fhss_utils/api.h>

namespace gr {
namespace fhss_utils {

/*!
 * \brief Estimates center frequency of a burst and re-centers it.
 * \ingroup fhss_utils
 *\n
 * Four Methods are provided for Center Frequency Estimation:\n
 * -Root-Mean_Squared: RMS estimate, this is the best for general use.\n
 * -Half Power: uses the center point in the cumulative sum of the PSD, this
 * method is not great for general use.\n
 * -Middle Out: uses the relationship between the peak and noise floor to
 * estimate the bandwidth, this method is best for high SNR signals. \n
 * -Coerce Only moves center freq to the nearest entry in the provided
 * channel_freqs list\n
 * Frequency coercion also be applied to any method if the list is provided.
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
    typedef std::shared_ptr<cf_estimate> sptr;

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

    /*!
     * \brief Set the expected minimum SNR (only used by middle out method).
     *
     * \param snr Expected minimum SNR.
     */
    virtual void set_snr_min(float snr) = 0;
    /*!
     * \brief Set the minimum threshold to use (only used by middle out method).
     *
     * \param thresh Minimum threshold to use relative to the peak.
     */
    virtual void set_thresh_min(float thresh) = 0;
};


/*!
 * \brief Center frequency estimation method options
 *
 */
enum cf_method { RMS = 0, HALF_POWER, MIDDLE_OUT, COERCE };


} // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_CF_ESTIMATE_H */
