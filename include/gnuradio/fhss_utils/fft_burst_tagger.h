/* -*- c++ -*- */
/*
 * Copyright 2018-2021 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_H
#define INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_H

#include <gnuradio/fhss_utils/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
namespace fhss_utils {

/*!
 * \brief This block is a fast energy detector that processes a stream of complex data and
 * will apply stream tags to identify the start, end, and approximate frequency of
 * detected bursts of energy.
 *
 * Internally this block processes data one FFT at a time. It maintains a dynamic noise
 * floor estimate for each bin over the prior `history_size` FFTs, and compares incoming
 * FFT bins to see if they are `threshold` dB above the moving noise floor estimate. Once
 * detected, potential bursts are observed to ensure they remain above the threshold for
 * at least `lookahead` FFTs, after which they are tracked as detections and tagged;
 * bursts that are shorter than this duration are likely to be charachterized as noise and
 * ignored.
 *
 * Burst tags can be applied before / after the detected start / end of burst using the
 * `burst_pre_len` and `burst_post_len` options respectively. End of burst tags are
 * appended once the value of the center burst bin or either adjacent bin is gone for
 * `burst_post_len` FFTs, so it is recommended to set this value to some nonzero number to
 * reduce the chance a burst end is tagged too early. If this value is set too long,
 * multiple bursts may be combined into one burst.
 *
 * Debug information including detected peak value files (saved in /tmp) and a complex PDU
 * suitable for display on the in-tree QT Time Sink of FFTs and dynamic threshold can be
 * endabled for debugging, but it is not recommended for normal operation as these take a
 * non-trivial amount of compute cycles. Deeper level debug information can be enabled
 * using compile options in the source.
 *
 * \ingroup fhss_utils
 *
 */
class FHSS_UTILS_API fft_burst_tagger : virtual public gr::block
{
public:
    typedef std::shared_ptr<fft_burst_tagger> sptr;

    /*!
     * \brief Creates a new instance of fhss_utils::fft_burst_tagger.
     *
     * @param center_freq - center frequency of data stream, unit Hz
     * @param fft_size - number of bins in the primary FFT
     * @param sample_rate - sample rate of incoming stream
     * @param burst_pre_len - number of FFTs before the burst to place the START tag
     * @param burst_post_len - number of FFTs after the burst to place the END tag
     * @param burst_width - estimated bandwidth of bursts in Hz
     * @param max_bursts - maximum number of bursts allowed simultaneously (0=default)
     * @param max_burst_len - bursts exceeding this length will be tagged immediately
     * @param threshold - detection threshold above dynamic noise average (dB)
     * @param history_size - number of FFTs to compute noise estimate over
     * @param lookahead - number of FFTs a burst must be present to be tagged
     * @param debug - true enables debug file and message output functionality
     * @return shared_ptr<fft_burst_tagger>
     */
    static sptr make(float center_freq,
                     int fft_size,
                     int sample_rate,
                     int burst_pre_len,
                     int burst_post_len,
                     int burst_width,
                     int max_bursts = 0,
                     int max_burst_len = 0,
                     float threshold = 7,
                     int history_size = 512,
                     int lookahead = 10,
                     bool debug = false);

    /**
     * Returns total number of bursts seen
     *
     * @return uint64_t - number of bursts
     */
    virtual uint64_t get_n_tagged_bursts() = 0;

    /**
     * Resets burst tagger
     */
    virtual void reset() = 0;

    /**
     * Sets threshold
     * Unit: dB
     *
     * @param threshold - threshold
     */
    virtual void set_threshold(float threshold) = 0;

    /**
     * Sets max burst bandwidth
     *
     * @param bw - bandwidth in Hz
     */
    virtual void set_max_burst_bandwidth(double bw) = 0;

    /**
     * Preloads the noise floor estimate with a uniform value so that bursts can be
     * detected immediately
     *
     * @param noise_density - noise density to preload
     * @param preload - must be set to true to take effect
     */
    virtual void preload_noise_floor(double noise_density, bool preload = false) = 0;
}; // end class fft_burst_tagger


} // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_H */
