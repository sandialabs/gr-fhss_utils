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


#ifndef INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_H
#define INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_H

#include <gnuradio/sync_block.h>
#include <fhss_utils/api.h>

namespace gr {
namespace fhss_utils {

/*!
 * \brief Tags Detected bursts
 * \ingroup fhss_utils
 *
 */
class FHSS_UTILS_API fft_burst_tagger : virtual public gr::sync_block
{
public:
    typedef boost::shared_ptr<fft_burst_tagger> sptr;

    /*!
     * \brief Creates a new instance of fhss_utils::fft_burst_tagger.
     *
     * @param center_frequency - center frequency of incoming stream, unit Hz
     * @param fft_size -
     * @param sample_rate - sample rate of incoming stream
     * @param burst_pre_len - XXX unit seconds
     * @param burst_post_len - XXX unit seconds
     * @param burst_width -XXX unit Hz
     * @param max_bursts - XXX
     * @param max_burst_len - XXX
     * @param threshold - XXX unit dB
     * @param history_size - XXX
     * @param lookahead - XXX unit seconds
     * @param debug - true enables debug functionality
     * @return shared_ptr<fft_burst_tagger>
     */
    static sptr make(float center_frequency,
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
     * Sets max burst bandwidth
     * Unit: Hz
     *
     * @param bw - bandwidth
     */
    virtual void set_max_burst_bandwidth(double bw) = 0;
}; // end class fft_burst_tagger


} // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_FFT_BURST_TAGGER_H */
