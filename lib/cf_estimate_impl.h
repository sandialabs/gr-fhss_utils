/* -*- c++ -*- */
/*
 * Copyright 2018, 2019, 2020 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_FHSS_UTILS_CF_ESTIMATE_IMPL_H
#define INCLUDED_FHSS_UTILS_CF_ESTIMATE_IMPL_H

#include <gnuradio/blocks/rotator.h>
#include <gnuradio/fft/fft.h>
#include <fhss_utils/cf_estimate.h>
#include <fhss_utils/constants.h>

namespace gr {
namespace fhss_utils {

class cf_estimate_impl : public cf_estimate
{
private:
    // message handlers
    void pdu_handler(pmt::pmt_t pdu);

    // this block supports a few different estimation methods
    int d_method;

    // fft tools
    void fft_setup(int power);
    void fft_cleanup();
    std::vector<gr::fft::fft_complex*> d_ffts;
    std::vector<float*> d_windows;
    std::vector<float> d_fft_mag2_gains;
    float* d_mags2;
    float* d_magnitude_shifted_f;

    // correction tools
    blocks::rotator d_rotate;
    std::vector<gr_complex> d_corrected_burst;

    // coerce tools
    std::vector<float> d_channel_freqs;
    float coerce_frequency(float center_frequency, float sample_rate);
    float rms(std::vector<float> mags2,
              std::vector<float> freq_axis,
              float center_frequency,
              float sample_rate);
    float half_power(std::vector<float> mags2);

    // BW and SNR estimation
    float rms_bw(std::vector<float> mags2,
                 std::vector<float> freq_axis,
                 float center_frequency);
    float estimate_pwr(std::vector<float> mags2,
                       std::vector<float> freq_axis,
                       float center_frequency,
                       float bandwidth,
                       float fft_mag2_gain);

public:
    /**
     * Constructor
     *
     * @param method - estimate method #cf_method
     * @param channel_freqs - channel freq list for coerce method
     */
    cf_estimate_impl(int method, std::vector<float> channel_freqs);

    /**
     * Deconstructor
     */
    ~cf_estimate_impl();

    // getters and setters

    /**
     * set coerce channel frequency list
     *
     * @param channel_freqs - channel freq list for coerce method
     */
    void set_freqs(std::vector<float> channel_freqs);

    /**
     * Set estimate method
     *
     * @param method - estimate method #cf_method
     */
    void set_method(int method);
};

} // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_CF_ESTIMATE_IMPL_H */
