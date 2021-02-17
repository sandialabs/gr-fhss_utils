/* -*- c++ -*- */
/*
 * Copyright 2018-2021 National Technology & Engineering Solutions of Sandia, LLC
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
    std::vector<gr::fft::fft_complex_fwd*> d_ffts;
    std::vector<float*> d_windows;
    std::vector<float> d_fft_mag2_gains;
    float* d_mags2;
    float* d_magnitude_shifted_f;

    // correction tools
    blocks::rotator d_rotate;
    std::vector<gr_complex> d_corrected_burst;

    // coerce tools
    std::vector<float> d_channel_freqs;

    /*!
     * \brief Coerce Frequency to nearest value in d_channel_freqs
     *
     * \param center_frequency Center Frequency of input data
     * \param sample_rate Sample Rate of input data
     */
    float coerce_frequency(float center_frequency, float sample_rate);

    /*!
     * \brief Return center frequency estimate using the RMS method
     *
     * \param mags2 Vector of magnitude^2 FFT of input data
     * \param freq_axis Frequency of each bin in the mags2 vector
     * \param center_frequency Center Frequency of input data
     * \param sample_rate Sample Rate of input data
     */
    float rms(std::vector<float> mags2,
              std::vector<float> freq_axis,
              float center_frequency,
              float sample_rate);

    /*!
     * \brief Return center frequency estimate using the Half Power method
     *
     * \param mags2 Vector of magnitude^2 FFT of input data
     */
    float half_power(std::vector<float> mags2);

    // BW and SNR estimation

    /*!
     * \brief Return bandwidth estimate using the RMS method
     *
     * \param mags2 Vector of magnitude^2 FFT of input data
     * \param freq_axis Frequency of each bin in the mags2 vector
     * \param center_frequency Center Frequency of input data
     */
    float rms_bw(std::vector<float> mags2,
                 std::vector<float> freq_axis,
                 float center_frequency);

    /*!
     * \brief Return estimated power in a burst
     *
     * \param mags2 Vector of magnitude^2 FFT of input data
     * \param freq_axis Frequency of each bin in the mags2 vector
     * \param center_frequency Center Frequency of input data
     * \param bandwidth Bandwidth of signal
     * \param fft_mag2_gain gain of FFT and windowing process
     */
    float estimate_pwr(std::vector<float> mags2,
                       std::vector<float> freq_axis,
                       float center_frequency,
                       float bandwidth,
                       float fft_mag2_gain);

public:
    cf_estimate_impl(int method, std::vector<float> channel_freqs);

    ~cf_estimate_impl() override;

    void set_freqs(std::vector<float> channel_freqs) override
    {
        d_channel_freqs = channel_freqs;
    };
    void set_method(int method) override { d_method = method; };
};

} // namespace fhss_utils
} // namespace gr

#endif /* INCLUDED_FHSS_UTILS_CF_ESTIMATE_IMPL_H */
