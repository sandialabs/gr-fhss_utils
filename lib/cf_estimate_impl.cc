/* -*- c++ -*- */
/*
 * Copyright 2020 gr-fhss_utils author.
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cf_estimate_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace fhss_utils {

cf_estimate::sptr cf_estimate::make(int method, std::vector<float> channel_freqs)
{
    return gnuradio::get_initial_sptr(new cf_estimate_impl(method, channel_freqs));
}

/*
 * The private constructor
 */
cf_estimate_impl::cf_estimate_impl(int method, std::vector<float> channel_freqs)
    : gr::block("cf_estimate",
                gr::io_signature::make(0, 0, 0),
                gr::io_signature::make(0, 0, 0)),
      d_method(method),
      d_channel_freqs(channel_freqs),
      d_gauss_sigma(7.0f / 8.0f),
      d_mags2(nullptr)
{
    message_port_register_in(PMT_IN);
    message_port_register_out(PMT_OUT);
    message_port_register_out(PMT_DEBUG);
    set_msg_handler(PMT_IN, boost::bind(&cf_estimate_impl::pdu_handler, this, _1));
    fft_setup(15);
}

/*
 * Our virtual destructor.
 */
cf_estimate_impl::~cf_estimate_impl() { fft_cleanup(); }

//////////////////////////////////////////////
// helper functions
//////////////////////////////////////////////

void cf_estimate_impl::fft_setup(int power)
{
    // initialize every power of fft up to "power" if they don't already exist
    for (int i = d_ffts.size(); i <= power; i++) {

        // init fft
        int fftsize = pow(2, i);
        d_ffts.push_back(new gr::fft::fft_complex(fftsize, true, 1));

        // init window
        d_windows.push_back(
            (float*)volk_malloc(sizeof(float) * fftsize, volk_get_alignment()));
        float two_sigma_squared = fftsize * d_gauss_sigma;
        two_sigma_squared *= 2 * two_sigma_squared;
        for (int j = 0; j < fftsize; j++) {
            float x = (-fftsize + 1) / 2.0f + j;
            d_windows[i][j] = std::exp((-x * x) / two_sigma_squared);
            // d_windows[i][j] = 1;
        }

        // init d_mags2
        if (d_mags2 != nullptr) {
            volk_free(d_mags2);
        }
        // VOLK 2.1 UPDATE ME
        d_mags2 = (float*)volk_malloc(sizeof(float) * fftsize, volk_get_alignment());
    }
}

void cf_estimate_impl::fft_cleanup()
{
    for (gr::fft::fft_complex* fft : d_ffts) {
        delete fft;
    }
    d_ffts.clear();
    for (float* win : d_windows) {
        volk_free(win);
    }
    d_windows.clear();
    if (d_mags2 != nullptr) {
        volk_free(d_mags2);
    }
    d_mags2 = nullptr;
}


//////////////////////////////////////////////
// message handlers functions
//////////////////////////////////////////////
void cf_estimate_impl::pdu_handler(pmt::pmt_t pdu)
{
    //////////////////////////////////
    // basic checks and pdu parsing
    //////////////////////////////////
    if (!pmt::is_pair(pdu)) {
        GR_LOG_WARN(d_logger, "PDU is not a pair, dropping\n");
        return;
    }

    pmt::pmt_t metadata = pmt::car(pdu);
    pmt::pmt_t pdu_data = pmt::cdr(pdu);

    if (!pmt::is_c32vector(pdu_data)) {
        GR_LOG_WARN(d_logger, "PDU is not c32vector, dropping\n");
        return;
    }

    if (!pmt::is_dict(metadata)) {
        GR_LOG_WARN(d_logger, "PDU metadata is not dict, dropping\n");
        return;
    }


    //////////////////////////////////
    // extract all needed data and metadata (sample_rate, center_frequency), optional
    // (relative_frequency)
    //////////////////////////////////
    if (!pmt::dict_has_key(metadata, PMT_CENTER_FREQUENCY) ||
        !pmt::dict_has_key(metadata, PMT_SAMPLE_RATE)) {
        GR_LOG_WARN(d_logger,
                    "cf_estimate needs 'center_frequency' and 'sample_rate' metadata\n");
        return;
    }
    double center_frequency =
        pmt::to_double(pmt::dict_ref(metadata, PMT_CENTER_FREQUENCY, pmt::PMT_NIL));
    double relative_frequency = pmt::to_double(
        pmt::dict_ref(metadata, PMT_RELATIVE_FREQUENCY, pmt::from_double(0.0)));
    double sample_rate =
        pmt::to_double(pmt::dict_ref(metadata, PMT_SAMPLE_RATE, pmt::PMT_NIL));

    // extract the data portion
    size_t burst_size;
    const gr_complex* data = pmt::c32vector_elements(pdu_data, burst_size);


    //////////////////////////////////
    // frequency analysis & PSD estimate
    //////////////////////////////////
    // fftsize is the data burst_size rounded up to a power of 2
    bool use_ceil = false;
    int fftpower = 0;
    size_t copy_size = 0;
    if (use_ceil) {
        fftpower = std::ceil(log2(burst_size));
        copy_size = burst_size;
    } else {
        fftpower = std::floor(log2(burst_size));
        copy_size = pow(2, fftpower);
    }
    int fftsize = pow(2, fftpower);
    fft_setup(fftpower);

    // copy in data to the right buffer, pad with zeros
    gr_complex* fft_in = d_ffts[fftpower]->get_inbuf();
    memset(fft_in, 0, sizeof(gr_complex) * fftsize);
    memcpy(fft_in, data, sizeof(gr_complex) * copy_size);

    // apply gaussian window in place
    volk_32fc_32f_multiply_32fc(fft_in, fft_in, d_windows[fftpower], copy_size);

    // perform the fft
    d_ffts[fftpower]->execute();

    // get magnitudes squared from fft output
    gr_complex* fft_out = d_ffts[fftpower]->get_outbuf();
    volk_32fc_magnitude_squared_32f(d_mags2, fft_out, fftsize);

    // fft shift and copy into STL object
    std::vector<float> mags2(d_mags2, d_mags2 + fftsize);
    std::rotate(mags2.begin(), mags2.begin() + fftsize / 2, mags2.end());

    // build the frequency axis, centered at center_frequency, in Hz
    double step_size = sample_rate / (float)fftsize;
    double start = center_frequency - (sample_rate / 2.0);
    std::vector<float> freq_axis;
    freq_axis.reserve(fftsize);
    for (size_t i = 0; i < fftsize; i++) {
        freq_axis.push_back(start + step_size * i);
    }

    // debug port publishes PSD
    message_port_pub(PMT_DEBUG,
                     pmt::cons(metadata, pmt::init_f32vector(fftsize, &mags2[0])));


    //////////////////////////////////
    // center frequency estimation
    //////////////////////////////////
    double shift = 0.0;
    // these methods return 'shift', a frequency correction from [-0.5,0.5]
    if (d_method == COERCE) {
        shift = this->coerce_frequency(center_frequency, sample_rate);
    } else if (d_method == RMS) {
        shift = this->rms(mags2, freq_axis, center_frequency, sample_rate);
    } else if (d_method == HALF_POWER) {
        shift = this->half_power(mags2);
    } else {
        GR_LOG_WARN(d_logger, "No valid method selected!\n");
        std::cout << " cf_estimate has no valid method selected! " << d_method
                  << std::endl;
    }

    //////////////////////////////////
    // correct the burst using the new center frequency
    //////////////////////////////////
    d_rotate.set_phase_incr(exp(gr_complex(0, -shift * 2.0 * M_PI)));
    d_rotate.set_phase(gr_complex(1, 0));
    d_corrected_burst.resize(burst_size);
    d_rotate.rotateN(&d_corrected_burst[0], data, burst_size);

    float cf_correction_hz = shift * sample_rate;
    center_frequency += cf_correction_hz;
    if (relative_frequency != 0)
        relative_frequency += cf_correction_hz;

    //////////////////////////////////
    // estimate Bandwidth and SNR
    //////////////////////////////////
    float bandwidth = rms_bw(mags2, freq_axis, center_frequency);
    float snr_db =
        snr_estimation(mags2, freq_axis, center_frequency, bandwidth, sample_rate);

    //////////////////////////////////
    // build the output pdu
    //////////////////////////////////
    metadata =
        pmt::dict_add(metadata, PMT_CENTER_FREQUENCY, pmt::from_double(center_frequency));
    if (relative_frequency != 0)
        metadata = pmt::dict_add(
            metadata, PMT_RELATIVE_FREQUENCY, pmt::from_double(relative_frequency));
    metadata = pmt::dict_add(metadata, PMT_BANDWIDTH, pmt::from_double(bandwidth));
    metadata = pmt::dict_add(metadata, PMT_SNRDB, pmt::from_double(snr_db));

    pmt::pmt_t out_data = pmt::init_c32vector(burst_size, d_corrected_burst);
    message_port_pub(PMT_OUT, pmt::cons(metadata, out_data));
}


//////////////////////////////////////////////
// center frequency estimation methods
//////////////////////////////////////////////
float cf_estimate_impl::half_power(std::vector<float> mags2)
{
    // calculate the total energy
    double energy = 0.0;
    for (auto& p : mags2) {
        energy += p;
    }

    // find center
    double half_power = energy / 2.0;
    double running_total = 0.0;
    size_t half_power_idx = 0;
    for (size_t i = 0; running_total < half_power; i++) {
        running_total += mags2[i];
        half_power_idx = i;
    }

    // convert index to frequency and return the correction/shift
    return ((double)half_power_idx / (double)mags2.size()) - 0.5;
}

float cf_estimate_impl::rms(std::vector<float> mags2,
                            std::vector<float> freq_axis,
                            float center_frequency,
                            float sample_rate)
{
    // python: cf_estimate = sum( [ f*(p**2) for f,p in zip(freq_axis, freq_mag)] ) /
    // energy
    double energy = 0.0;
    for (auto& p : mags2) {
        energy += p;
    }

    // calculate the top integral: integrate(f*PSD(f)^2, df)
    double top_integral = 0.0;
    for (size_t i = 0; i < mags2.size(); i++) {
        top_integral += freq_axis[i] * mags2[i];
    }

    // normalize by total energy
    return ((top_integral / energy) - center_frequency) / sample_rate;
}

float cf_estimate_impl::coerce_frequency(float center_frequency, float sample_rate)
{
    // we can only coerce the frequencies if we have a list of good frequencies
    if (d_channel_freqs.size() == 0) {
        GR_LOG_WARN(d_logger, "No channel freqs to coerce!\n");
        return 0.0;
    }

    // find the closest listed frequency to this center_frequency
    float dist = abs(d_channel_freqs[0] - center_frequency);
    float channel_freq = d_channel_freqs[0];
    for (size_t i = 1; i < d_channel_freqs.size(); ++i) {
        float temp_dist = abs(d_channel_freqs[i] - center_frequency);
        if (temp_dist < dist) {
            dist = temp_dist;
            channel_freq = d_channel_freqs[i];
        }
    }

    // return the shift or correction amount
    return (channel_freq - center_frequency) / sample_rate;

} /* end coerce_frequency */


//////////////////////////////////////////////
// BW and SNR estimation methods
//////////////////////////////////////////////
float cf_estimate_impl::rms_bw(std::vector<float> mags2,
                               std::vector<float> freq_axis,
                               float center_frequency)
{
    // python: bw_rms = np.sqrt( sum( [(abs(f-cf_estimate)**2) * (p**2) for f,p in
    // zip(freq_axis, freq_mag)] ) / energy)
    double energy = 0.0;
    for (auto& p : mags2) {
        energy += p;
    }

    // calculate the top integral, integrate((f-cf)^2 * PSD(f)^2, df)
    double top_integral = 0.0;
    for (size_t i = 0; i < mags2.size(); i++) {
        top_integral += (std::pow(freq_axis[i] - center_frequency, 2) * mags2[i]);
    }

    // normalize by total energy and take sqrt
    return std::sqrt(top_integral / energy);
}

float cf_estimate_impl::snr_estimation(std::vector<float> mags2,
                                       std::vector<float> freq_axis,
                                       float center_frequency,
                                       float bandwidth,
                                       float sample_rate)
{
    // average burst power - center_frequency to +- bandwidth/2
    double start_freq = center_frequency - bandwidth / 2.0;
    double stop_freq = center_frequency + bandwidth / 2.0;
    double signal_power = 0.0;
    double noise_power = 0.0;
    for (size_t i = 0; i < mags2.size(); i++) {
        if ((freq_axis[i] > start_freq) && (freq_axis[i] < stop_freq)) {
            signal_power += mags2[i];
        } else {
            noise_power += mags2[i];
        }
    }

    // normalize, ratio, and convert to dB
    signal_power /= bandwidth;
    noise_power /= (sample_rate - bandwidth);
    return 10 * std::log10(signal_power / noise_power);
}

//////////////////////////////////////////////
// getters/setters
//////////////////////////////////////////////
void cf_estimate_impl::set_freqs(std::vector<float> channel_freqs)
{
    d_channel_freqs = channel_freqs;
}

void cf_estimate_impl::set_method(int method) { d_method = method; }

} /* namespace fhss_utils */
} /* namespace gr */
