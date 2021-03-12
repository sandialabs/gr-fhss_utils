/* -*- c++ -*- */
/*
 * Copyright 2018, 2019, 2020 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 * Copyright 2021 Jacob Gilbert
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>
#include "cf_estimate_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace fhss_utils {

const int MAX_FFT_POWER = 8;    // largest FFT this block will use
const int MIN_FFT_POWER = 5;    // smallest FFT this block will use
const int MIN_NFFTS = 4;        // minimum number of FFTs this block will average
const int MIN_BURST_SIZE = pow(2, MIN_FFT_POWER) * MIN_NFFTS;
// PDUs should be large enough to take MIN_NFFTS of size pow(2,MIN_FFT_POWER)


/*
 *    TODO    there is WAY too much use of copy constructors for std::vectors
 *    TODO    in the various helper functions implemneted in this block...
 */


cf_estimate::sptr cf_estimate::make(int method, std::vector<float> channel_freqs, float snr_min)
{
    return gnuradio::get_initial_sptr(new cf_estimate_impl(method, channel_freqs, snr_min));
}

/*
 * The private constructor
 */
cf_estimate_impl::cf_estimate_impl(int method, std::vector<float> channel_freqs, float snr_min)
    : gr::block("cf_estimate",
                gr::io_signature::make(0, 0, 0),
                gr::io_signature::make(0, 0, 0)),
      d_method(method),
      d_channel_freqs(channel_freqs),
      d_snr_min(snr_min)
{
    message_port_register_in(PMTCONSTSTR__in());
    message_port_register_out(PMTCONSTSTR__out());
    message_port_register_out(PMTCONSTSTR__debug());
    set_msg_handler(PMTCONSTSTR__in(),
                    [this](pmt::pmt_t msg) { this->cf_estimate_impl::pdu_handler(msg); });

    fft_setup(MAX_FFT_POWER);

    if (d_channel_freqs.size() == 0 && d_method == COERCE) {
        GR_LOG_WARN(d_logger,
                    "CF Estimator operating in COERCE mode with empty channel frequency "
                    "list; no CF correction will be applied!");
    }
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
    for (int i = 0; i <= power; i++) {

        /* initialize fft */
        int fftsize = pow(2, i);
        d_ffts.push_back(new gr::fft::fft_complex(fftsize, true, 1));

        /*
         * initialize windows and calculate window gain, which in this case specifically
         * refers to the non-coherent gain which is the rms value of the window weights:
         *
         *   G_{nc} = \sqrt{\frac{\Sigma w[n]^2}{N}}
         *
         * FIXME GR3.9: use the new GAUSSIAN/TUKEY built in fft::window kernels to
         * generate the window
         */
        bool USE_GAUSSIAN_WINDOW(true);

        d_windows.push_back(
            (float*)volk_malloc(sizeof(float) * fftsize, volk_get_alignment()));

        double g = 0;
        // this gaussian window is narrow to create a well localized spectral peak
        if (USE_GAUSSIAN_WINDOW) {
            float two_sigma_squared = fftsize * 1.0 / 32.0; // former d_guass_sigma = 7/8
            two_sigma_squared *= 2 * two_sigma_squared;
            for (int j = 0; j < fftsize; j++) {
                float x = (-fftsize + 1) / 2.0f + j;
                d_windows[i][j] = std::exp((-x * x) / two_sigma_squared);
                g += (d_windows[i][j] * d_windows[i][j]);
            }
        } else {
            std::vector<float> window =
                fft::window::build(fft::window::WIN_BLACKMAN, fftsize, 0);
            for (auto j = 0; j < fftsize; j++) {
                d_windows[i][j] = window[j];
                g += (d_windows[i][j] * d_windows[i][j]);
            }
        }
        // scale the window to compensate for FFT size and window rms gain:
        //   gain_rms^2 * fftsize = sqrt(win_gain / fftsize)^2 * fftsize = win_gain
        for (auto j = 0; j < fftsize; j++) {
            d_windows[i][j] /= g*4;     // FIXME: magic number so power/noise are correct
        }
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
    if (!pmt::dict_has_key(metadata, PMTCONSTSTR__center_frequency()) ||
        !pmt::dict_has_key(metadata, PMTCONSTSTR__sample_rate())) {
        GR_LOG_WARN(d_logger,
                    "cf_estimate needs 'center_frequency' and 'sample_rate' metadata\n");
        return;
    }
    double center_frequency = pmt::to_double(
        pmt::dict_ref(metadata, PMTCONSTSTR__center_frequency(), pmt::PMT_NIL));
    double relative_frequency = pmt::to_double(pmt::dict_ref(
        metadata, PMTCONSTSTR__relative_frequency(), pmt::from_double(0.0)));
    double sample_rate =
        pmt::to_double(pmt::dict_ref(metadata, PMTCONSTSTR__sample_rate(), pmt::PMT_NIL));
    double noise_density_db = pmt::to_double(
        pmt::dict_ref(metadata, PMTCONSTSTR__noise_density(), pmt::from_double(NAN)));

    // extract the data portion
    size_t burst_size;
    const gr_complex* data = pmt::c32vector_elements(pdu_data, burst_size);

    if (burst_size < MIN_BURST_SIZE) {
        GR_LOG_INFO(d_logger,
                    boost::format("Burst of length %d too small (min = %d), dropping.") %
                        burst_size % MIN_BURST_SIZE);
    }

    //////////////////////////////////
    // frequency analysis & PSD estimate
    //////////////////////////////////
    // fftsize is the data burst_size rounded down to a power of 2
    int fftpower = std::floor(log2(burst_size * 1.0 / MIN_NFFTS));

    if (fftpower > MAX_FFT_POWER) {
        fftpower = MAX_FFT_POWER;
    }

    size_t fftsize = pow(2, fftpower);
    size_t nffts = std::floor(1.0 * burst_size / fftsize);
    size_t copy_size = nffts * fftsize;

    // the offset ensures the center part of the burst is used
    uint32_t offset = std::max((size_t)0, (burst_size - copy_size) / 2);

    std::vector<float> fftm2(fftsize, 0); // stores the result of each fft
    std::vector<float> mags2(fftsize, 0); // stores the accumulated result of each fft

    for (int ii = 0; ii < nffts; ii++) {

        // TODO: verify VOLK alignment for these calls below....doesnt look like it

        // copy in data to the right buffer, pad with zeros
        gr_complex* fft_in = d_ffts[fftpower]->get_inbuf();
        memset(fft_in, 0, sizeof(gr_complex) * fftsize);
        memcpy(fft_in, data + offset + ii * fftsize, sizeof(gr_complex) * fftsize);

        // apply gain-compensated window in place and get mag^2 FFT output
        volk_32fc_32f_multiply_32fc(fft_in, fft_in, d_windows[fftpower], fftsize);
        d_ffts[fftpower]->execute();
        volk_32fc_magnitude_squared_32f(
            &fftm2[0], d_ffts[fftpower]->get_outbuf(), fftsize);

        // accumulate
        //std::transform(mags2.begin(), mags2.end(), fftm2.begin(), mags2.begin(), std::plus<float>());
        for (int jj=0;jj<fftsize;jj++) mags2[jj] += fftm2[jj];
    }

    // normalize accumulated fft bins
    //volk_32f_s32f_multiply_32f(&mags2[0], &mags2[0], 1.0 / nffts, fftsize);
    for (int ii=0;ii<fftsize;ii++) mags2[ii] /= nffts;

    // fft shift
    std::rotate(mags2.begin(), mags2.begin() + fftsize / 2, mags2.end());

    // if the noise floor was not in the metadata, set it to the smallest value
    double bin_size = sample_rate / (float)fftsize;
    float noise_floor_db;
    if (std::isnan(noise_density_db)) {
        noise_floor_db = 10 * log10(*std::min_element(std::begin(mags2), std::end(mags2)));
    } else {
        noise_floor_db = noise_density_db + 10 * log10(bin_size);
    }

    // build the frequency axis, centered at center_frequency, in Hz
    double start = center_frequency - (sample_rate / 2.0);
    std::vector<float> freq_axis;
    freq_axis.reserve(fftsize);
    for (size_t i = 0; i < fftsize; i++) {
        freq_axis.push_back(start + bin_size * i);
    }

        // generate the debug message TODO: make conditional
    d_bug.clear();
    d_bug.reserve(fftsize);
    for (auto ii = 0; ii < fftsize; ii++)
        d_bug.push_back(gr_complex(10 * log10(mags2[ii]), noise_floor_db));

    //////////////////////////////////
    // center frequency and bandwidth estimation
    //////////////////////////////////
    float bandwidth = rms_bw(mags2, freq_axis, center_frequency);

    double shift = 0.0;
    // these methods return 'shift', a frequency correction from [-0.5,0.5]
    if (d_method == RMS) {
        shift = this->rms(mags2, freq_axis, center_frequency, sample_rate);
    } else if (d_method == HALF_POWER) {
        shift = this->half_power(mags2);
    } else if (d_method == MIDDLE_OUT) {
        // this is going to update the bandwdith value too
        shift = this->middle_out(mags2, freq_axis, noise_floor_db, bandwidth);
    } else {
        // in COERCE mode no estimation is performed
    }

    // if a frequency coercion list has been provided, apply that
    shift +=
        this->coerce_frequency((center_frequency + (shift * sample_rate)), sample_rate);

    // debug port publishes PSD
    message_port_pub(PMTCONSTSTR__debug(),
                     pmt::cons(metadata, pmt::init_c32vector(fftsize, &d_bug[0])));


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
    // estimate SNR and build the output PDU
    //////////////////////////////////
    metadata = pmt::dict_add(
        metadata, PMTCONSTSTR__center_frequency(), pmt::from_double(center_frequency));
    if (relative_frequency != 0)
        metadata = pmt::dict_add(metadata,
                                 PMTCONSTSTR__relative_frequency(),
                                 pmt::from_double(relative_frequency));

    metadata =
        pmt::dict_add(metadata, PMTCONSTSTR__bandwidth(), pmt::from_double(bandwidth));

    if (noise_density_db != NAN) {
        float pwr_db = estimate_pwr(mags2, freq_axis, center_frequency, bandwidth);
        float snr_db = pwr_db - (noise_density_db + 10 * log10(bandwidth));
        metadata =
            pmt::dict_add(metadata, PMTCONSTSTR__pwr_db(), pmt::from_double(pwr_db));
        metadata =
            pmt::dict_add(metadata, PMTCONSTSTR__snr_db(), pmt::from_double(snr_db));
    }

    pmt::pmt_t out_data = pmt::init_c32vector(burst_size, d_corrected_burst);
    message_port_pub(PMTCONSTSTR__out(), pmt::cons(metadata, out_data));
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

/*
 * this computes the bandwidth as defined by the portion of the signal that is
 * continuously greater than halfway betwen the peak-8dB and the noise floor
 */
float cf_estimate_impl::middle_out(std::vector<float> mags2,
                                   std::vector<float> freq_axis,
                                   float noise_floor_db,
                                   float &bandwidth)
{
    auto peak = std::max_element(std::begin(mags2), std::end(mags2));
    int peak_idx = std::distance(std::begin(mags2), peak);

    // determine the bandwidth threshold; linear-scale mag^2
    float threshold = (10 * log10(*peak) - d_snr_min + noise_floor_db) / 2.0;
    threshold = pow(10, threshold / 10.0);
    //threshold = pow(10, (noise_floor_db-0.1) / 10.0);

    // now determine the first index below the threshold in either direction
    auto high(peak);
    auto low(peak);
    for (; high < mags2.end()-1; high++)
        if (*high <= threshold) break;
    for (; low > mags2.begin(); low--)
        if (*low <= threshold) break;
    int high_idx = peak_idx + std::distance(peak, high);
    int low_idx = peak_idx + std::distance(peak, low);

    bandwidth = freq_axis[high_idx] - freq_axis[low_idx];
    float shift = (high_idx + low_idx) *0.5 - (mags2.size() - 1) * 0.5;
    shift /= mags2.size();

    for (int i = low_idx; i <= high_idx; i++)
        d_bug[i] = std::real(d_bug[i]) + 1j * 10*log10(threshold);
    d_bug[peak_idx] = std::real(d_bug[peak_idx]) + 1j * 10*log10(*peak);

    return shift;
}

/*
 * FIXME: using the frequency axis is an odd way to do this... in general the logic here
 * could use some cleanup
 */
float cf_estimate_impl::estimate_pwr(std::vector<float> mags2,
                                     std::vector<float> freq_axis,
                                     float center_frequency,
                                     float bandwidth)
{
    double start_freq = center_frequency - bandwidth / 2.0;
    double stop_freq = center_frequency + bandwidth / 2.0;
    double signal_power = 0.0;
    for (auto i = 0; i < mags2.size(); i++) {
        if ((freq_axis[i] > start_freq) && (freq_axis[i] < stop_freq)) {
            signal_power += mags2[i];
        }
    }

    // return the total signal power
    return (10 * log10(signal_power));
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
