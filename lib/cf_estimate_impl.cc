/* -*- c++ -*- */
/*
 * Copyright 2018-2021 National Technology & Engineering Solutions of Sandia, LLC
 * (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
 * retains certain rights in this software.
 * Copyright 2021 Jacob Gilbert
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cf_estimate_impl.h"
#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>
#include <gnuradio/fhss_utils/constants.h>
#include <gnuradio/io_signature.h>

namespace gr {
namespace fhss_utils {

const int MAX_FFT_POWER = 8;    // largest FFT this block will use
const int MIN_FFT_POWER = 5;    // smallest FFT this block will use
const int MIN_NFFTS = 4;        // minimum number of FFTs this block will average
const int MIN_BURST_SIZE = pow(2, MIN_FFT_POWER) * MIN_NFFTS;
// PDUs should be large enough to take MIN_NFFTS of size pow(2,MIN_FFT_POWER)


cf_estimate::sptr cf_estimate::make(int method, std::vector<float> channel_freqs)
{
    return gnuradio::make_block_sptr<cf_estimate_impl>(method, channel_freqs);
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
      d_snr_min(10.0),
      d_thresh_min(-25.0)
{
    message_port_register_in(PMTCONSTSTR__in());
    message_port_register_out(PMTCONSTSTR__out());
    message_port_register_out(PMTCONSTSTR__debug());
    set_msg_handler(PMTCONSTSTR__in(),
                    [this](pmt::pmt_t msg) { this->pdu_handler(msg); });
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
        d_ffts.push_back(new gr::fft::fft_complex_fwd(fftsize, 1));

        /* initialize normalized FFT windows and pre-compensate for FFT size */
        d_windows.push_back(
            (float*)volk_malloc(sizeof(float) * fftsize, volk_get_alignment()));

        // this gaussian window is narrow and produces a well defined spectral peak
        float sigma = fftsize * 1.0 / 32.0;
        std::vector<float> win(fftsize, 0);
        for (int j = 0; j < fftsize; j++) {
            float x = (-fftsize + 1) / 2.0f + j;
            win[j] = std::exp((-x * x) / (2 * sigma * sigma));
        }
        // scale the window to compensate for FFT size and window rms gain:
        double pwr_acc = 0.0;
        for (auto x: win) pwr_acc += x*x/fftsize;
        for (auto j=0; j < fftsize; j++) {
            d_windows[i][j] = win[j] / (std::sqrt(pwr_acc) * fftsize);
        }
        // std::vector<float> window =
        //     fft::window::build(fft::window::WIN_GAUSSIAN, fftsize, sigma, true);
        // //    fft::window::build(fft::window::WIN_BLACKMAN, fftsize, 0, true);
        // for (int j = 0; j < fftsize; j++) {
        //     d_windows[i][j] = window[j] / fftsize;
        // }
    }
}

void cf_estimate_impl::fft_cleanup()
{
    for (gr::fft::fft_complex_fwd* fft : d_ffts) {
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
    if (!pmt::is_pdu(pdu)) {
        GR_LOG_WARN(d_logger, "Message is not a PDU, dropping\n");
        return;
    }

    pmt::pmt_t meta = pmt::car(pdu);
    pmt::pmt_t pdu_data = pmt::cdr(pdu);

    //////////////////////////////////
    // extract all needed data and metadata (sample_rate, center_frequency)
    //////////////////////////////////
    if (!pmt::dict_has_key(meta, PMTCONSTSTR__center_frequency()) ||
        !pmt::dict_has_key(meta, PMTCONSTSTR__sample_rate())) {
        GR_LOG_WARN(d_logger,
                    "cf_estimate needs 'center_frequency' and 'sample_rate' metadata\n");
        return;
    }
    double center_frequency = pmt::to_double(
        pmt::dict_ref(meta, PMTCONSTSTR__center_frequency(), pmt::PMT_NIL));
    double relative_frequency = pmt::to_double(pmt::dict_ref(
        meta, PMTCONSTSTR__relative_frequency(), pmt::from_double(0.0)));
    double sample_rate =
        pmt::to_double(pmt::dict_ref(meta, PMTCONSTSTR__sample_rate(), pmt::PMT_NIL));
    double noise_density_db = pmt::to_double(
        pmt::dict_ref(meta, PMTCONSTSTR__noise_density(), pmt::from_double(NAN)));
    // TODO: scale the coarse_bandwidth by a parameter? A fixed value? Not at all?
    double coarse_bandwidth = pmt::to_double(
        pmt::dict_ref(meta, PMTCONSTSTR__bandwidth(), pmt::from_double(sample_rate)));

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

    // use the initial bandwidth estimate to determine which FFT bins to process
    // for finer center frequency and bandwidth estimates
    size_t start_bin = 0;
    // if the coarse bandwidth is at least the sample rate, use the whole FFT
    if (coarse_bandwidth >= sample_rate) {
        start_bin = 0;
    // if the bandwidth is not positive, use the two middle bins
    } else if (coarse_bandwidth <= 0.0) {
        start_bin = (size_t) (bin_size / 2.0 - 1.0);
    }
    else {
        start_bin = (size_t)((sample_rate - coarse_bandwidth) / (2.0 * bin_size));
    }
    //while (center_frequency - freq_axis[start_bin + 1] > coarse_bandwidth / 2.0) {
    //    start_bin++;
    //}

    //////////////////////////////////
    // center frequency and bandwidth estimation
    //////////////////////////////////
    float shift = 0.0;     // frequency correction factor from [-0.5,0.5]
    float bandwidth = 0.0;  // bandwidth estimate in hz
    bool meas_err = false;  // indicates unreliable measurement

    if (d_method == MIDDLE_OUT) {
        meas_err |= middle_out(mags2, bin_size, noise_floor_db, bandwidth, shift, start_bin);
    } else {
        meas_err |= rms_bw(mags2, freq_axis, center_frequency, bandwidth, start_bin);
        if (d_method == RMS) {
            meas_err |= rms_cf(mags2, freq_axis, center_frequency, sample_rate, shift, start_bin);
        } else if (d_method == HALF_POWER) {
            meas_err |= half_power_cf(mags2, shift, start_bin);
        } else {
        // in COERCE mode no frequency estimation is performed
        }
    }

    // if a frequency coercion list has been provided, apply that
    float coerce_shift;
    meas_err |= coerce_frequency((center_frequency + (shift * sample_rate)), sample_rate, coerce_shift);
    shift += coerce_shift;

    // debug port publishes PSD
    message_port_pub(PMTCONSTSTR__debug(),
                     pmt::cons(meta, pmt::init_c32vector(fftsize, &d_bug[0])));

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
    meta = pmt::dict_add(
        meta, PMTCONSTSTR__center_frequency(), pmt::from_double(center_frequency));
    if (relative_frequency != 0)
        meta = pmt::dict_add(meta,
                                 PMTCONSTSTR__relative_frequency(),
                                 pmt::from_double(relative_frequency));

    meta = pmt::dict_add(meta, PMTCONSTSTR__bandwidth(), pmt::from_double(bandwidth));

    if (noise_density_db != NAN) {
        float pwr_db;
        meas_err |= estimate_pwr(mags2, freq_axis, center_frequency, bandwidth, pwr_db);
        float snr_db = pwr_db - (noise_density_db + 10 * log10(bandwidth));
        meta = pmt::dict_add(meta, PMTCONSTSTR__pwr_db(), pmt::from_double(pwr_db));
        meta = pmt::dict_add(meta, PMTCONSTSTR__snr_db(), pmt::from_double(snr_db));
    }
    if (meas_err) {
        meta = pmt::dict_add(meta, pmt::mp("meas_err"), pmt::PMT_T);
    }

    pmt::pmt_t out_data = pmt::init_c32vector(burst_size, d_corrected_burst);
    message_port_pub(PMTCONSTSTR__out(), pmt::cons(meta, out_data));
}


bool cf_estimate_impl::coerce_frequency(float center_frequency, float sample_rate, float &shift)
{
    // we can only coerce the frequencies if we have a list of good frequencies
    if (d_channel_freqs.size() == 0) {
        shift = 0.0;
        return false;
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
    shift = (channel_freq - center_frequency) / sample_rate;
    return false;
} /* end coerce_frequency */

//////////////////////////////////////////////
// center frequency estimation methods
//////////////////////////////////////////////
bool cf_estimate_impl::rms_cf(const std::vector<float> &mags2,
                               const std::vector<float> &freq_axis,
                               float center_frequency,
                               float sample_rate,
                               float &shift,
                               size_t start_bin)
{
    // python: cf_estimate = sum( [ f*(p**2) for f,p in zip(freq_axis, freq_mag)] ) /
    // energy
    double energy = 0.0;
    //for (auto& p : mags2) {
    //    energy += p;
    //}

    // calculate the top integral: integrate(f*PSD(f)^2, df)
    double top_integral = 0.0;
    for (size_t i = start_bin; i < mags2.size() - start_bin; i++) {
        energy += mags2[i];
        top_integral += freq_axis[i] * mags2[i];
    }

    // normalize by total energy
    shift = ((top_integral / energy) - center_frequency) / sample_rate;

    return false; // TODO: flag if error or low confidence
} /* end rms_cf */

bool cf_estimate_impl::half_power_cf(const std::vector<float> &mags2, float &shift, size_t start_bin)
{
    // calculate the total energy
    double energy = 0.0;
    for (size_t i = start_bin; i < mags2.size() - start_bin; i++) {
        energy += mags2[i];
    }
    //for (auto& p : mags2) {
    //    energy += p;
    //}

    // find center
    double half_power = energy / 2.0;
    double running_total = 0.0;
    size_t half_power_idx = 0;
    for (size_t i = start_bin; running_total < half_power; i++) {
        running_total += mags2[i];
        half_power_idx = i;
    }

    // convert index to frequency and return the correction/shift
    shift = ((double)half_power_idx / (double)mags2.size()) - 0.5;

    return false; // TODO: flag if error or low confidence
} /* end half_power_cf */

bool cf_estimate_impl::rms_bw(const std::vector<float> &mags2,
                               const std::vector<float> &freq_axis,
                               float center_frequency,
                               float &bandwidth,
                               size_t start_bin)
{
    // python: bw_rms = np.sqrt( sum( [(abs(f-cf_estimate)**2) * (p**2) for f,p in
    // zip(freq_axis, freq_mag)] ) / energy)

    //for (auto& p : mags2) {
    //    energy += p;
    //}

    // calculate the top integral, integrate((f-cf)^2 * PSD(f)^2, df)
    double energy = 0.0;
    double top_integral = 0.0;
    for (size_t i = start_bin; i < mags2.size() - start_bin; i++) {
        energy += mags2[i];
        top_integral += (std::pow(freq_axis[i] - center_frequency, 2) * mags2[i]);
    }

    // normalize by total energy and take sqrt
    bandwidth = std::sqrt(top_integral / energy);

    return false; // TODO: flag if error or low confidence
} /* end rms_bw */

/*
 * this computes the bandwidth as defined by the portion of the signal that is
 * continuously greater than halfway betwen the peak-8dB and the noise floor
 */
bool cf_estimate_impl::middle_out(const std::vector<float> &mags2,
                                   float bin_resolution,
                                   float noise_floor_db,
                                   float &bandwidth,
                                   float &shift,
                                   size_t start_bin)
{
    auto peak = std::max_element(std::begin(mags2) + start_bin, std::end(mags2) - start_bin);
    int peak_idx = std::distance(std::begin(mags2), peak);

    // determine the bandwidth threshold; linear-scale mag^2
    float peak_db = 10 * log10(*peak);
    float threshold = (peak_db - d_snr_min + noise_floor_db) / 2.0;
    if (threshold < noise_floor_db) {
        threshold = noise_floor_db;
    } else if (threshold < (peak_db + d_thresh_min)) {
        threshold = peak_db + d_thresh_min;
    }

    bool meas_err(threshold <= noise_floor_db);
    threshold = pow(10, threshold / 10.0);

    // now determine the first index below the threshold in either direction
    auto hi_bin(peak);
    auto lo_bin(peak);
    for (; hi_bin < mags2.end()-1; hi_bin++)
        if (*hi_bin <= threshold) break;
    for (; lo_bin > mags2.begin(); lo_bin--)
        if (*lo_bin <= threshold) break;
    int hi_idx = peak_idx + std::distance(peak, hi_bin);
    int lo_idx = peak_idx + std::distance(peak, lo_bin);

    for (int i = lo_idx; i <= hi_idx; i++)
        d_bug[i] = std::real(d_bug[i]) + 1j * 10*log10(threshold);
    d_bug[peak_idx] = std::real(d_bug[peak_idx]) + 1j * 10*log10(*peak);

    // bandwidth gets compensated for the linear-interpolated distance to the threshold
    // unless the burst extends to the edge of the PSD
    float bw_interp_offset = 0.0;
    if (lo_idx > 0) {
        bw_interp_offset += (threshold - *lo_bin) / (*(lo_bin + 1) - *lo_bin);
    }
    if (hi_idx < (mags2.size() - 1)) {
        bw_interp_offset += (threshold - *hi_bin) / (*(hi_bin - 1) - *hi_bin);
    }
    bandwidth = (hi_idx - lo_idx - bw_interp_offset) * bin_resolution;

    if (bandwidth <=0) { // this should never happen, but don't return a negative BW
        meas_err = true;
        bandwidth = 0.1;
    }

    // do linear interpolation to compensate for disparate magnitudes at threshold bins
    float sbi_factor = 0;  // sub-bin interpolation factor, %FFT bin to adjust high/low bin by
    if (*hi_bin > *lo_bin) {
        // if the higher bin is closer to the threshold, we will be shifting lower bin up in freq
        sbi_factor = (*hi_bin - *lo_bin) / (*(lo_bin + 1) - *lo_bin);
    } else {
        // if the lower bin is closer to the threshold, we will be shifting higher bin down in freq
        sbi_factor = -1 * (*lo_bin - *hi_bin) / (*(hi_bin - 1) - *hi_bin);
    }
    shift = ((hi_idx + lo_idx + sbi_factor) * 0.5 - (mags2.size() - 1) * 0.5) / mags2.size();

    return meas_err;
} /* end middle_out */

/*
 * FIXME: using the frequency axis is an odd way to do this... in general the logic here
 * could use some cleanup
 */
bool cf_estimate_impl::estimate_pwr(const std::vector<float> &mags2,
                                     const std::vector<float> &freq_axis,
                                     float center_frequency,
                                     float bandwidth,
                                     float &power)
{
    double start_freq = center_frequency - bandwidth / 2.0;
    double stop_freq = center_frequency + bandwidth / 2.0;
    double signal_power = 0.0;
    for (auto i = 0; i < mags2.size(); i++) {
        if ((freq_axis[i] > start_freq) && (freq_axis[i] < stop_freq)) {
            signal_power += mags2[i];
        }
    }
    // if nothing was summed above, set the power to the peak bin value so its not zero
    if(signal_power == 0) {
        power = *std::max_element(std::begin(mags2), std::end(mags2));
        return true;  // flag measurement as an error
    }

    // return the total signal power
    power = 10 * log10(signal_power);
    return false;
} /* end estimate_pwr */

} /* namespace fhss_utils */
} /* namespace gr */
