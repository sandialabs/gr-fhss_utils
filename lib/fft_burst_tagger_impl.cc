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

#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>
#include <gnuradio/io_signature.h>
#include <fhss_utils/constants.h>

#include "fft_burst_tagger_impl.h"

#include <volk/volk.h>

#include <inttypes.h>
#include <stdio.h>

#ifdef __AVX__
#include <immintrin.h>
#endif

namespace gr {
namespace fhss_utils {

fft_burst_tagger::sptr fft_burst_tagger::make(float center_freq,
                                              int fft_size,
                                              int sample_rate,
                                              int burst_pre_len,
                                              int burst_post_len,
                                              int burst_width,
                                              int max_bursts,
                                              int max_burst_len,
                                              float threshold,
                                              int history_size,
                                              int lookahead,
                                              bool debug)
{
    return gnuradio::get_initial_sptr(new fft_burst_tagger_impl(center_freq,
                                                                fft_size,
                                                                sample_rate,
                                                                burst_pre_len,
                                                                burst_post_len,
                                                                burst_width,
                                                                max_bursts,
                                                                max_burst_len,
                                                                threshold,
                                                                history_size,
                                                                lookahead,
                                                                debug));
}

/*
 * The private constructor
 */
fft_burst_tagger_impl::fft_burst_tagger_impl(float center_freq,
                                             int fft_size,
                                             int sample_rate,
                                             int burst_pre_len,
                                             int burst_post_len,
                                             int burst_width,
                                             int max_bursts,
                                             int max_burst_len,
                                             float threshold,
                                             int history_size,
                                             int lookahead,
                                             bool debug)
    : gr::block("fft_burst_tagger",
                gr::io_signature::make(1, 1, sizeof(gr_complex)),
                gr::io_signature::make(1, 1, sizeof(gr_complex))),
      d_center_freq(center_freq),
      d_sample_rate(sample_rate),
      d_fft_size(fft_size),
      d_burst_pre_len(burst_pre_len),
      d_lookahead(lookahead),
      d_burst_id(0),
      d_pre_burst_id(0),
      d_n_tagged_bursts(0),
      d_abs_fft_index(0),
      d_max_burst_len(max_burst_len),
      d_fft(NULL),
      d_history_size(history_size),
      d_history_primed(false),
      d_history_index(0),
      d_burst_post_len(burst_post_len),
      d_debug(debug),
      d_pub_debug(false),
      d_burst_debug_file(NULL),
      d_mask_owners(d_fft_size)

{
    const int nthreads = 1;
    set_output_multiple(d_fft_size);

    /*
     * The fine FFT is used to make a slightly more precise bandwidth and frequency
     * estimate at the start of the detected burst than the coarse FFT that is
     * computed on every sample of data. There will be d_lookahead samples of data
     * between the start of energy (pre burst) and point of declaring a burst to use.
     */
    d_fine_fft_size =
        d_fft_size * std::min(16, (int)pow(2, int(log(d_lookahead) / log(2))));
    d_fft = new fft::fft_complex(d_fft_size, true, nthreads);
    d_fine_fft = new fft::fft_complex(d_fine_fft_size, true, nthreads);

#ifdef __USE_MKL__
    MKL_LONG status =
        DftiCreateDescriptor(&m_fft, DFTI_SINGLE, DFTI_COMPLEX, 1, d_fft_size);
    status = DftiSetValue(m_fft, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiCommitDescriptor(m_fft);
    status =
        DftiCreateDescriptor(&m_fine_fft, DFTI_SINGLE, DFTI_COMPLEX, 1, d_fine_fft_size);
    status = DftiSetValue(m_fine_fft, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiCommitDescriptor(m_fine_fft);
#endif

    /*
     * We need to keep a number of samples in the buffer to ensure we are not emitting
     * data that might contain bursts that have not been tagged yet. Because bursts are
     * tagged after being active for  `d_lookahead+1`  FFTs, and that tag is placed
     * `d_burst_pre_len`  samples ahead of the detected start of the burst, we always need
     * have this amount of samples available at the start of our input buffer. This also
     * means that we start processing this many items into the work function, but emit
     * data starting back at the beginning of the work function.
     */
    d_work_history_nffts = (d_lookahead + d_burst_pre_len + 1);
    set_history(d_work_history_nffts * d_fft_size + 1);
    d_work_sample_offset = d_work_history_nffts * d_fft_size;

    // setup the normalized FFT windows and pre-compensate for the FFT size
    // TODO: these may not be the best window to use for this application...
    d_window_f = (float*)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
    std::vector<float> win =
        fft::window::build(fft::window::WIN_BLACKMAN, d_fft_size, 0);
    double pwr_acc = 0.0;
    for (auto x: win) pwr_acc += x*x/d_fft_size;
    for (auto ii=0; ii < d_fft_size; ii++) {
        d_window_f[ii] = win[ii] / (std::sqrt(pwr_acc) * d_fft_size);
    }
    d_fine_window_f =
        (float*)volk_malloc(sizeof(float) * d_fine_fft_size, volk_get_alignment());
    std::vector<float> fwin =
        fft::window::build(fft::window::WIN_BLACKMAN, d_fine_fft_size, 0);
    pwr_acc = 0.0;
    for (auto x: fwin) pwr_acc += x*x/d_fine_fft_size;
    for (auto ii=0; ii < d_fine_fft_size; ii++) {
        d_fine_window_f[ii] = fwin[ii] / (std::sqrt(pwr_acc) * d_fine_fft_size);
    }

    d_baseline_sum_f =
        (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
    d_magnitude_shifted_f =
        (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
    d_fine_magnitude_shifted_f =
        (float *)volk_malloc(sizeof(float) * d_fine_fft_size, volk_get_alignment());
    d_relative_magnitude_f =
        (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
    d_burst_mask_i =
        (uint32_t *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
    d_burst_mask_j =
        (uint32_t *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());

    d_bin_width_db = 10 * log10(d_sample_rate / d_fft_size);
    d_rel_mag_hist = 1;
    d_rel_hist_index = 0;
    d_relative_history_f =
        (float*)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());

    memset(d_baseline_sum_f, 0, sizeof(float) * d_fft_size);
    memset(d_magnitude_shifted_f, 0, sizeof(float) * d_fft_size);
    memset(d_fine_magnitude_shifted_f, 0, sizeof(float) * d_fft_size);
    memset(d_relative_magnitude_f, 0, sizeof(float) * d_fft_size);
    extra = 0;

    memset(d_burst_mask_i, ~0, sizeof(uint32_t) * d_fft_size);
    memset(d_burst_mask_j, ~0, sizeof(uint32_t) * d_fft_size);

    d_threshold = pow(10, threshold / 10);

    if (d_debug) {
        GR_LOG_INFO(d_logger,
                    boost::format("threshold=%f, d_threshold=%f (%f/%d)") % threshold %
                        d_threshold % (d_threshold * d_history_size) % d_history_size);
    }

    d_peaks.resize(d_fft_size);
    d_current_peaks = 0;
    d_bin_averages.resize(d_fft_size, moving_average(d_history_size));

    // this is the maximum number of simultaneous bursts
    if (max_bursts) {
        d_max_bursts = max_bursts;
    } else {
        // Consider the data to be invalid and reset the detector if more than 80% of all
        // channels are in use, this is useful when saturation occurrs
        d_max_bursts = (sample_rate / burst_width) * 0.8;
    }

    for (size_t i = 0; i < d_fft_size; i++)
        d_mask_owners[i].uid = i;

    // Area to ignore around an already found signal in FFT bins, rounded up to multiple
    // of two, which also helps compensate for the window width
    d_burst_width = 2 * std::ceil(1.0 * burst_width / (2.0 * sample_rate / fft_size));
    GR_LOG_INFO(d_logger, boost::format("bursts width is +/- %d FFT bins") % (d_burst_width / 2));

    d_filter_bandwidth = 0;

    if (d_debug) {
        GR_LOG_INFO(d_logger, boost::format("d_max_bursts=%d") % d_max_bursts);
    }

    if (d_debug) {
        d_burst_debug_file = fopen("/tmp/fft_burst_tagger-bursts.log", "w");
    }

    message_port_register_out(PMTCONSTSTR__debug());
}

/*
 * Our virtual destructor.
 */
fft_burst_tagger_impl::~fft_burst_tagger_impl()
{
    GR_LOG_INFO(d_logger, boost::format("Tagged %lu bursts") % d_n_tagged_bursts);
    delete d_fft;
    delete d_fine_fft;
    volk_free(d_window_f);
    volk_free(d_baseline_sum_f);
    volk_free(d_relative_magnitude_f);
    volk_free(d_relative_history_f);
    volk_free(d_magnitude_shifted_f);
    volk_free(d_fine_window_f);
    volk_free(d_fine_magnitude_shifted_f);
    volk_free(d_burst_mask_i);
    volk_free(d_burst_mask_j);
    if (d_burst_debug_file) {
        fclose(d_burst_debug_file);
    }
#ifdef __USE_MKL__
    MKL_LONG status = DftiFreeDescriptor(&m_fft);
    status = DftiFreeDescriptor(&m_fine_fft);
#endif
}

bool fft_burst_tagger_impl::start()
{
    GR_LOG_INFO(d_logger, boost::format("total time: %f") % d_total_timer.elapsed());

    // check if there are any blocks subscribed to the `debug` port - if not, disable output
    if(pmt::is_pair(pmt::dict_ref(d_message_subscribers, PMTCONSTSTR__debug(), pmt::PMT_NIL))) {
        GR_LOG_INFO(d_logger, "Debug Message Port is connected - Enabling Output");
        d_pub_debug = true;
    } else {
        GR_LOG_INFO(d_logger, "Debug Message Port is not connected - Disabling Output");
        d_pub_debug = false;
    }
    return true;
}

bool fft_burst_tagger_impl::stop()
{
#ifdef DO_TIMER
    GR_LOG_INFO(d_logger, boost::format("total time: %f") % d_total_timer.elapsed());
    GR_LOG_INFO(d_logger, boost::format("rel mag: %f") % d_rel_mag_timer.elapsed());
    GR_LOG_INFO(d_logger, boost::format("fft time: %f") % d_fft_timer.elapsed());
    GR_LOG_INFO(d_logger,
                boost::format("update pb time: %f") % d_update_pb_timer.elapsed());
    GR_LOG_INFO(d_logger,
                boost::format("update ab time: %f") % d_update_ab_timer.elapsed());
    GR_LOG_INFO(d_logger,
                boost::format("remove tb time: %f") % d_remove_tb_timer.elapsed());
    GR_LOG_INFO(d_logger, boost::format("extract time: %f") % d_extract_timer.elapsed());
    GR_LOG_INFO(d_logger, boost::format("delete time: %f") % d_delete_timer.elapsed());
    GR_LOG_INFO(d_logger, boost::format("new pb time: %f") % d_new_pb_timer.elapsed());
    GR_LOG_INFO(d_logger, boost::format("new burst time: %f") % d_new_b_timer.elapsed());
    GR_LOG_INFO(d_logger, boost::format("update cb: %f") % d_update_cb_timer.elapsed());
    GR_LOG_INFO(d_logger, boost::format("other timer: %f") % d_other.elapsed());
#endif

    GR_LOG_INFO(d_logger, boost::format("saw %lu ffts") % d_abs_fft_index);
    GR_LOG_INFO(d_logger, boost::format("extra = %zu") % extra);

    return true;
}

/*
 *  Compute the magnitude of current data normalized by nearby magnitudes
 */
bool fft_burst_tagger_impl::compute_relative_magnitude(void)
{
    // Only compute if we have sufficient data in our circular buffer
    d_rel_mag_timer.start();
    if (!d_history_primed) {
        return false;
    }

    volk_32f_x2_divide_32f(
        d_relative_magnitude_f, d_magnitude_shifted_f, d_baseline_sum_f, d_fft_size);
    d_rel_mag_timer.end();
    return true;
}

void fft_burst_tagger_impl::update_circular_buffer(void)
{
    // We only update the average if there is no burst going on at the moment
    // TODO: also block if there are any pre_bursts. Is that what we want?
    if (!d_history_primed) {
        for (size_t i = 0; i < d_fft_size; i++) {
            d_baseline_sum_f[i] = d_bin_averages[i].add(d_magnitude_shifted_f[i]);
        }
        d_history_index++;

        if (d_history_index == d_history_size) {
            d_history_primed = true;
        }
    } else if ((d_abs_fft_index & 0x3) == 0) {
        for (size_t i = 0; i < d_fft_size; i++) {
            if (d_burst_mask_i[i] != 0 && d_burst_mask_j[i] != 0) {
                d_baseline_sum_f[i] = d_bin_averages[i].add(d_magnitude_shifted_f[i]);
            }
        }
        memcpy(d_burst_mask_j, d_burst_mask_i, d_fft_size * sizeof(uint32_t));
    }
    // Copy the relative magnitude history into circular buffer
    memcpy(d_relative_history_f + d_rel_hist_index * d_fft_size,
           d_relative_magnitude_f,
           sizeof(float) * d_fft_size);
    if (d_rel_mag_hist > 0)
        d_rel_hist_index = (d_rel_hist_index + 1) % d_rel_mag_hist;
}

void fft_burst_tagger_impl::update_active_bursts(void)
{
    auto b = std::begin(d_bursts);
    while (b != std::end(d_bursts)) {
        if (d_relative_magnitude_f[b->center_bin - 1] > d_threshold ||
            d_relative_magnitude_f[b->center_bin] > d_threshold ||
            d_relative_magnitude_f[b->center_bin + 1] > d_threshold) {
            b->last_active = d_abs_fft_index;
        }
        ++b;
    }
}

void fft_burst_tagger_impl::reset()
{
    std::lock_guard<std::mutex> lock(d_work_mutex);
    _reset();
}

void fft_burst_tagger_impl::_reset()
{
    GR_LOG_WARN(d_logger, "=========== RESETTING BURST DETECTOR ===========");
    // close and tag all current bursts
    // d_peaks.clear();
    d_current_peaks = 0;
    d_pre_bursts.clear();
    d_new_bursts.clear();
    for (burst b : d_bursts) {
        b.stop = d_abs_fft_index;
        if (b.valid)
            d_gone_bursts.push_back(b);
    }
    d_bursts.clear();
    for (size_t i = 0; i < d_fft_size; i++)
        d_mask_owners[i].clear();

    tag_gone_bursts(d_abs_fft_index * d_fft_size);

    // reset baseline
    d_history_index = 0;
    d_history_primed = false;
    memset(d_baseline_sum_f, 0, sizeof(float) * d_fft_size);
    memset(d_burst_mask_i, ~0, sizeof(uint32_t) * d_fft_size);
    memset(d_burst_mask_j, ~0, sizeof(uint32_t) * d_fft_size);
}

void fft_burst_tagger_impl::delete_gone_bursts(void)
{
    auto b = std::begin(d_bursts);

    while (b != std::end(d_bursts)) {
        if ((b->last_active + d_burst_post_len) < d_abs_fft_index ||
            (d_max_burst_len && d_abs_fft_index - b->start > d_max_burst_len)) {
            b->stop = d_abs_fft_index;
            if (b->valid)
                d_gone_bursts.push_back(*b);
            remove_ownership(*b);
            b = d_bursts.erase(b);
        } else {
            ++b;
        }
    }
}

bool fft_burst_tagger_impl::check_prev_magnitude(size_t bin)
{
    // This is a cicular buffer, order doesn't really matter though.
    size_t cIndex = 0;
    size_t start_bin = std::max(bin - d_burst_width / 2, (size_t)0);
    size_t stop_bin = std::min(bin + d_burst_width / 2, (size_t)(d_fft_size - 1));
    for (size_t i = 0; i < d_rel_mag_hist; i++) {
        bool found = false;
        for (size_t b = start_bin; b <= stop_bin; b++) {
            if (d_relative_history_f[cIndex + b] > d_threshold) {
                found = true;
            }
        }
        cIndex += d_fft_size;
        if (found == false)
            return false;
    }
    return true;
}

void fft_burst_tagger_impl::create_new_potential_bursts(void)
{
    // A potential burst is one that we are not yet sure about.  The burst needs to exist
    // for a sufficent amount of time in order for us to keep it.
    bool allow[d_fft_size];
    memset(allow, 1, sizeof(bool) * d_fft_size);
    for (size_t i = 0; i < d_current_peaks; i++) {
        peak& p = d_peaks[i];
        if (d_burst_mask_i[p.bin] && allow[p.bin]) {
            if (check_prev_magnitude(p.bin)) {
                pre_burst b;
                b.center_bin = p.bin;
                b.start_bin = std::max(b.center_bin - d_burst_width / 2, 0);
                b.stop_bin = std::min(b.center_bin + d_burst_width / 2, d_fft_size - 1);

                b.peak_count = 1 + d_rel_mag_hist;
                b.max_relative_magnitude = p.relative_magnitude;
                size_t start_bin = b.start_bin;
                size_t stop_bin = b.stop_bin;

                // find the largest peak within d_burst_width/2 bins on either side of the
                // last detection
                for (size_t j = 0; j < d_rel_mag_hist; j++) {
                    for (int k = start_bin; k < stop_bin; k++) {
                        if (d_relative_history_f[k] > b.max_relative_magnitude) {
                            b.max_relative_magnitude = d_relative_history_f[k];
                            b.center_bin = k;
                        }
                    }
                }

                // The burst might have started a few FFTs earlier so offset by d_rel_mag_hist
                long fftidx = d_abs_fft_index - d_burst_pre_len - d_rel_mag_hist;
                // when using noise floor preload, the burst start should not be negative
                b.start = std::max(fftidx, 0L);
                b.id = d_pre_burst_id++;
                b.thresh_count = 1 + d_rel_mag_hist;

                d_pre_bursts.push_back(b);
                add_ownership(b);

                if (d_burst_debug_file) {
                    fprintf(
                        d_burst_debug_file, "%" PRIu64 ",%d,x\n", b.start, b.center_bin);
                }
            } else {
                size_t start_bin = std::max((int)p.bin - d_burst_width / 2, 0);
                size_t stop_bin =
                    std::min((int)p.bin + d_burst_width / 2, (d_fft_size - 1));
                for (size_t j = start_bin; j <= stop_bin; j++) {
                    allow[j] = 0;
                }
            }
        }
    }
    // TODO: move this
    if (d_max_bursts > 0 && d_bursts.size() > d_max_bursts) {
        _reset();
    }
}

void fft_burst_tagger_impl::add_ownership(const pre_burst& b)
{
    for (size_t i = b.start_bin; i <= b.stop_bin; i++) {
        d_burst_mask_i[i] = 0;
        d_burst_mask_j[i] = 0;
        if (d_mask_owners[i].push_back(b.id)) {
            GR_LOG_WARN(d_logger,
                        boost::format("Owners::push_back - Trying to add too many points "
                                      "for  bin id=%zu, size=%zu, burst=%zu") %
                            d_mask_owners[i].uid % d_mask_owners[i]._size % b.id);
            GR_LOG_WARN(d_logger,
                        boost::format("Owner bins are %zu, %zu, %zu %zu") %
                            d_mask_owners[i].ids[0] % d_mask_owners[i].ids[1] %
                            d_mask_owners[i].ids[2] % d_mask_owners[i].ids[3]);
        }
    }
}

void fft_burst_tagger_impl::update_ownership(const pre_burst& pb, const burst& b)
{
    if (pb.start_bin != b.start_bin || pb.stop_bin != b.stop_bin) {
        // There is a chance that we will have two pre-bursts that
        for (size_t i = pb.start_bin; i <= pb.stop_bin; i++) {
            if (d_mask_owners[i].size() == 1) {
                d_burst_mask_i[i] = ~0;
            }
            if (d_mask_owners[i].erase(pb.id)) {
                GR_LOG_ERROR(
                    d_logger,
                    boost::format("Owners::Remove - Couldn't find id to erase. This "
                                  "should never happen - bin id = %zu, burst=%zu") %
                        d_mask_owners[i].uid % pb.id);
            }
        }
        for (size_t i = b.start_bin; i <= b.stop_bin; i++) {
            d_burst_mask_i[i] = 0;
            d_burst_mask_j[i] = 0;
            if (d_mask_owners[i].push_back(b.id)) {
                GR_LOG_WARN(d_logger,
                            boost::format("Owners::push_back - Trying to add too many "
                                          "points for  bin id=%zu, size=%zu, burst=%zu") %
                                d_mask_owners[i].uid % d_mask_owners[i]._size % b.id);
                GR_LOG_WARN(d_logger,
                            boost::format("Owner bins are %zu, %zu, %zu %zu") %
                                d_mask_owners[i].ids[0] % d_mask_owners[i].ids[1] %
                                d_mask_owners[i].ids[2] % d_mask_owners[i].ids[3]);
            }
        }
    } else {
        for (size_t i = b.start_bin; i <= b.stop_bin; i++) {
            if (d_mask_owners[i].update(pb.id, b.id)) {
                GR_LOG_ERROR(d_logger,
                             "Owners::Update - Couldn't find id to update. This should "
                             "never happen");
            }
        }
    }
}

void fft_burst_tagger_impl::remove_ownership(const pre_burst& b)
{
    for (size_t i = b.start_bin; i <= b.stop_bin; i++) {
        if (d_mask_owners[i].size() == 1) {
            d_burst_mask_i[i] = ~0;
        }
        if (d_mask_owners[i].erase(b.id)) {
            GR_LOG_ERROR(d_logger,
                         boost::format("Owners::Remove - Couldn't find id to erase. This "
                                       "should never happen - bin id = %zu, burst=%zu") %
                             d_mask_owners[i].uid % b.id);
        }
    }
}

void fft_burst_tagger_impl::remove_ownership(const burst& b)
{
    for (size_t i = b.start_bin; i <= b.stop_bin; i++) {
        if (d_mask_owners[i].size() == 1) {
            d_burst_mask_i[i] = ~0;
        }
        if (d_mask_owners[i].erase(b.id)) {
            GR_LOG_ERROR(d_logger,
                         boost::format("Owners::Remove - Couldn't find id to erase. This "
                                       "should never happen - bin id = %zu, burst=%zu") %
                             d_mask_owners[i].uid % b.id);
        }
    }
}

/*
 * Convert potential bursts into burst once they have existed for a sufficient amount of
 * time.
 */
void fft_burst_tagger_impl::create_new_bursts(const gr_complex* input, int fft)
{
    auto b = d_pre_bursts.begin();
    while (b != d_pre_bursts.end()) {

        // if d_lookahead peaks have been measured, convert the pre_burst into a burst
        if (b->peak_count >= d_lookahead + 1) {
            burst new_b;
            new_b.start = b->start;
            new_b.last_active = d_abs_fft_index;

            // use the max peak's bin and magnitude for the burst
            // take the last N samples up to the end of the current FFT.
            int offset = (fft + 1) * d_fft_size - d_fine_fft_size;
            volk_32fc_32f_multiply_32fc(d_fine_fft->get_inbuf(),
                                        &input[offset],
                                        d_fine_window_f,
                                        d_fine_fft_size);

#ifdef __USE_MKL2__
            DftiComputeForward(
                m_fine_fft, d_fine_fft->get_inbuf(), d_fine_fft->get_outbuf());
#else
            d_fine_fft->execute();
#endif

            // Get the mag squared of the fft and fft shift
            size_t fftd2 = d_fine_fft_size / 2;
            volk_32fc_magnitude_squared_32f(
                d_fine_magnitude_shifted_f + fftd2, d_fine_fft->get_outbuf(), fftd2);
            volk_32fc_magnitude_squared_32f(
                d_fine_magnitude_shifted_f, d_fine_fft->get_outbuf() + fftd2, fftd2);
            // Search around the peak to try to estimate the bandwidth
            size_t factor = d_fine_fft_size / d_fft_size;
            // Ensure that a center bin - burst_width < 0 becomes a 1, not a very large
            // positive number.
            size_t search_start =
                (size_t)std::max((int)factor * (b->center_bin - d_burst_width / 2), 1);
            // Ensure that search_stop is not 0 by performing ceil() division on
            // (d_burst_width / 2) note: int(true) == 1 in C++
            size_t search_stop = std::min(factor * (b->center_bin + (d_burst_width / 2) +
                                                    int(d_burst_width % 2 != 0)),
                                          (size_t)(d_fine_fft_size - 2));

            // Normalize the magnitude by the noise floor
            for (size_t i = search_start; i < search_stop; i++) {
                d_fine_magnitude_shifted_f[i] /=
                    d_baseline_sum_f[i / factor] / d_lookahead;
            }

            float* max_entry = std::max_element(d_fine_magnitude_shifted_f + search_start,
                                                d_fine_magnitude_shifted_f + search_stop);
            // Smooth to get max value.  This is safe because of the bounds on search
            // start and stop.
            float max_value = .5 * max_entry[0] + .25 * (max_entry[-1] + max_entry[1]);
            // How high is this above the noise floor?
            float snr_est = max_value;
            float* n = d_baseline_sum_f + b->center_bin;

            size_t minIndex = search_start;
            size_t maxIndex = search_stop;
            // 6 db from peak
            float thresh = max_value * .0625;
            // float thresh = (max_value - noise_floor) * 0.0625;
            while (minIndex <= search_stop &&
                   d_fine_magnitude_shifted_f[minIndex] < thresh) {
                minIndex++;
            }
            while (maxIndex >= search_start &&
                   d_fine_magnitude_shifted_f[maxIndex] < thresh) {
                maxIndex--;
            }
            new_b.bandwidth = (maxIndex - minIndex + 1) * d_sample_rate / d_fine_fft_size;
            new_b.center_freq =
                ((ssize_t)(maxIndex + minIndex) - (ssize_t)d_fine_fft_size) / 2 *
                d_sample_rate / d_fine_fft_size;

            new_b.center_bin = (maxIndex + minIndex) / 2 / factor;
            new_b.start_bin = std::max(new_b.center_bin - d_burst_width / 2, 0);
            new_b.stop_bin =
                std::min(new_b.center_bin + d_burst_width / 2, d_fft_size - 1);
            float max_relative_magnitude = b->max_relative_magnitude;

            // now that we know the center bin, make sure this burst doesn't already exist
            // TODO: improve this, maybe use the mask
            bool already_exists = false;
            for (burst existing_b : d_bursts) {
                if (abs(new_b.center_bin - existing_b.center_bin) <=
                    (d_burst_width + 1) / 2) {
                    already_exists = true;
                    break;
                }
            }
            if (already_exists) {
                remove_ownership(*b);
                b = d_pre_bursts.erase(b);
                continue;
            }

            new_b.id = d_burst_id++;

            // normalize the magnitude
            new_b.magnitude = 10 * log10(max_relative_magnitude * d_history_size);
            double noise = 0;
            for (int i = new_b.start_bin; i <= new_b.stop_bin; i++) {
                noise += d_baseline_sum_f[i];
            }
            // normalize by total number of bins that form the sum
            new_b.noise = noise / (new_b.stop_bin - new_b.start_bin + 1);

            // Allow for us to filter out the data if needed
            if (d_filter_bandwidth > 0 && new_b.bandwidth > d_filter_bandwidth) {
                // keep tracking the burst but don't write out the tags so we don't
                // tune/filter/decimate
                new_b.valid = false;
            } else {
                d_new_bursts.push_back(new_b);
                new_b.valid = true;
            }
            d_bursts.push_back(new_b);

            update_ownership(*b, new_b);
            b = d_pre_bursts.erase(b);
            continue;
        }
        ++b;
    }
}

void fft_burst_tagger_impl::update_potential_bursts(void)
{
    float* c_magnitude = d_relative_magnitude_f;
    float threshold_low = d_threshold * .25;
    auto b = d_pre_bursts.begin();
    while (b != d_pre_bursts.end()) {

        // initialize a peak to store the largest new peak in
        peak p;
        size_t start_bin = b->start_bin;
        size_t stop_bin = b->stop_bin;
        p.bin = stop_bin;
        p.relative_magnitude = c_magnitude[stop_bin];

        // find the largest peak within d_burst_width/2 bins on either side of the first
        // detection
        for (int i = start_bin; i < stop_bin; i++) {
            if (c_magnitude[i] > p.relative_magnitude) {
                p.relative_magnitude = c_magnitude[i];
                p.bin = i;
            }
        }

        // if the largest peak is below the threshold, drop the pre_burst
        if (p.relative_magnitude < threshold_low) {
            // erase

            // Think about not creating a potential burst until we see it a few times.  In
            // my dataset, 60% of pb vanish after 1 frame, and 90% after 4 frames. I could
            // also change this to be a hash table, with potentially cheaper add/delete.
            if (b->peak_count == 4)
                extra++;

            remove_ownership(*b);
            b = d_pre_bursts.erase(b);
            continue;
        } else {
            if (p.relative_magnitude > d_threshold) {
                b->thresh_count++;
            } else if (b->thresh_count < b->peak_count * .75) {
                remove_ownership(*b);
                b = d_pre_bursts.erase(b);
                continue;
            }
        }

        // otherwise, add the peak to the list of measurements
        b->peak_count++;
        if (p.relative_magnitude > b->max_relative_magnitude) {
            b->max_relative_magnitude = p.relative_magnitude;
            b->center_bin = p.bin;
        }
        ++b;
    }
}

/*
 * Remove peaks from consideration that we are already tracking.  This function is
 * called after we have updated bursts, so we only want to see new ones present. Volk
 * function requires ints as inputs, but we are really just zeroing out and bin that
 * we are already tracking
 */
void fft_burst_tagger_impl::remove_currently_tracked_bursts(void)
{
    volk_32i_x2_and_32i((int32_t*)d_relative_magnitude_f,
                        (int32_t*)d_relative_magnitude_f,
                        (int32_t*)d_burst_mask_i,
                        d_fft_size);
}

inline void swap_peaks(std::vector<peak>& p, size_t i, size_t j)
{
    peak temp = p[i];
    p[i] = p[j];
    p[j] = temp;
}

void fft_burst_tagger_impl::extract_peaks(void)
{
    const uint32_t start_bin = d_burst_width / 2;
    const uint32_t end_bin = d_fft_size - d_burst_width / 2;
    size_t length = end_bin - start_bin;
    size_t length8 = 0;
    size_t index = 0;

#ifdef __AVX__
    length8 = length >> 3;
    float* c_magnitude = d_relative_magnitude_f + start_bin;
    // Copy d_threshold to every position in 256 bit array
    __m256 threshold = _mm256_set1_ps(d_threshold);

    for (size_t j = 0; j < length8; j++, c_magnitude += 8) {
        // Unalgined load.  About ~1% slower than aligned load, not worth changing
        __m256 magnitude = _mm256_loadu_ps(c_magnitude);
        // Sub data from threshold
        __m256 diff = _mm256_sub_ps(threshold, magnitude);
        // negative numbers have a sign bit of 1.  Pull the sign bit off all 8 values and
        // store in an int.
        int which = _mm256_movemask_ps(diff);

        if (which) {
            size_t loc = start_bin + (j << 3);
            while (which) {
                if (which & 0x1) {
                    d_peaks[index].bin = loc;
                    d_peaks[index++].relative_magnitude = d_relative_magnitude_f[loc];
                }
                loc++;
                which >>= 1;
            }
        }
    }
#endif
    // Do the trailer and get the non multiple of 8 tail.
    for (size_t bin = (length8 << 3) + start_bin; bin < end_bin; bin++) {
        float rel_mag = d_relative_magnitude_f[bin];
        if (rel_mag > d_threshold) {
            d_peaks[index].bin = bin;
            d_peaks[index++].relative_magnitude = rel_mag;
        }
    }
    // We only need to sort if the last time we had peaks.  Otherwise there can't be any
    // that will pass this time.
    bool doSort = d_current_peaks > 0 && d_rel_mag_hist > 0;
    d_current_peaks = index;

    // sort by descending relative magnitude
    // Special case a few things here. To save time
    // Case 0,1: Just return
    // Case 2: swap if needed
    // Case 3 is a little more complicated, but not bad
    if (!doSort)
        return;
    if (index == 2) {
        if (d_peaks[0].relative_magnitude < d_peaks[1].relative_magnitude) {
            swap_peaks(d_peaks, 0, 1);
        }
    } else if (index == 3) {
        if (d_peaks[1].relative_magnitude > d_peaks[0].relative_magnitude) {
            swap_peaks(d_peaks, 0, 1);
        }
        if (d_peaks[2].relative_magnitude > d_peaks[0].relative_magnitude) {
            swap_peaks(d_peaks, 0, 2);
        }
        if (d_peaks[2].relative_magnitude > d_peaks[1].relative_magnitude) {
            swap_peaks(d_peaks, 2, 1);
        }
    } else if (index > 3) {
        std::sort(
            d_peaks.begin(), d_peaks.begin() + index, [](const peak& a, const peak& b) {
                return a.relative_magnitude > b.relative_magnitude;
            });
    }
}

void fft_burst_tagger_impl::save_peaks_to_debug_file(char* filename)
{
    FILE* file = fopen(filename, "a");
    for (size_t index = 0; index < d_current_peaks; index++) {
        const peak& p = d_peaks[index];
        fprintf(file, "%" PRIu64 ",%d,x\n", d_abs_fft_index, p.bin);
    }
    fclose(file);
}

void fft_burst_tagger_impl::tag_new_bursts(void)
{
    for (burst b : d_new_bursts) {
        pmt::pmt_t key = PMTCONSTSTR__new_burst();
        float relative_frequency = (b.center_freq) / d_sample_rate;

        /* NOISE FLOOR ESTIMATION
         *
         * The burst `noise` field contains a sum of d_history_size FFT values averaged
         * across the number of bins the burst exists in. We are now reporting the noise
         * as a density in dBFS/Hz so this needs to be normalized by the bandwidth per
         * bin. This is done using precomputed dB values for the bin width in dB.
         */
        double noise_density = 10 * log10(b.noise) - d_bin_width_db;

        pmt::pmt_t value = pmt::make_dict();
        value = pmt::dict_add(value, PMTCONSTSTR__burst_id(), pmt::from_uint64(b.id));
        value = pmt::dict_add(value,
                              PMTCONSTSTR__relative_frequency(),
                              pmt::from_float(relative_frequency));
        value = pmt::dict_add(
            value, PMTCONSTSTR__center_frequency(), pmt::from_float(d_center_freq));
        value =
            pmt::dict_add(value, PMTCONSTSTR__magnitude(), pmt::from_float(b.magnitude));
        value = pmt::dict_add(
            value, PMTCONSTSTR__sample_rate(), pmt::from_float(d_sample_rate));
        value = pmt::dict_add(
            value, PMTCONSTSTR__noise_density(), pmt::from_float(noise_density));
        value =
            pmt::dict_add(value, PMTCONSTSTR__bandwidth(), pmt::from_float(b.bandwidth));

        add_item_tag(0, b.start * d_fft_size, key, value);
    }
    d_new_bursts.clear();
}

void fft_burst_tagger_impl::tag_gone_bursts(int noutput_items)
{
    auto b = std::begin(d_gone_bursts);

    while (b != std::end(d_gone_bursts)) {
        uint64_t output_index = b->stop * d_fft_size;

        if (nitems_written(0) <= output_index &&
            output_index < nitems_written(0) + noutput_items) {
            pmt::pmt_t key = PMTCONSTSTR__gone_burst();
            pmt::pmt_t value = pmt::make_dict();
            value =
                pmt::dict_add(value, PMTCONSTSTR__burst_id(), pmt::from_uint64(b->id));

            add_item_tag(0, output_index, key, value);
            d_n_tagged_bursts++;

            b = d_gone_bursts.erase(b);
        } else {
            ++b;
        }
    }
}

void fft_burst_tagger_impl::publish_debug()
{
    if(!d_pub_debug) {
        return;
    }

    std::vector<gr_complex> dbg(d_fft_size);
    // generate output vector, real value is the current FFT, imaginary is the threshold
    // which allows this to be plotted on the time sink
    for (auto jj = 0; jj < d_fft_size; jj++) {
        dbg[jj] = 10 * log10(d_magnitude_shifted_f[jj]) +
                  1j * 10 * log10(1e-30 + 1.0 * d_baseline_sum_f[0 + jj]);
    }

    message_port_pub(PMTCONSTSTR__debug(),
                     pmt::cons(pmt::PMT_NIL, pmt::init_c32vector(d_fft_size, &dbg[0])));
}

uint64_t fft_burst_tagger_impl::get_n_tagged_bursts() { return d_n_tagged_bursts; }

void fft_burst_tagger_impl::preload_noise_floor(double noise_density, bool preload)
{
    // the preload boolean allows this function to be called but also bypassed (GRC thing...)
    if (preload) {
        GR_LOG_INFO(d_logger, boost::format("initializing noise denisty to %.2f dB") % noise_density);
        double noise_density_linear = pow(10, noise_density / 10.0);
        for (auto i = 0; i < d_fft_size; i++) {
            for (auto j = 0; j < d_history_size; j++) {
                d_bin_averages[i].add(noise_density_linear);
            }
        }
        d_history_primed = true;
    } else {
        GR_LOG_INFO(d_logger, "preload called but bypassed, running normally");
    }
}

int fft_burst_tagger_impl::general_work(int noutput_items,
                                        gr_vector_int& ninput_items,
                                        gr_vector_const_void_star& input_items,
                                        gr_vector_void_star& output_items)
{
    std::lock_guard<std::mutex> lock(d_work_mutex);
    d_total_timer.start();

    const gr_complex* in = (const gr_complex*)input_items[0];
    const gr_complex* in_new = in + d_work_sample_offset; // first new input sample
    gr_complex* out = (gr_complex*)output_items[0];

    int nffts = std::max(
        (std::min((int)(ninput_items[0] - d_work_sample_offset), noutput_items)) /
            d_fft_size,
        0);
    int nitems = nffts * d_fft_size;

    int produced = nitems;
    // no data is produced until we work through the history
    if (d_abs_fft_index < d_work_history_nffts) {
        produced = 0;
        nffts = std::min(nffts, (int)(d_work_history_nffts - d_abs_fft_index));
        nitems = nffts * d_fft_size;
    }

    // we always consume however many FFTs were processed, after starting we consume this
    // many also
    int consumed = nitems;

    std::vector<tag_t> freq_tags;
    get_tags_in_window(freq_tags,
                       0,
                       d_work_sample_offset,
                       d_work_sample_offset + nitems,
                       PMTCONSTSTR__rx_freq());
    if (freq_tags.size() > 0) {
        d_center_freq = pmt::to_double(freq_tags[freq_tags.size() - 1].value);
    }

    // loop through input data processing it one FFT at a time
    for (int ii = 0; ii < nffts; ii++) {

        // Apply a blackman window and take fft
        d_fft_timer.start();
        volk_32fc_32f_multiply_32fc(
            d_fft->get_inbuf(), &in_new[ii * d_fft_size], d_window_f, d_fft_size);
#ifdef __USE_MKL__
        DftiComputeForward(m_fft, d_fft->get_inbuf(), d_fft->get_outbuf());
#else
        d_fft->execute();
#endif

        // Get the mag squared of the fft and fft shift
        size_t fftd2 = d_fft_size / 2;
        volk_32fc_magnitude_squared_32f(
            d_magnitude_shifted_f + fftd2, d_fft->get_outbuf(), fftd2);
        volk_32fc_magnitude_squared_32f(
            d_magnitude_shifted_f, d_fft->get_outbuf() + fftd2, fftd2);
        d_fft_timer.end();

        // this will be skipped if the noise floor estimate is not valid
        if (compute_relative_magnitude()) {

            /* checks pre burst list to update bursts that are still present */
            d_update_pb_timer.start();
            update_potential_bursts();
            d_update_pb_timer.end();

            /* checks burst list to update bursts that are still present */
            d_update_ab_timer.start();
            update_active_bursts();
            d_update_ab_timer.end();

            if (d_debug) {
                extract_peaks();
                save_peaks_to_debug_file((char*)"/tmp/fft_burst_tagger-peaks.log");
            }

            /* prunes the list of current bursts */
            d_remove_tb_timer.start();
            remove_currently_tracked_bursts();
            d_remove_tb_timer.end();

            d_extract_timer.start();
            extract_peaks();
            d_extract_timer.end();

            if (d_debug) {
                save_peaks_to_debug_file(
                    (char*)"/tmp/fft_burst_tagger-peaks-filtered.log");
            }

            /* if any bursts are gone, this will remove them */
            d_delete_timer.start();
            delete_gone_bursts();
            d_delete_timer.end();

            /* checks FFT data to see if any new pre bursts have been identified */
            d_new_pb_timer.start();
            create_new_potential_bursts();
            d_new_pb_timer.end();

            /* checks pre burst list to see if any should be converted to bursts */
            d_new_b_timer.start();
            create_new_bursts(in_new, ii);
            d_new_b_timer.end();
        }

        /* updates the dynamic noise floor estimate */
        d_update_cb_timer.start();
        update_circular_buffer();
        d_update_cb_timer.end();

        /* uncomment this to see every single FFT and corresponding noise average; this
         * will make online processing effectively unusable but can be very helpful to
         * debug internal behavior when operating from file.
         */
        // if (d_debug) { publish_debug(); usleep(25000); }

        d_abs_fft_index++;
    }

    // TODO: these may write tags into the future, is that a problem?
    d_other.start();
    tag_new_bursts();
    tag_gone_bursts(produced);

    publish_debug();

    d_other.end();
    d_total_timer.end();

    if (nffts == 0) {
        usleep(100);
    }

    consume_each(consumed);

    if (produced) {
        memcpy(out, in, sizeof(gr_complex) * produced);
    }
    return produced;
}

} /* namespace fhss_utils */
} /* namespace gr */
