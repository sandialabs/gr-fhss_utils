
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include <gnuradio/math.h>
#include <volk/volk.h>
#include "fine_burst_measure_impl.h"
#include <gnuradio/fft/window.h>

namespace gr {
  namespace fhss_utils {

    // The bounded histogram assumes that min and max values come from the data and thus will never be exceeded.  If you cannot
    // guarantee that the values will fall in this range, then you must use the unbounded version to prevent out of bounds memory
    // access
    void slow_histogram_bounded(float* data, size_t length, uint32_t* hist, float min_value, float max_value, size_t num_bins) {
      memset(hist, 0, num_bins*sizeof(hist[0]));
      float sub = min_value;
      float scale = (max_value - min_value) / (num_bins - 1);
      for (size_t i = 0; i < length; i++) {
        // Round to the nearest bin
        size_t bin = scale*(data[i] - sub) + 0.5;
        hist[bin]++;
      }
    }
    
    void slow_histogram_unbounded(float* data, size_t length, uint32_t* hist, float min_value, float max_value, size_t num_bins) {
      memset(hist, 0, num_bins*sizeof(hist[0]));
      float sub = min_value;
      float scale = 1.0 / ((max_value - min_value) / (num_bins - 1));
      for (size_t i = 0; i < length; i++) {
        // Round to the nearest bin
        float cdata = std::max(min_value, std::min(max_value, data[i]));
        size_t bin = scale*(cdata - sub) + 0.5;
        hist[bin]++;
      }
    }

    /*void fast_histogram(float* data, size_t length, float sub, float scale) {
      // Make this a volk kernel
      memset(d_histogram, 0, 4*d_histogram_bins*sizeof(d_histogram[0]));
      size_t length4 = (length >> 2) << 2;
      uint32_t* hist0 = d_histogram;
      uint32_t* hist1 = d_histogram + d_histogram_bins;
      uint32_t* hist2 = d_histogram + 2*d_histogram_bins;
      uint32_t* hist3 = d_histogram + 3*d_histogram_bins;

      // Allow for alignment - Loop over the first N elements to ensure we can do aligned loads.

      for (size_t i = 0; i < length4; i+=4) {
        // Round to the nearest bin
        // Use intrinsics to:
        // load
        // sub
        // mult
        // add
        // round (mm_cvtps_epi32)
        // store in array
        // inc histograms
        size_t bin = scale*(data[i] - sub) + 0.5;
        hist0[bin]++;
        bin = scale*(data[i+1] - sub) + 0.5;
        hist1[bin]++;
        bin = scale*(data[i+2] - sub) + 0.5;
        hist2[bin]++;
        bin = scale*(data[i+3] - sub) + 0.5;
        hist3[bin]++;
      }

      // Change to use intrinsic add functions 0+1, 2+3, and then those together
      
      for (size_t i = 0; i < d_histogram_bins; i++) {
        d_histogram[i] += hist1[i] + hist2[i] + hist3[i];
      }

      // Add in any tail elements in case length is not a mult of 4
      for ( ; i < length; i++) {
        size_t bin = scale*(data[i] - sub) + 0.5;
        d_histogram[bin]++;
    }*/

    fine_burst_measure::sptr
      fine_burst_measure::make(float freq_resolution, uint32_t analyze_samples, float threshold)
    {
      return gnuradio::get_initial_sptr
        (new fine_burst_measure_impl(freq_resolution, analyze_samples, threshold));
    }

    /*
     * The private constructor
     */
    fine_burst_measure_impl::fine_burst_measure_impl(float freq_resolution, uint32_t analyze_samples, float threshold)
      : gr::block("fine_burst_measure",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0)),
      d_freq_resolution(freq_resolution),
      d_analyze_samples(analyze_samples),
      d_threshold(threshold),
      d_phase_diff_f(nullptr),
      d_tmp_c(nullptr),
      d_min_freq(-M_PI),
      d_max_freq(M_PI)
    {
      message_port_register_in(pmt::mp("pdu_in"));
      message_port_register_out(pmt::mp("pdu_out"));
      set_msg_handler(pmt::mp("pdu_in"), boost::bind(&fine_burst_measure_impl::pdu_handler, this, _1));
      d_fft_size = d_analyze_samples;
      
      d_tmp_c = (gr_complex *)volk_malloc(d_analyze_samples * sizeof(gr_complex), volk_get_alignment());
      d_phase_diff_f = (float *)volk_malloc(d_analyze_samples * sizeof(float), volk_get_alignment());
      d_magnitude_f = (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
      d_magnitude_shifted_f = (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
      d_magnitude_smooth_f = (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());

      d_fft = new fft::fft_complex(d_fft_size, true, 1);
        
      d_window_f = (float *)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
      std::vector<float> window = fft::window::build(fft::window::WIN_BLACKMAN, d_fft_size, 0);
      memcpy(d_window_f, &window[0], sizeof(float) * d_fft_size);
      
      fid = fopen("/home/rfec/burst_data.dat", "wb");
      fid2 = fopen("/home/rfec/burst_hist.dat", "wb");
    }

    fine_burst_measure_impl::~fine_burst_measure_impl() {
      volk_free(d_tmp_c);
      volk_free(d_phase_diff_f);
      volk_free(d_magnitude_f);
      volk_free(d_magnitude_shifted_f);
      volk_free(d_magnitude_smooth_f);
      volk_free(d_window_f);
      delete d_fft;
      fclose(fid);
      fclose(fid2);
    }

    void fine_burst_measure_impl::resize_bins(size_t numBins) {
      if (d_histogram.size() < numBins) {
        d_histogram.resize(numBins);
        d_histogram_smooth.resize(numBins);
      }
    }

    fsk_info fine_burst_measure_impl::findFreqOffset(float* input, size_t length, const std::vector<float>& freqs) {
      // Cluster each point to its nearest freq
      /*printf("freqs: ");
      for (size_t i = 0; i < freqs.size(); i++) {
        printf("%f ", freqs[i]);
      }
      printf("\n");*/
      std::vector<std::vector<float> > data;
      data.resize(freqs.size());
      for (size_t i = 0; i < length; i++) {
        size_t bestMatch = 0;
        float bestDiff = std::abs(input[i] - freqs[0]);
        for (size_t j = 1; j < freqs.size(); j++) {
          if (std::abs(input[i] - freqs[j]) < bestDiff) {
            bestMatch = j;
            bestDiff = std::abs(input[i] - freqs[j]);
          }
        }
        if (bestDiff < .1)
        data[bestMatch].push_back(input[i]);
      }

      // Find the average of each bestMatch
      std::vector<float> fcs;
      for (size_t i = 0; i < freqs.size(); i++) {
        float sum = 0;
        for (size_t j = 0; j < data[i].size(); j++) {
          sum += data[i][j];
        }        
        //printf("peak %zu: %f, %zu\n", i, sum, data[i].size());
        if (data[i].size() > 0) fcs.push_back(sum / data[i].size());
      }

      // Find the average of the fcs
      float shift = 0;
      for (size_t i = 0; i < fcs.size(); i++) {
        shift += fcs[i];
      }
      if (fcs.size() > 0) shift /= fcs.size();
      fsk_info info;
      info.shift = shift;
      for (size_t i = 0; i < fcs.size(); i++) {
        info.freqs.push_back(fcs[i] - shift);
      }
      
      //printf("shift = %f\n", shift);
      return info;
        
    }

    float fine_burst_measure_impl::findPower(const gr_complex* in, const std::vector<float>& freqs, size_t inSize) {
      size_t cpySize = std::min(inSize, d_fft_size);
      volk_32fc_32f_multiply_32fc(d_fft->get_inbuf(), in, d_window_f, cpySize);
      //memcpy(d_fft->get_inbuf(), in, cpySize*sizeof(gr_complex));
      memset(&d_fft->get_inbuf()[cpySize], 0, sizeof(gr_complex)*(d_fft_size - cpySize));
      d_fft->execute();
      volk_32fc_magnitude_squared_32f(d_magnitude_f, d_fft->get_outbuf(), d_fft_size);
      memcpy(&d_magnitude_shifted_f[0], &d_magnitude_f[d_fft_size/2], sizeof(float) * d_fft_size/2);
      memcpy(&d_magnitude_shifted_f[d_fft_size/2], &d_magnitude_f[0], sizeof(float) * d_fft_size/2);
      
      float power = 0;
      for (size_t i = 0; i < d_fft_size; i++) {
        power += d_magnitude_shifted_f[i];
      }
      return power / inSize;
    }

    float fine_burst_measure_impl::findOccupiedBW(float power, size_t inSize) {
      // Search for a 99.5% occupied bandwidth
      power *= inSize;
      float thresh99 = power*.005;
      float sum99 = 0;
      float bw99 = 0;
      for (size_t i = 0; i < d_fft_size / 2; i++) {
        if (sum99 < thresh99) sum99 += d_magnitude_shifted_f[i] + d_magnitude_shifted_f[d_fft_size - i - 1];
        else {
          bw99 = (d_fft_size - 2*i) / float(d_fft_size);
          break;
        }
      }
      return bw99;

    }
      


    // What params do I need in order to run this???
    // For each burst, I need to know sample_rate, freq_offset.
    // Let's allow for different input sample_rates because it isn't hard to do.
    // The only thing we need is to specify is the desired bin width in Hz.

    // The histogram will go from -fs/2 to fs/2.  How many bins should there be? 
    // I am expected to be at about 4x the symbol rate.

    void fine_burst_measure_impl::pdu_handler(pmt::pmt_t pdu) {
      if (!pmt::is_pair(pdu)) {
        GR_LOG_WARN(d_logger, "PDU not a pair, dropping");
        return;
      }

      pmt::pmt_t metadata = pmt::car(pdu);
      pmt::pmt_t data = pmt::cdr(pdu);

      if (!pmt::is_dict(metadata)) {
        GR_LOG_WARN(d_logger, "PDU metadata not a dictionary, dropping");
        return;
      }

      float center_frequency = pmt::to_float(pmt::dict_ref(metadata, pmt::mp("center_frequency"), pmt::PMT_NIL));
      float relative_frequency = pmt::to_float(pmt::dict_ref(metadata, pmt::mp("relative_frequency"), pmt::PMT_NIL));
      float sample_rate = pmt::to_float(pmt::dict_ref(metadata, pmt::mp("sample_rate"), pmt::PMT_NIL));

      if (!pmt::is_c32vector(data)) {
        GR_LOG_WARN(d_logger, "PDU data not complex, dropping");
        return;
      }

      size_t num_bins = std::ceil(sample_rate / d_freq_resolution);
      resize_bins(num_bins);

      size_t burst_size = pmt::length(data);
      const gr_complex* burst = (const gr_complex*) pmt::c32vector_elements(data, burst_size);
      int n_samps(std::min(uint32_t(burst_size*.8), d_analyze_samples));
      int start = (burst_size - n_samps) / 2;
          
      // Create a histogram of phase differences from the middle of the data
      volk_32fc_x2_multiply_conjugate_32fc(&d_tmp_c[0], &burst[start+1], &burst[start], n_samps);
      float maxV = 0;
      for(int i = 0; i < n_samps; i++) {
        d_phase_diff_f[i] = gr::fast_atan2f(std::imag(d_tmp_c[i]), std::real(d_tmp_c[i]));
      }
      // atan2 returns a number between -pi and pi, so we have our min and max freq
      slow_histogram_unbounded(d_phase_diff_f, n_samps, d_histogram.data(), d_min_freq, d_max_freq, num_bins);

      uint64_t bid = pmt::to_uint64(pmt::dict_ref(metadata, pmt::mp("burst_id"), pmt::PMT_NIL));
      // Do a moving average on the histogram to smooth it.
      for (size_t i = 1; i < num_bins-1; i++) {
        d_histogram_smooth[i] = .25*(d_histogram[i-1] + d_histogram[i+1]) + .5*d_histogram[i];
      }
      d_histogram_smooth[0] = d_histogram[0];
      d_histogram_smooth[num_bins-1] = d_histogram[num_bins-1];

      //Search for any peaks.
      d_peak_list.clear();
      uint32_t* hh = &d_histogram_smooth[0];
      for (size_t i = 1; i < num_bins-1; i++) {
        if (hh[i] > hh[i-1] && hh[i] > hh[i+1]) {
          // We have a potential peak.
          d_peak_list.push_back({i, hh[i]});
        }
      }
      std::sort(d_peak_list.begin(), d_peak_list.end(),
          [] (const peak &a, const peak&b) {
          return a.value > b.value;
      });
      // We are looking for 2 FSK
      // Our first peak is our best.  Don't consider any peak with a height less than 1/10
      if (d_peak_list.size() > 0) {
        float thresh = d_peak_list[0].value * .1;
        for (size_t i = 1; i < d_peak_list.size(); i++) {
          if (d_peak_list[i].value < thresh) {
            d_peak_list.resize(i);
          }
        }
        //printf("%zu: peaks = %zu\n", bid, d_peak_list.size());
        float bestScore = 0;
        size_t bestIndex = 0;
        fsk_info info;
        float bandwidth;
        if (d_peak_list.size() == 1) {
          info = findFreqOffset(d_phase_diff_f, n_samps, {bin2Freq(d_peak_list[0].center_bin, num_bins)});
        } else if (d_peak_list.size() > 1) {
          // Find the best second peak
          for (size_t i = 1; i < d_peak_list.size(); i++) {
            float score = std::abs(d_peak_list[0].center_bin - d_peak_list[i].center_bin) * d_peak_list[i].value;
            if (score > bestScore) {
              bestScore = score;
              bestIndex = i;
            }
          }
          info = findFreqOffset(d_phase_diff_f, n_samps, {bin2Freq(d_peak_list[0].center_bin, num_bins), bin2Freq(d_peak_list[bestIndex].center_bin, num_bins)});
        }
        //printf("%zu: %f\n", bid, info.shift);

        // Apply fine frequency correction
        d_r.set_phase_incr(exp(gr_complex(0, -info.shift)));
        d_r.set_phase(gr_complex(1, 0));
        d_tmp_burst.resize(burst_size);
        d_r.rotateN(&d_tmp_burst[0], burst, burst_size);
        center_frequency += (info.shift / M_PI / 2.0) * sample_rate;
        relative_frequency += (info.shift / M_PI / 2.0) * sample_rate;

        float power = findPower(burst + start, info.freqs, n_samps);
        bandwidth = findOccupiedBW(power, n_samps);
        
        pmt::pmt_t pdu_vector = pmt::init_c32vector(burst_size, d_tmp_burst);
        metadata = pmt::dict_add(metadata, pmt::mp("center_frequency"), pmt::mp(center_frequency));
        metadata = pmt::dict_add(metadata, pmt::mp("relative_frequency"), pmt::mp(relative_frequency));
        metadata = pmt::dict_add(metadata, pmt::mp("bandwidth"), pmt::from_float(bandwidth * sample_rate));
        metadata = pmt::dict_add(metadata, pmt::mp("burst_power"), pmt::from_float(10*log10(power)));
        if (info.freqs.size() == 2) {
          //metadata = pmt::dict_add(metadata, pmt::mp("mark_freq"), pmt::from_float(std::max(info.freqs[0], info.freqs[1]) * sample_rate + center_frequency));
          //metadata = pmt::dict_add(metadata, pmt::mp("space_freq"), pmt::from_float(std::min(info.freqs[0], info.freqs[1]) * sample_rate + center_frequency));
          metadata = pmt::dict_add(metadata, pmt::mp("mark_freq"), pmt::from_float(std::max(info.freqs[0], info.freqs[1]) * sample_rate/M_PI + relative_frequency));
          metadata = pmt::dict_add(metadata, pmt::mp("space_freq"), pmt::from_float(std::min(info.freqs[0], info.freqs[1]) * sample_rate/M_PI + relative_frequency));
        }
          
        
        //std::cout << metadata << std::endl;
        
        pmt::pmt_t out_msg = pmt::cons(metadata, pdu_vector);
        message_port_pub(pmt::mp("pdu_out"), out_msg);
      
        
      }
      

      //printf("%zu: Num peaks = %zu, start = %zu, n_samps = %zu, v = %f\n", bid, d_peak_list.size(), start, n_samps, d_magnitude_smooth_f[1024]);
      /*fwrite(&bid, sizeof(bid), 1, fid);
      fwrite(&burst_size, sizeof(burst_size), 1, fid);
      fwrite(&d_tmp_burst[0], sizeof(burst[0]), burst_size, fid);
      //fwrite(&burst[0], sizeof(burst[0]), burst_size, fid);
      fwrite(&bid, sizeof(bid), 1, fid2);
      //fwrite(&d_fft_size, sizeof(d_fft_size), 1, fid2);
      //fwrite(&d_magnitude_smooth_f[0], 4, d_fft_size, fid2);
      fwrite(&num_bins, sizeof(num_bins), 1, fid2);
      fwrite(&d_histogram_smooth[0], 4, num_bins, fid2);*/

      /* output pmt needs to have the following fields
          start_time   - fbt
          duration     - fbt
          absolute center freq  - fbm
          relative center freq (to collection center) - fbm
          sample rate    - tbtp
          bandwidth/mark to space diff
          baud rate - pcr
          snr/evm - pcr
      */

    }

  }
}
      
      
      
