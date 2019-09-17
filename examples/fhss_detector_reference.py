#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Burst Detector Reference
# Author: J. Gilbert
# Description: Reference Application
# GNU Radio version: 3.7.13.5
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from PyQt4 import Qt
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import qtgui
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.qtgui import Range, RangeWidget
from optparse import OptionParser
import fhss_utils
import pdu_utils
import sip
import sys
import time
from gnuradio import qtgui


class fhss_detector_reference(gr.top_block, Qt.QWidget):

    def __init__(self, burst_width=int(500e3), cfo_start_offset=0, cfo_threshold=0.5, cfo_time_to_average=0.0005, decimation=32, fft_size=256, hist_time=0.004, lookahead_time=0.0005, max_burst_time=0.5, min_burst_time=0.001, output_attenuation=40, output_cutoff=0.5, output_trans_width=0.05, post_burst_time=0.00008, pre_burst_time=0.00008, samp_rate=int(30.72e6)):
        gr.top_block.__init__(self, "Burst Detector Reference")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Burst Detector Reference")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "fhss_detector_reference")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())


        ##################################################
        # Parameters
        ##################################################
        self.burst_width = burst_width
        self.cfo_start_offset = cfo_start_offset
        self.cfo_threshold = cfo_threshold
        self.cfo_time_to_average = cfo_time_to_average
        self.decimation = decimation
        self.fft_size = fft_size
        self.hist_time = hist_time
        self.lookahead_time = lookahead_time
        self.max_burst_time = max_burst_time
        self.min_burst_time = min_burst_time
        self.output_attenuation = output_attenuation
        self.output_cutoff = output_cutoff
        self.output_trans_width = output_trans_width
        self.post_burst_time = post_burst_time
        self.pre_burst_time = pre_burst_time
        self.samp_rate = samp_rate

        ##################################################
        # Variables
        ##################################################
        self.threshold = threshold = 10
        self.gain = gain = 40
        self.fir_taps = fir_taps = firdes.low_pass_2(1, 1, 0.5*output_cutoff, 0.5*output_trans_width, output_attenuation)
        self.decim_taps = decim_taps = firdes.low_pass_2(1, 1, output_cutoff/decimation, output_trans_width/decimation, output_attenuation)
        self.center_freq = center_freq = 915e6

        ##################################################
        # Blocks
        ##################################################
        self._threshold_range = Range(6, 25, 1, 10, 200)
        self._threshold_win = RangeWidget(self._threshold_range, self.set_threshold, 'Threshold', "counter_slider", float)
        self.top_grid_layout.addWidget(self._threshold_win, 2, 2, 1, 1)
        for r in range(2, 3):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(2, 3):
            self.top_grid_layout.setColumnStretch(c, 1)
        self._gain_range = Range(0, 70, 1, 40, 200)
        self._gain_win = RangeWidget(self._gain_range, self.set_gain, 'Gain', "counter_slider", float)
        self.top_grid_layout.addWidget(self._gain_win, 2, 1, 1, 1)
        for r in range(2, 3):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(1, 2):
            self.top_grid_layout.setColumnStretch(c, 1)
        self._center_freq_range = Range(100e6, 2000e6, 1e6, 915e6, 200)
        self._center_freq_win = RangeWidget(self._center_freq_range, self.set_center_freq, 'Center Frequency', "counter_slider", float)
        self.top_grid_layout.addWidget(self._center_freq_win, 2, 0, 1, 1)
        for r in range(2, 3):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 1):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.uhd_usrp_source_1 = uhd.usrp_source(
        	",".join(("", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_source_1.set_samp_rate(samp_rate)
        self.uhd_usrp_source_1.set_center_freq(center_freq, 0)
        self.uhd_usrp_source_1.set_gain(gain, 0)
        self.uhd_usrp_source_1.set_antenna('TX/RX', 0)
        self.uhd_usrp_source_1.set_auto_dc_offset(True, 0)
        self.uhd_usrp_source_1.set_auto_iq_balance(True, 0)
        (self.uhd_usrp_source_1).set_min_output_buffer(2048000)
        self.qtgui_waterfall_sink_x_0_0 = qtgui.waterfall_sink_c(
        	1024, #size
        	firdes.WIN_BLACKMAN_hARRIS, #wintype
        	0, #fc
        	samp_rate, #bw
        	"Input Spectrogram", #name
                1 #number of inputs
        )
        self.qtgui_waterfall_sink_x_0_0.set_update_time(0.01)
        self.qtgui_waterfall_sink_x_0_0.enable_grid(False)
        self.qtgui_waterfall_sink_x_0_0.enable_axis_labels(True)

        if not True:
          self.qtgui_waterfall_sink_x_0_0.disable_legend()

        if "complex" == "float" or "complex" == "msg_float":
          self.qtgui_waterfall_sink_x_0_0.set_plot_pos_half(not True)

        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        colors = [0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]
        for i in xrange(1):
            if len(labels[i]) == 0:
                self.qtgui_waterfall_sink_x_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_waterfall_sink_x_0_0.set_line_label(i, labels[i])
            self.qtgui_waterfall_sink_x_0_0.set_color_map(i, colors[i])
            self.qtgui_waterfall_sink_x_0_0.set_line_alpha(i, alphas[i])

        self.qtgui_waterfall_sink_x_0_0.set_intensity_range(-120, 0)

        self._qtgui_waterfall_sink_x_0_0_win = sip.wrapinstance(self.qtgui_waterfall_sink_x_0_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_waterfall_sink_x_0_0_win, 0, 0, 1, 3)
        for r in range(0, 1):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 3):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.qtgui_waterfall_sink_x_0 = qtgui.waterfall_sink_c(
        	128, #size
        	firdes.WIN_BLACKMAN_hARRIS, #wintype
        	0, #fc
        	samp_rate/decimation/2, #bw
        	"Burst Spectrogram", #name
                0 #number of inputs
        )
        self.qtgui_waterfall_sink_x_0.set_update_time(0.005)
        self.qtgui_waterfall_sink_x_0.enable_grid(False)
        self.qtgui_waterfall_sink_x_0.enable_axis_labels(True)

        if not True:
          self.qtgui_waterfall_sink_x_0.disable_legend()

        if "msg_complex" == "float" or "msg_complex" == "msg_float":
          self.qtgui_waterfall_sink_x_0.set_plot_pos_half(not True)

        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        colors = [0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]
        for i in xrange(1):
            if len(labels[i]) == 0:
                self.qtgui_waterfall_sink_x_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_waterfall_sink_x_0.set_line_label(i, labels[i])
            self.qtgui_waterfall_sink_x_0.set_color_map(i, colors[i])
            self.qtgui_waterfall_sink_x_0.set_line_alpha(i, alphas[i])

        self.qtgui_waterfall_sink_x_0.set_intensity_range(-140, 10)

        self._qtgui_waterfall_sink_x_0_win = sip.wrapinstance(self.qtgui_waterfall_sink_x_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_waterfall_sink_x_0_win, 1, 0, 1, 1)
        for r in range(1, 2):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 1):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0 = qtgui.time_sink_f(
        	102400, #size
        	1, #samp_rate
        	"Soft Symbols", #name
        	0 #number of inputs
        )
        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_update_time(0.01)
        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_y_axis(-1.5, 1.5)

        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.enable_tags(-1, True)
        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.000001, .001, 0, "")
        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.enable_autoscale(False)
        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.enable_grid(True)
        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.enable_stem_plot(False)

        if not True:
          self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.disable_legend()

        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["black", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "blue"]
        styles = [0, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        markers = [0, -1, -1, -1, -1,
                   -1, -1, -1, -1, -1]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]

        for i in xrange(1):
            if len(labels[i]) == 0:
                self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_0_0_0_0_0_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0_0_0_0_0_0_0_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_time_sink_x_0_0_0_0_0_0_0_0_win, 1, 2, 1, 1)
        for r in range(1, 2):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(2, 3):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.qtgui_time_sink_x_0_0_0_0_0 = qtgui.time_sink_f(
        	102400, #size
        	samp_rate/decimation, #samp_rate
        	"FM Demodulation", #name
        	0 #number of inputs
        )
        self.qtgui_time_sink_x_0_0_0_0_0.set_update_time(0.01)
        self.qtgui_time_sink_x_0_0_0_0_0.set_y_axis(-2, 2)

        self.qtgui_time_sink_x_0_0_0_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_0_0_0_0_0.enable_tags(-1, True)
        self.qtgui_time_sink_x_0_0_0_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.000001, .001, 0, "")
        self.qtgui_time_sink_x_0_0_0_0_0.enable_autoscale(False)
        self.qtgui_time_sink_x_0_0_0_0_0.enable_grid(False)
        self.qtgui_time_sink_x_0_0_0_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0_0_0_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_0_0_0_0_0.enable_stem_plot(False)

        if not True:
          self.qtgui_time_sink_x_0_0_0_0_0.disable_legend()

        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "blue"]
        styles = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
                   -1, -1, -1, -1, -1]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]

        for i in xrange(1):
            if len(labels[i]) == 0:
                self.qtgui_time_sink_x_0_0_0_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_time_sink_x_0_0_0_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0_0_0_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0_0_0_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0_0_0_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0_0_0_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0_0_0_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_0_0_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0_0_0_0_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_time_sink_x_0_0_0_0_0_win, 1, 1, 1, 1)
        for r in range(1, 2):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(1, 2):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.pdu_utils_pdu_split_0_0 = pdu_utils.pdu_split(False)
        self.pdu_utils_pdu_fir_filter_1 = pdu_utils.pdu_fir_filter(2, (fir_taps))
        self.pdu_utils_pdu_fine_time_measure_0 = pdu_utils.pdu_fine_time_measure(pre_burst_time, post_burst_time, 10, 15)
        self.pdu_utils_pdu_clock_recovery_0_0 = pdu_utils.pdu_clock_recovery(True)
        self.pdu_utils_pdu_clock_recovery_0 = pdu_utils.pdu_clock_recovery(False)
        self.fhss_utils_tagged_burst_to_pdu_0 = fhss_utils.tagged_burst_to_pdu(decimation, (decim_taps), min_burst_time, max_burst_time, 0.0, 1.0, 1.0, samp_rate, 7)
        self.fhss_utils_pdu_quadrature_demod_cf_0 = fhss_utils.pdu_quadrature_demod_cf(([1]), 1)
        self.fhss_utils_fine_burst_measure_0 = fhss_utils.fine_burst_measure(4000, 2048, 0.5)
        self.fhss_utils_fft_burst_tagger_0 = fhss_utils.fft_burst_tagger(center_freq, fft_size, samp_rate, int(round((float(samp_rate)/fft_size)*pre_burst_time)), int(round((float(samp_rate)/fft_size)*post_burst_time)), burst_width, 0, 0, threshold, int(round((float(samp_rate)/fft_size)*hist_time)), int(round((float(samp_rate)/fft_size)*lookahead_time)), False)
        (self.fhss_utils_fft_burst_tagger_0).set_min_output_buffer(2048000)
        self.blocks_message_debug_0_1 = blocks.message_debug()



        ##################################################
        # Connections
        ##################################################
        self.msg_connect((self.fhss_utils_fine_burst_measure_0, 'pdu_out'), (self.pdu_utils_pdu_fir_filter_1, 'pdu_in'))
        self.msg_connect((self.fhss_utils_pdu_quadrature_demod_cf_0, 'fpdus'), (self.pdu_utils_pdu_clock_recovery_0, 'pdu_in'))
        self.msg_connect((self.fhss_utils_pdu_quadrature_demod_cf_0, 'fpdus'), (self.pdu_utils_pdu_clock_recovery_0_0, 'pdu_in'))
        self.msg_connect((self.fhss_utils_pdu_quadrature_demod_cf_0, 'fpdus'), (self.qtgui_time_sink_x_0_0_0_0_0, 'in'))
        self.msg_connect((self.fhss_utils_tagged_burst_to_pdu_0, 'cpdus'), (self.fhss_utils_fine_burst_measure_0, 'pdu_in'))
        self.msg_connect((self.pdu_utils_pdu_clock_recovery_0, 'pdu_out'), (self.pdu_utils_pdu_split_0_0, 'pdu_in'))
        self.msg_connect((self.pdu_utils_pdu_clock_recovery_0, 'pdu_out'), (self.qtgui_time_sink_x_0_0_0_0_0_0_0_0, 'in'))
        self.msg_connect((self.pdu_utils_pdu_fine_time_measure_0, 'pdu_out'), (self.fhss_utils_pdu_quadrature_demod_cf_0, 'cpdus'))
        self.msg_connect((self.pdu_utils_pdu_fine_time_measure_0, 'pdu_out'), (self.qtgui_waterfall_sink_x_0, 'in'))
        self.msg_connect((self.pdu_utils_pdu_fir_filter_1, 'pdu_out'), (self.pdu_utils_pdu_fine_time_measure_0, 'pdu_in'))
        self.msg_connect((self.pdu_utils_pdu_split_0_0, 'dict'), (self.blocks_message_debug_0_1, 'print'))
        self.connect((self.fhss_utils_fft_burst_tagger_0, 0), (self.fhss_utils_tagged_burst_to_pdu_0, 0))
        self.connect((self.uhd_usrp_source_1, 0), (self.fhss_utils_fft_burst_tagger_0, 0))
        self.connect((self.uhd_usrp_source_1, 0), (self.qtgui_waterfall_sink_x_0_0, 0))

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "fhss_detector_reference")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_burst_width(self):
        return self.burst_width

    def set_burst_width(self, burst_width):
        self.burst_width = burst_width

    def get_cfo_start_offset(self):
        return self.cfo_start_offset

    def set_cfo_start_offset(self, cfo_start_offset):
        self.cfo_start_offset = cfo_start_offset

    def get_cfo_threshold(self):
        return self.cfo_threshold

    def set_cfo_threshold(self, cfo_threshold):
        self.cfo_threshold = cfo_threshold

    def get_cfo_time_to_average(self):
        return self.cfo_time_to_average

    def set_cfo_time_to_average(self, cfo_time_to_average):
        self.cfo_time_to_average = cfo_time_to_average

    def get_decimation(self):
        return self.decimation

    def set_decimation(self, decimation):
        self.decimation = decimation
        self.set_decim_taps(firdes.low_pass_2(1, 1, self.output_cutoff/self.decimation, self.output_trans_width/self.decimation, self.output_attenuation))
        self.qtgui_waterfall_sink_x_0.set_frequency_range(0, self.samp_rate/self.decimation/2)
        self.qtgui_time_sink_x_0_0_0_0_0.set_samp_rate(self.samp_rate/self.decimation)

    def get_fft_size(self):
        return self.fft_size

    def set_fft_size(self, fft_size):
        self.fft_size = fft_size

    def get_hist_time(self):
        return self.hist_time

    def set_hist_time(self, hist_time):
        self.hist_time = hist_time

    def get_lookahead_time(self):
        return self.lookahead_time

    def set_lookahead_time(self, lookahead_time):
        self.lookahead_time = lookahead_time

    def get_max_burst_time(self):
        return self.max_burst_time

    def set_max_burst_time(self, max_burst_time):
        self.max_burst_time = max_burst_time

    def get_min_burst_time(self):
        return self.min_burst_time

    def set_min_burst_time(self, min_burst_time):
        self.min_burst_time = min_burst_time

    def get_output_attenuation(self):
        return self.output_attenuation

    def set_output_attenuation(self, output_attenuation):
        self.output_attenuation = output_attenuation
        self.set_fir_taps(firdes.low_pass_2(1, 1, 0.5*self.output_cutoff, 0.5*self.output_trans_width, self.output_attenuation))
        self.set_decim_taps(firdes.low_pass_2(1, 1, self.output_cutoff/self.decimation, self.output_trans_width/self.decimation, self.output_attenuation))

    def get_output_cutoff(self):
        return self.output_cutoff

    def set_output_cutoff(self, output_cutoff):
        self.output_cutoff = output_cutoff
        self.set_fir_taps(firdes.low_pass_2(1, 1, 0.5*self.output_cutoff, 0.5*self.output_trans_width, self.output_attenuation))
        self.set_decim_taps(firdes.low_pass_2(1, 1, self.output_cutoff/self.decimation, self.output_trans_width/self.decimation, self.output_attenuation))

    def get_output_trans_width(self):
        return self.output_trans_width

    def set_output_trans_width(self, output_trans_width):
        self.output_trans_width = output_trans_width
        self.set_fir_taps(firdes.low_pass_2(1, 1, 0.5*self.output_cutoff, 0.5*self.output_trans_width, self.output_attenuation))
        self.set_decim_taps(firdes.low_pass_2(1, 1, self.output_cutoff/self.decimation, self.output_trans_width/self.decimation, self.output_attenuation))

    def get_post_burst_time(self):
        return self.post_burst_time

    def set_post_burst_time(self, post_burst_time):
        self.post_burst_time = post_burst_time

    def get_pre_burst_time(self):
        return self.pre_burst_time

    def set_pre_burst_time(self, pre_burst_time):
        self.pre_burst_time = pre_burst_time

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_source_1.set_samp_rate(self.samp_rate)
        self.qtgui_waterfall_sink_x_0_0.set_frequency_range(0, self.samp_rate)
        self.qtgui_waterfall_sink_x_0.set_frequency_range(0, self.samp_rate/self.decimation/2)
        self.qtgui_time_sink_x_0_0_0_0_0.set_samp_rate(self.samp_rate/self.decimation)

    def get_threshold(self):
        return self.threshold

    def set_threshold(self, threshold):
        self.threshold = threshold

    def get_gain(self):
        return self.gain

    def set_gain(self, gain):
        self.gain = gain
        self.uhd_usrp_source_1.set_gain(self.gain, 0)


    def get_fir_taps(self):
        return self.fir_taps

    def set_fir_taps(self, fir_taps):
        self.fir_taps = fir_taps
        self.pdu_utils_pdu_fir_filter_1.set_taps((self.fir_taps))

    def get_decim_taps(self):
        return self.decim_taps

    def set_decim_taps(self, decim_taps):
        self.decim_taps = decim_taps

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.uhd_usrp_source_1.set_center_freq(self.center_freq, 0)


def argument_parser():
    description = 'Reference Application'
    parser = OptionParser(usage="%prog: [options]", option_class=eng_option, description=description)
    parser.add_option(
        "", "--burst-width", dest="burst_width", type="intx", default=int(500e3),
        help="Set Burst Width [Hz] [default=%default]")
    parser.add_option(
        "", "--cfo-start-offset", dest="cfo_start_offset", type="eng_float", default=eng_notation.num_to_str(0),
        help="Set CFO Start Offset [s] [default=%default]")
    parser.add_option(
        "", "--cfo-threshold", dest="cfo_threshold", type="eng_float", default=eng_notation.num_to_str(0.5),
        help="Set CFO threshold [default=%default]")
    parser.add_option(
        "", "--cfo-time-to-average", dest="cfo_time_to_average", type="eng_float", default=eng_notation.num_to_str(0.0005),
        help="Set CFO Average Time [s] [default=%default]")
    parser.add_option(
        "", "--decimation", dest="decimation", type="intx", default=32,
        help="Set Decimation [default=%default]")
    parser.add_option(
        "", "--fft-size", dest="fft_size", type="intx", default=256,
        help="Set FFT Size [default=%default]")
    parser.add_option(
        "", "--hist-time", dest="hist_time", type="eng_float", default=eng_notation.num_to_str(0.004),
        help="Set History Time [s] [default=%default]")
    parser.add_option(
        "", "--lookahead-time", dest="lookahead_time", type="eng_float", default=eng_notation.num_to_str(0.0005),
        help="Set Lookahead Time [s] [default=%default]")
    parser.add_option(
        "", "--max-burst-time", dest="max_burst_time", type="eng_float", default=eng_notation.num_to_str(0.5),
        help="Set Max Burst Time [s] [default=%default]")
    parser.add_option(
        "", "--min-burst-time", dest="min_burst_time", type="eng_float", default=eng_notation.num_to_str(0.001),
        help="Set Min Burst Time [s] [default=%default]")
    parser.add_option(
        "", "--output-attenuation", dest="output_attenuation", type="eng_float", default=eng_notation.num_to_str(40),
        help="Set Output Attenuation [default=%default]")
    parser.add_option(
        "", "--output-cutoff", dest="output_cutoff", type="eng_float", default=eng_notation.num_to_str(0.5),
        help="Set Output Cutoff [cycles/samp] [default=%default]")
    parser.add_option(
        "", "--output-trans-width", dest="output_trans_width", type="eng_float", default=eng_notation.num_to_str(0.05),
        help="Set Output Trans. Width [cycles/samp] [default=%default]")
    parser.add_option(
        "", "--post-burst-time", dest="post_burst_time", type="eng_float", default=eng_notation.num_to_str(0.00008),
        help="Set Post Burst Time [s] [default=%default]")
    parser.add_option(
        "", "--pre-burst-time", dest="pre_burst_time", type="eng_float", default=eng_notation.num_to_str(0.00008),
        help="Set Pre Burst Time [s] [default=%default]")
    parser.add_option(
        "-r", "--samp-rate", dest="samp_rate", type="intx", default=int(30.72e6),
        help="Set Sample Rate [default=%default]")
    return parser


def main(top_block_cls=fhss_detector_reference, options=None):
    if options is None:
        options, _ = argument_parser().parse_args()

    from distutils.version import StrictVersion
    if StrictVersion(Qt.qVersion()) >= StrictVersion("4.5.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls(burst_width=options.burst_width, cfo_start_offset=options.cfo_start_offset, cfo_threshold=options.cfo_threshold, cfo_time_to_average=options.cfo_time_to_average, decimation=options.decimation, fft_size=options.fft_size, hist_time=options.hist_time, lookahead_time=options.lookahead_time, max_burst_time=options.max_burst_time, min_burst_time=options.min_burst_time, output_attenuation=options.output_attenuation, output_cutoff=options.output_cutoff, output_trans_width=options.output_trans_width, post_burst_time=options.post_burst_time, pre_burst_time=options.pre_burst_time, samp_rate=options.samp_rate)
    tb.start()
    tb.show()

    def quitting():
        tb.stop()
        tb.wait()
    qapp.connect(qapp, Qt.SIGNAL("aboutToQuit()"), quitting)
    qapp.exec_()


if __name__ == '__main__':
    main()
