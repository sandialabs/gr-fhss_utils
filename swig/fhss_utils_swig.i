/* -*- c++ -*- */

#define FHSS_UTILS_API

%include "gnuradio.i"           // the common stuff

//load generated python docstrings
%include "fhss_utils_swig_doc.i"

%{
#include "fhss_utils/cf_estimate.h"
#include "fhss_utils/fft_burst_tagger.h"
#include "fhss_utils/tagged_burst_to_pdu.h"
%}

%include "fhss_utils/cf_estimate.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, cf_estimate);
%include "fhss_utils/fft_burst_tagger.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, fft_burst_tagger);

%include "fhss_utils/tagged_burst_to_pdu.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, tagged_burst_to_pdu);
