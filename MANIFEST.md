title: The FHSS_UTILS OOT Module
brief: Tools for processing frequency hopping spread spectrum signals,
tags: # Tags are arbitrary, but look at CGRAN what other authors are using
  - FHSS
author:
  - Jacob Gilbert
copyright_owner:
  - Jacob Gilbert
license:
gr_supported_version: v3.7, v3.8
#repo: # Put the URL of the repository here, or leave blank for default
#website: <module_website> # If you have a separate project website, put it here
#icon: <icon_url> # Put a URL to a square image here that will be used as an icon on CGRAN
---
This GNU Radio module contains tools for processing frequency hopping spread spectrum signals, as well as known and unknown bursty signals in general. Blocks derived from the gr-iridium project exist to detect narrowband bursts within wideband signals and downconvert and center them. Metadata is tracked through this process enabling reconstruction of where the bursts originated in time and frequency.
