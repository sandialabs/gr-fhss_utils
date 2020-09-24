![snl](docs/figures/snl.png "Sandia National Laboratories")

## GNU Radio FHSS Utilities

This GNU Radio module contains tools for processing frequency hopping spread spectrum signals. Blocks derived from the gr-iridium project exist to detect narrowband bursts within wideband signals and downconvert and center them. Metadata is tracked through this process enabling reconstruction of where the bursts originated in time and frequency. Another set of blocks exists to baseband all bursts within a high fidelity signal capture which is useful for reverse engineering of FHSS datasets.

---

### General Concept of High Fidelity FHSS Signal Dehopping

The dataset dehopper blocks were designed to quickly allow for good accuracy dehopping of high fidelity FHSS FSK recordings. This is accomplished by a two-stage dehopping process by which a coarse FFT is taken and peak values are taken by a simple sample-and-hold block when an amplitude threshold is crossed, then a second stage does fine frequency correction by taking an instantaneous frequency average. This works well for FSK signals but requires some work for other signals.

![fhss_dehopper.png](docs/figures/fhss_dehopper.png "Simple FHSS Dehopping Flowgraph")

This has been implemented without use of OOT DSP, though the module has python hier blocks that are installed with the build. Alternately there are GRC hier blocks located in the _gr-fhss_utils/examples/hier_blocks/_ directory

---

An overview of the Sandia National Laboratories Utilities GR Modules can be found in the README for https://github.com/sandialabs/gr-pdu_utils
