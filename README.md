# NuView

Browser-based event display for liquid argon TPC neutrino data, starting with
the MicroBooNE public samples. Part of the Auraeon project family:
interactive visualizations of physical systems.

## Quick start

1. Download a MicroBooNE open-sample HDF5 file (NoWire variant is sufficient):
   https://doi.org/10.5281/zenodo.7261921
2. Extract an event:
   `python extract_event.py nue_NoWire_18.h5 1299 event_1299.json`
3. Open `index.html` in a browser and drop the JSON onto it.

Three wire-plane views (wire number vs drift time), hits colored by
deposited charge. No build step, no dependencies beyond h5py + numpy.

## Roadmap

- Truth mode: color hits by true particle via edep_table join
- Reco mode: color by Pandora PFP / slice
- Event browser across a file; 3D truth-trajectory companion view

## Acknowledgment

We acknowledge the MicroBooNE Collaboration for making publicly available the
data sets employed in this work: DOI 10.5281/zenodo.7261921 (CC-BY 4.0).
See https://github.com/uboone/OpenSamples for official documentation.
Code: MIT license.
