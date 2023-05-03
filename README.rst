JWST/MIRI Imaging/Coronagraphy Flux Calibration
===============================================

Code to measure flux calibration stars and compute the
flux calibration for JWST/MIRI Imaging and Coronagraphic filters.

Main goal is to fully automate this work to allow for straightforward
recalculations from the raw data to the final calibration factors for each
filter.  This enables tests of the sensitivity of the flux calibration
to changes in the jwst pipeline or any of the steps in the flux calibration
factor calculation.

Contributors
------------
Karl Gordon

License
-------

This code is licensed under a 3-clause BSD style license (see the
``LICENSE`` file).

Data
----

Observations are organized into three main directories for the three types
of stars (HotStars, ADwarfs, and SolarAnalogs).  In each directory, subdirectories
for each filter (e.g., F1500W) have the raw (uncal) data.  The automated
reductions done in step 1 below create subdirectories in the filter directories
for each observation for each star (set1, set2, ...).

Details
-------

1. Reduce data with jwst pipeline: pipeline_one_filter.py.
   Runs all three stages of pipeline for a one type of star and one filter.
   There are some differences from the defaults for MIRI.
   `pipeline_all` runs everything for all MIRI imaging/coronagraphic filters.

2. Compute WebbPSFs: calc_webbpsfs.py.
   Runs WebbPSFs with parameters for MIRI

3. Compute aperture corrections: plot_many_apertures.py.
   Uses observations of a bright stars and PSFs from WebbPSFs to compute aperture
   corrections.  Uses the bright star observations for the inner region and
   the WebbPSFs for the outer regions.
   `calc_all_ees` runs all MIRI imaging filters for BD+60 1753 observations.

4. Compute the calibration factors: calc_calcfactors.py
   Uses the results of 1 and 3 to calculate the calibration factors for all
   observed absflux stars for one filter for all three types (if present).
   
