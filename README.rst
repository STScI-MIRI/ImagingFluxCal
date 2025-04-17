JWST/MIRI Imaging/Coronagraphy Flux Calibration
===============================================

Code to measure flux calibration stars and compute the
flux calibration for JWST/MIRI Imaging and Coronagraphic filters.

Main goal is to fully automate this work to allow for straightforward
recalculations from the raw data to the final calibration factors for each
filter.  This enables tests of the sensitivity of the flux calibration
to changes in the jwst pipeline or any of the steps in the flux calibration
factor calculation.

Citation
--------

The results of this work are published in
`Gordon et al. 2025, AJ, 169, 6 <https://ui.adsabs.harvard.edu/abs/2025AJ....169....6G>`_.
Please cite this paper if this repository is used.

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

1. Reduce data with jwst pipeline: `pipeline_one_filter.py`.
   Runs all three stages of pipeline for a one type of star and one filter.
   There are some differences from the defaults for MIRI.
   `pipeline_all` runs everything for all MIRI imaging/coronagraphic filters.

2. Compute WebbPSFs: `calc_webbpsfs.py`.
   Runs WebbPSFs with specfic parameters for MIRI.  Needs to be done for each filter.
   Critial for coronagraphy to set the normalize="exit_pupil" to ensure the 
   WebbPSF is normalized to a total of 1.

3. Compute encircled energies and individual observation aperture corrections.
   Uses `calc_encircled_energy.py` to compute the encircled energies
   Uses observations of a bright stars and PSFs from WebbPSFs to compute aperture
   corrections.  Uses the bright star observations for the inner region and
   the WebbPSFs for the outer regions.
   `calc_all_ees` runs all MIRI imaging filters for a defined set of observations.
   This script varies the fwhmfac where the observed/webbpsf PSFs are merged.
   `calc_all_ees_onlyneeded` only runs for the fwhmfacs that are used in the apcor
   ref file and paper plots.

4. Create the aperture correction reference file: `create_apcor_reffile.py`
   Uses the individual aperture corrections computed in step 3 to compute the
   aperture corrections with uncertainties (where possible) to create the
   apcor reference file in the correct format.
   Manually specified set of files that have been visual vetted from the results
   of step 3.

5. Measure the flux in a fixed aperture: aper_one_filter.py.
   Measures the flux using an aperture and background annulus.
   `aper_all.py` does this for all MIRI imaging and coronagraphic filters.

6. Compute the sensitivity loss using the repeatability observations of 
   BD+60 1753.  First run `plot_repeatability.py` to fit the observations with
   an exponential model for all the Imager filters.  This produces data files
   with the fit coeffiecients.  Then run `generate_coron_repeat_from_imager.py`
   to create similar fit coefficient data files for the Coronagraphic filters 
   based on interpolating the Imager results based on filter pivot wavelength.

7. Compute the model flux densities in all the filters: `model_fluxes.py --grieke`

8. Measure the subarray dependence: `plot_subtrans_allobs.py`.
   Adjust the "Adopted" points manually till they match the scatter in the ratio
   of the calibration factors.  Copy the values relative to FULL printed to 
   the screen to the `subarr_cor` variable in the `plot_calfacs` function.

   Before this plot is made, the subarray dependence based on the dedicated subarray
   observations need to be measured by running `plot_subtrans.py --filter=F770W`
   and the same again for `F1280W`.  In addition, the subarray dependence based
   on the ratio of calibration factors from all the observations needs to measured
   using `plot_multi_calfacs.py --grieke --xaxisval=subarr --noignore`.

9. Compute the calibration factors: `calc_calfactors.py`
   Uses the results of 4, 5, and 6 to calculate the calibration factors for all
   observed absflux stars for one filter for all three types (if present).
   Produces a table giving the calibration factors for each observation.
   Produces plots of calibration factors versus model flux, time, well depth,
   etc.
   Example: `calc_calfactors.py --grieke --applytime --subarrcor --nocurval`

   A meta plot for all the filters can be created using `plot_multi_calfacs.py`.

10. Compute a set of photom reference files accounting for the time dependent
   sensitivity losses using `create_photom_reffile.py`

Figures
-------

1. Montage of observed imager PSFs: `Plotting/plot_example_images.py`

2. Montage of observed coronagraphic PSFs: `Plotting/plot_example_images.py --coron`

3. Encirciled energy plot: `Plotting/plot_encircled_energy.py`
   By default makes a figure showing the encircled energies for all 9
   imaging filters, with a veritical offset between them.  Call with
   `--coron` to get the coronagraphic plot.

4. Example of photometry technique: `aper_one_filter.py`
   One of the standard outputs.

5. Time dependent sensitivity loss: `plot_repeatability.py`

6. Subarray dependence: `plot_subtrans_allobs.py`.
   `plot_subtrans.py --filter=F770W` and `plot_subtrans.py --filter=F1280W` need
   to be run first to get the values for the dedicated subarray transfer observations.

7. Source dependencies:
   `calc_calfactors.py --filter=F1280W --sourcemulti --grieke --subarrcor --applytime --nocurval`

8. Source type dependence: `plot_srctype_allobs.py`

9. Detector dependencies:
   `calc_calfactors.py --filter=F1280W --detmulti --grieke --subarrcor --applytime --nocurval`

10. Change in delivered calibration versus time: `Plotting/plot_del_photom_vs_time.py`

Appendix Figures:

`plot_multi_calfacs.py --xaxisval=mflux --grieke --subarrcor`

Replace `mflux` with desired xaxis value.

Tables
------

1, 2, 3. Observation details: `Tables/create_obstable.py`
   Output to screen.

4. Aperture corrections: `create_apcor_reffile.py`
   Portion of output to the screen.

5. Output from creating Subarray dependence.

6. Output from `create_photom_reffile.py`.


Since paper
-----------

Work continues.

1. Plotting relative photometry versus dither position.  `Plotting/plot_dither_pos.py`
   Results published as a JWST report.
