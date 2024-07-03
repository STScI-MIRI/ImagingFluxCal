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

1. Reduce data with jwst pipeline: `pipeline_one_filter.py`.
   Runs all three stages of pipeline for a one type of star and one filter.
   There are some differences from the defaults for MIRI.
   `pipeline_all` runs everything for all MIRI imaging/coronagraphic filters.

2a. Compute WebbPSFs: `calc_webbpsfs.py`.
   Runs WebbPSFs with specfic parameters for MIRI.  Needs to be done for each filter.
   Critial for coronagraphy to set the normalize="exit_pupil" to ensure the 
   WebbPSF is normalized to a total of 1.

2b. Add the cruciform to the F560W and F770W WebbPSFs: `add_cruciform_to_webbpsfs`.
   Adds the cruciform artifact to the WebbPSFs based on a simple model.  Done for 
   each filter.

3. Compute encircled energies and individual observation aperture corrections.
   Uses `calc_encircled_energy.py` to compute the encircled energies
   Uses observations of a bright stars and PSFs from WebbPSFs to compute aperture
   corrections.  Uses the bright star observations for the inner region and
   the WebbPSFs for the outer regions.
   `calc_all_ees` runs all MIRI imaging filters for a defined set of observations.
   This script varies the fwhmfac where the observed/webbpsf PSFs are merged.

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

8. Compute the calibration factors: `calc_calfactors.py`
   Uses the results of 4, 5, and 6 to calculate the calibration factors for all
   observed absflux stars for one filter for all three types (if present).
   Produces a table giving the calibration factors for each observation.
   Produces plots of calibration factors versus model flux, time, well depth,
   etc.
   A meta plot for all the filters can be created using `plot_multi_calfacs.py`.
   Example: `calc_calfactors.py --grieke --applytime --nocurval --noignore`

9. Compute a set of photom reference files accounting for the time dependent
   sensitivity losses using `python create_photom_reffile.py`

Figures
-------

1. Montage of observed PSFs: `Plotting/plot_example_images.py`
   By default makes a figure for the 9 imaging filters.  Use the option
   `--coron` to make the figure with the 4 coronagraphic filters.

2. Encirciled energy plot: `Plotting/plot_encircled_energy.py`
   By default makes a figure showing the encircled energies for all 9
   imaging filters, with a veritical offset between them.  Call with
   `--coron` to get teh coronagraphic plot.

3. Example of photometry technique: `aper_one_filter.py`
   One of the standard outputs.

4. Example of model spectra and photometry: `plot_example_model_fluxes.py`

5. Time dependent sensitivity loss: `plot_repeatability.py`

6. Subarray dependence: `plot_subtrans_allobs.py`.
   `plot_subtrans.py --filter=770W` and `plot_subtrans.py --filter=F1280W` need
   to be run first to get the values for the dedicated subarray transfer observations.

7. Source dependencies:
   `calc_calfactors.py -filter=F1280W --sourcemulti --grieke --subarrcor --applytime --noignore --nocurval`

8. Source type dependence: `plot_srctype_allobs.py`

9. Detector dependencies:
   `calc_calfactors.py -filter=F1280W --detmulti --grieke --subarrcor --applytime --noignore --nocurval`

Appendix Figures:

`plot_multi_calfacs --xaxisval=mflux --griek --subarrcor`

Replace `mflux` with desired xaxis value.

Tables
------

1, 2, 3. Observation details: `Tables/create_obstable.py`
   Output to screen.

4. Aperture corrections: `create_apcor_reffile.py`
   Portion of output to the screen.

5. Output from creating Subarray dependence.

6. Output from `create_photom_reffile.py`.