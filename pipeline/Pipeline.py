from pipeline import Config, Reducer, Cataloguer, ShiftFinder, FluxFinder
from pipeline import DataAnalyser, Utilities, VariableDetector

from datetime import datetime
import os
import numpy as np

HOME = os.path.expanduser("~")

## Run for a single field defined by a config object
def run(config, show_plots=False, show_errors=False, solve_astrometry=True, skip_existing_images=True):

    start_time = datetime.now()
    print("[JOB] Started {} at {}".format(config.image_prefix, start_time.strftime("%H:%M:%S")))
    
    n_sources, time = images_to_light_curves(config, start_time,
            skip_existing_images, solve_astrometry=solve_astrometry)

    _results_table, output_time = run_existing(config, n_sources, time,
            show_plots=show_plots, show_errors=show_errors)

    
    _ = Utilities.finished_job("everything for {}".format(config.image_prefix), start_time)


## Run only the variable detection for a single field
## Assumes we already have shifts, times, light curves, adjusted light curves
def run_analysis(config, show_plots=False, show_errors=False, assume_already_adjusted=True):
    start_time = datetime.now()

    c = Cataloguer(config)

    catalogue_set_number = 1
    catalogue_image_number = 1

    catalogue_image_path = os.path.join(config.image_dir,
            config.image_format_str
            .format(catalogue_set_number, catalogue_image_number))

    print("[Pipeline] Cataloguing image {}".format(catalogue_image_path))
    n_sources = c.generate_catalogue(catalogue_image_path, solve=False, write_file=False)

    cataloguer_time = Utilities.finished_job("cataloguing stars", start_time)

    run_existing(config, n_sources, cataloguer_time,
            show_plots=show_plots, show_errors=show_errors,
            assume_already_adjusted=assume_already_adjusted)

    _ = Utilities.finished_job("everything for {}".format(config.image_prefix), start_time)


def images_to_light_curves(config, start_time, skip_existing_images=True, solve_astrometry=True):

    ## Reducer
    ## Takes raw images, subtracts bias and divides by flat field

    r = Reducer(config, "No filter")        ## Only "No filter" for Trius
    r.reduce(skip_existing=skip_existing_images)    ## Skip images that have already been reduced

    reducer_time = Utilities.finished_job("reducing images", start_time)
    
    ## Cataloguer
    ## Creates a catalogue of stars found in the given image
    c = Cataloguer(config)

    catalogue_set_number = 1
    catalogue_image_number = 1

    catalogue_image_path = os.path.join(config.image_dir,
            config.image_format_str
            .format(catalogue_set_number, catalogue_image_number))

    print("[Pipeline] Cataloguing image {}".format(catalogue_image_path))
    n_sources = c.generate_catalogue(catalogue_image_path, solve=solve_astrometry)
    c.generate_image_times()

    cataloguer_time = Utilities.finished_job("cataloguing stars", reducer_time)
    

    ## Moving object finder
    #mof = MovingObjectFinder.MovingObjectFinder(Constants.folder)
    #mof.find_moving_objects()
     
     
    ## ShiftFinder
    ## Gets the shift of each star for each image in the series
    print("[Pipeline] Finding shifts in each image")
    sf = ShiftFinder(config, n_sources)
    sf.generate_shifts()
    reference_ids = sf.get_reference_ids()
     
    shift_finder_time = Utilities.finished_job("finding shifts", cataloguer_time)
     

    ## FluxFinder
    ff = FluxFinder(config, n_sources)
    print("[Pipeline] Making initial light curves")

    ## Find the flux of each star in each image then create a light curve
    ## Write the light curves to file
    ff.make_light_curves()

    light_curve_time = Utilities.finished_job("making light curves", shift_finder_time)

    return n_sources, light_curve_time


def run_existing(config, n_sources, start_time,
        show_plots=False, show_errors=False, assume_already_adjusted=False):
    """
    Run for a field and assume the light curves are already present.
    Very useful for fast debugging the later stages of the pipeline like the variable finder etc.

    Parameters
    ----------

    assume_already_adjusted: bool
        Assume we already have adjusted light curves.
        Set to False if you're debugging the adjustment process.

    """

    ## FluxFinder
    ff = FluxFinder(config, n_sources)
    
    ## Assume already have unadjusted light curves

    ## Sometimes we want to remake adjusted light curves
    if not assume_already_adjusted:

        print("[Pipeline] Creating average light curve")
        da = DataAnalyser(config, adjusted=False)
        mean, std, med, n_positive = da.get_means_and_stds()
        source_ids = da.get_source_ids()
        da.plot_means_and_stds()
        
        vd = VariableDetector(config, source_ids, mean, std, med, n_positive, adjusted=False)
        exclude_ids = vd.std_dev_search(config.avg_exclude_threshold, 1e4, 0, 0)
        avg_ids = da.get_ids_for_avg(exclude_ids)

        da.make_avg_curve(avg_ids)
        make_avg_curve_time = Utilities.finished_job("making average curve", start_time)

        ## 'adjusts' light curves by dividing by average
        print("[Pipeline] Adjusting")
        ff.create_adjusted_light_curves(source_ids, std)

        adjustment_time = Utilities.finished_job("adjusting light curves", make_avg_curve_time)
    else:
        adjustment_time = start_time



    ## Now we have adjusted light curves
    ## New DataAnalyser for non-adjusted
    da = DataAnalyser(config, adjusted=True)
    means_adj, stds_adj, medians_adj, n_positive_adj = da.get_means_and_stds()
    source_ids = da.get_source_ids()
    da.plot_means_and_stds()

    print("[Pipeline] Getting variables post-adjustment")
    vd = VariableDetector(config, source_ids, means_adj, stds_adj,
            medians_adj, n_positive_adj, adjusted=True)

    n_measures = config.n_sets * config.set_size
    variable_ids_s = vd.std_dev_search(config.variability_threshold,
            config.variability_max, config.min_signal_to_noise, n_measures)
    variable_ids_a = vd.amplitude_search(config.amplitude_score_threshold)
    variable_ids_a = vd.filter_variables(variable_ids_a)

    variable_ids = np.unique(np.concatenate((variable_ids_a, variable_ids_s)))
    #variable_ids = variable_ids_a

    post_time = Utilities.finished_job("post-adjustment", adjustment_time)

    print("[Pipeline] Plotting variable curves")
    ff.plot_avg_light_curve(config.avg_curve_path, show=show_plots, show_errors=show_errors)
    ff.plot_given_light_curves(variable_ids, adjusted=True, show=show_plots, show_errors=show_errors)
    #ff.plot_adjusted_comparison(variable_ids, plot_dir=config.testing_dir,
    #        show=show_plots, show_errors=show_errors)

    print("[Pipeline] Outputting results")
    results_table = da.output_results(variable_ids, vd)
    #results_table = da.output_results(source_ids, vd)
    ff.create_thumbnails(results_table)
    print("[Pipeline] Variables found: total {}; std {}; amp {}"
            .format(len(variable_ids), len(variable_ids_s), len(variable_ids_a)))

    output_time = Utilities.finished_job("outputting results", post_time)
    
    return results_table, output_time


