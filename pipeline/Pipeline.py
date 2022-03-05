from pipeline import Config, Reducer, Cataloguer, ShiftFinder, FluxFinder, DataAnalyser, Utilities

from datetime import datetime
import os

HOME = os.path.expanduser("~")

## Run for a single field defined by a config object
def run(config, show_plots=False, show_errors=False, solve_astrometry=True, skip_existing_images=True):

    start_time = datetime.now()
    print("[JOB] Started {} at {}".format(config.image_prefix, start_time.strftime("%H:%M:%S")))
    
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

    
    ## DataAnalyser
    da = DataAnalyser(config)
    
    print("[Pipeline] Creating average light curve")
    da.create_avg_curve()

    make_avg_curve_time = Utilities.finished_job("making average curve", light_curve_time)

    ## 'adjusts' light curves by dividing by average
    print("[Pipeline] Adjusting")
    ff.create_adjusted_light_curves()

    adjustment_time = Utilities.finished_job("adjusting light curves", make_avg_curve_time)

    print("[Pipeline] Getting variables post-adjustment")
    variable_ids = da.get_variables()

    _ = Utilities.finished_job("post-adjustment", adjustment_time)

    print("[Pipeline] Plotting variable curves")
    ff.plot_avg_light_curve(config.avg_curve_path, show=show_plots, show_errors=show_errors)
    ff.plot_given_light_curves(variable_ids, adjusted=True, show=show_plots, show_errors=show_errors)

    print("[Pipeline] Outputting results")
    da.output_results()
    da.create_thumbnails(ff)
    
    _ = Utilities.finished_job("everything for {}".format(config.image_prefix), start_time)


def run_existing(config, show_plots=False, show_errors=False, solve_astrometry=False,
        assume_already_adjusted=True):
    """
    Run for a field and assume the light curves are already present.
    Very useful for fast debugging the later stages of the pipeline like the variable finder etc.

    Parameters
    ----------

    assume_already_adjusted: bool
        Assume we already have adjusted light curves.
        Set to False if you're debugging the adjustment process.

    """

    start_time = datetime.now()
    print("[JOB] Started {} at {}".format(config.image_prefix, start_time.strftime("%H:%M:%S")))
    
    ## Reducer
    ## Assume images already reduced
    
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
    ## Assume already have image time file
    #c.generate_image_times()

    cataloguer_time = Utilities.finished_job("cataloguing stars", start_time)
     
     
    ## ShiftFinder
    ## Assume already have shifts file

    ## FluxFinder
    ff = FluxFinder(config, n_sources)
    
    ## Assume already have unadjusted light curves
    
    ## DataAnalyser
    da = DataAnalyser(config)

    ## Sometimes we want to remake adjusted light curves
    if not assume_already_adjusted:
        print("[Pipeline] Creating average light curve")
        da.create_avg_curve()

        make_avg_curve_time = Utilities.finished_job("making average curve", cataloguer_time)

        ## 'adjusts' light curves by dividing by average
        print("[Pipeline] Adjusting")
        ff.create_adjusted_light_curves()

        adjustment_time = Utilities.finished_job("adjusting light curves", make_avg_curve_time)
    else:
        adjustment_time = cataloguer_time

    print("[Pipeline] Getting variables post-adjustment")
    variable_ids = da.get_variables()

    _ = Utilities.finished_job("post-adjustment", adjustment_time)

    print("[Pipeline] Plotting variable curves")
    ff.plot_avg_light_curve(config.avg_curve_path, show=show_plots, show_errors=show_errors)
    ff.plot_given_light_curves(variable_ids, adjusted=True, show=show_plots, show_errors=show_errors)

    print("[Pipeline] Outputting results")
    da.output_results()
    da.create_thumbnails(ff)
    
    _ = Utilities.finished_job("everything for {}".format(config.image_prefix), start_time)
