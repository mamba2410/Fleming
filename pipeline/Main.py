import Utilities
import Constants
import Reducer
import Cataloguer
import ShiftFinder
import FluxFinder
import DataAnalyser
import MovingObjectFinder
import StreakFinder

from datetime import datetime
import os

def main():

    HOME = os.path.expanduser("~")

    start_time = datetime.now()
    print("Started at {}".format(start_time.strftime("%H:%M:%S")))
    
    ## Config
    ## Config object for a run of data
    ## See `Constants.py` for available options and default values
    config = Constants.Config(
        image_dir = os.path.expanduser("~/mnt/data/jgtdata/l137_0/0121"),
        image_prefix = "l137_0",
        has_sets = True,
        n_sets = 7,
        set_size = 50,
    )
    
    ## Reducer
    ## Takes raw images, subtracts bias and divides by flat field
    r = Reducer.Reducer(config, "No filter")    ## Only "No filter" for Trius
    r.reduce(skip_existing=True)    ## Skip images that have already been reduced

    reducer_time = Utilities.finished_job("reducing images", start_time)
    
    
    ## Cataloguer
    ## Creates a catalogue of stars found in the given image
    c = Cataloguer.Cataloguer(config)

    catalogue_set_number = 1
    catalogue_image_number = 1

    catalogue_image_path = os.path.join(config.image_dir,
            config.image_format_str
            .format(catalogue_set_number, catalogue_image_number))

    n_sources = c.generate_catalogue(catalogue_image_path, solve=True)
    c.generate_image_times()

    cataloguer_time = Utilities.finished_job("cataloguing stars", reducer_time)
    

    ## Moving object finder
    #mof = MovingObjectFinder.MovingObjectFinder(Constants.folder)
    #mof.find_moving_objects()
     
     
    ## ShiftFinder
    ## Gets the shift of each star for each image in the series
    sf = ShiftFinder.ShiftFinder(config, n_sources)
    sf.generate_shifts()
    reference_ids = sf.get_reference_ids()
     
    shift_finder_time = Utilities.finished_job("finding shifts", cataloguer_time)
     

    ## FluxFinder
    ff = FluxFinder.FluxFinder(config, n_sources)

    ## Find the flux of each star in each image then create a light curve
    ## Write the light curves to file
    ff.make_light_curves()

    light_curve_time = Utilities.finished_job("making light curves", shift_finder_time)

    
    ## DataAnalyser
    da = DataAnalyser.DataAnalyser(config)
    
    print("[Main] Creating average light curve")
    da.create_avg_curve()

    make_avg_curve_time = Utilities.finished_job("making average curve", light_curve_time)

    ## 'adjusts' light curves by dividing by average
    print("[Main] Adjusting")
    ff.create_adjusted_light_curves()

    adjustment_time = Utilities.finished_job("adjusting light curves", make_avg_curve_time)

    print("[Main] Getting variables post-adjustment")
    variable_ids = da.get_variables()

    _ = Utilities.finished_job("post-adjustment", adjustment_time)

    print("[Main] Plotting variable curves")
    ff.plot_avg_light_curve(config.avg_curve_path, show=False)
    ff.plot_all_light_curves(variable_ids, adjusted=True, show=False)

    print("[Main] Outputting results")
    da.output_results()
    da.create_thumbnails(ff)
    
    _ = Utilities.finished_job("everything", start_time)


if __name__ == "__main__":
    main()
