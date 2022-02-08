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

    start_time = datetime.now()
    print("Started at {}".format(start_time.strftime("%H:%M:%S")))
    
    ## Config
    ## Config object for a run of data
    ## See `Constants.py` for available options and default values
    config = Constants.Config(
        image_dir = "/home/callum/mnt/data/jgtdata/l137_0",
        image_prefix = "l137_0",
        has_sets = True,
        set_size = 50,
        n_sets = 7
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

    c.catalogue(catalogue_image_path)
    n_sources = c.get_n_sources()

    cataloguer_time = Utilities.finished_job("cataloguing stars", reducer_time)
    

    ## Moving object finder
    #mof = MovingObjectFinder.MovingObjectFinder(Constants.folder)
    #mof.find_moving_objects()
     
     
    ## ShiftFinder
    ## Gets the shift of each star for each image in the series
    sf = ShiftFinder.ShiftFinder(config)
    sf.get_all_shifts()
     
    shift_finder_time = Utilities.finished_job("finding shifts", cataloguer_time)
     

    ## FluxFinder
    ff = FluxFinder.FluxFinder(config, n_sources)

    ## Find the flux of each star in each image then create a light curve
    ## Write the light curves to file
    ff.find_all_fluxes()
    ff.make_light_curves()

    light_curve_time = Utilities.finished_job("making light curves", shift_finder_time)
    
    ## DataAnalyser
    da = DataAnalyser.DataAnalyser(config)
    
    print("[Main] Getting variables")
    da.get_means_and_stds(adjusted=False)
    variable_ids = da.get_variables(ff, adjusted=False)
    da.plot_means_and_stds(adjusted=False)

    find_variables_time =Utilities.finished_job("finding un-adjusted variables", light_curve_time)

    ## Create average light curve of all viable stars
    print("[Main] Making average curve")
    avg_ids = da.get_ids_for_avg()
    ff.make_avg_curve(avg_ids)

    ## 'adjusts' light curves by dividing by average
    print("[Main] Adjusting")
    ff.divide_by_average()

    adjustment_time = Utilities.finished_job("adjusting light curves", find_variables_time)

    avg_fname = "{}_avg{}".format(config.image_prefix, config.standard_file_extension)
    avg_path = os.path.join(config.workspace_dir, avg_fname)

    print("[Main] Plotting variable curves")
    ff.plot_avg_light_curve(avg_path, adjusted=True, show=False)
    ff.plot_all_light_curves(variable_ids, adjusted=True, show=False)

    print("[Main] Getting variables post-adjustment")
    da.get_means_and_stds(adjusted=True)
    _ = da.get_variables(ff, adjusted=True)
    da.plot_means_and_stds(adjusted=True)

    _ = Utilities.finished_job("post-adjustment", adjustment_time)

    print("[Main] Outputting results")
    da.output_results()
    da.create_thumbnails(ff)
    
    _ = Utilities.finished_job("everything", start_time)


if __name__ == "__main__":
    main()
