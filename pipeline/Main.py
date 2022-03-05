import pipeline.Utilities
import pipeline.Constants
import pipeline.Reducer
import pipeline.Cataloguer
import pipeline.ShiftFinder
import pipeline.FluxFinder
import pipeline.DataAnalyser
import pipeline.MovingObjectFinder
import pipeline.StreakFinder

from datetime import datetime
import os

HOME = os.path.expanduser("~")

def main():

    start_time = datetime.now()
    print("[MAIN] Started everything at {}".format(start_time.strftime("%H:%M:%S")))

    field_details = [
            # image_dir             image_prefix    n_sets  flat_prefix bias_prefix
            #["~/mnt/jgt/2019/1028", "l135",         6,      "noflat",   "bias"],
            #["~/mnt/jgt/2021/1210Trius", "l135_5",  4,      "dflat7",   "bias2"],
            #["~/mnt/jgt/2022/0105", "l136",         9,      "dflat",    "bias2"],
            #["~/mnt/jgt/2022/0117", "l136_5",       7,      "dflat",    "bias4"],
            #["~/mnt/jgt/2022/0121", "l137_0",       7,      "dflat",    "bias_shutter"],
            #["~/mnt/jgt/2022/0124", "l137_5",       7,      "dflat",    "bias"],
            #["~/mnt/jgt/2022/0301", "l138_0",       9,      "noflat",   "bias_end"],
            #["~/mnt/windows/Users/callu/tmp/0301", "l138_0", 9, "noflat", "bias_end", ".fits"],
            ["~/mnt/jgt/2019/0218", "l140_0",       7,      "noflat",   "BiasLast"],
            #["~/mnt/jgt/2019/0221", "l140_5",       7,      "noflat",   "Bias_end_"],
            #["~/mnt/jgt/2020/0206Trius", "l141",    7,      "flat",     "bias2"],
            ["~/mnt/jgt/2020/0212", "l141_5",       8,      "flat",     "bias2"],
            ]

    for raw_image_dir, image_prefix, n_sets, flat_prefix, bias_prefix in field_details:

        config = Constants.Config(
            raw_image_dir = os.path.expanduser(raw_image_dir),
            image_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images"),
            image_prefix = image_prefix,
            n_sets = n_sets,
            flat_prefix = flat_prefix,
            bias_prefix = bias_prefix,
        )

        run(config)

    _ = Utilities.finished_job("everything", start_time)


def main_single():

    config = Constants.Config(
        raw_image_dir = os.path.expanduser("~/mnt/uni/tmp_moving/0301"),
        image_prefix = "l138_0",
        n_sets = 9,
        bias_prefix = "bias",
        fits_extension = ".fits",
        fits_date_format = "%Y.%m.%dT%H:%M:%S.%f",
        has_filter_in_header = False,
        variability_threshold = 0.5,
    )

    run(config)


def run(config):

    start_time = datetime.now()
    print("[JOB] Started {} at {}".format(config.image_prefix, start_time.strftime("%H:%M:%S")))
    
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
    
    _ = Utilities.finished_job("everything for {}".format(config.image_prefix), start_time)


if __name__ == "__main__":
    main()
    #main_single()
