import Utilities
import Constants
import Reducer
import Cataloguer
import ShiftFinder
import FluxFinder
import DataAnalyser
import PeriodFinder
import MovingObjectFinder
import StreakFinder

from datetime import datetime
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

def setup():
    print("[Setup] Removing results and light curves")
    try:
        shutil.rmtree("./workspace/results/")
        shutil.rmtree("./workspace/rejects/")
        shutil.rmtree("./workspace/adjusted_light_curves/")
    except:
        pass

def main():

    start_time = datetime.now()
    print("Started at {}".format(start_time.strftime("%H:%M:%S")))
    
    ## Config
    ## Config object for a run of data
    ## See `Constants.py` for available options and default values
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
    
    ## Reducer
    ## Takes raw images, subtracts bias and divides by flat field
    r = Reducer.Reducer(config, "No filter")    ## Only "No filter" for Trius
    #r.reduce(skip_existing=True)    ## Skip images that have already been reduced
    #r.reduce(skip_existing=False)    ## Always create new images

    reducer_time = Utilities.finished_job("reducing images", start_time)
    
    ## Cataloguer
    ## Creates a catalogue of stars found in the given image
    c = Cataloguer.Cataloguer(config)

    catalogue_set_number = 1
    catalogue_image_number = 1

    catalogue_image_path = os.path.join(config.image_dir,
            config.image_format_str
            .format(catalogue_set_number, catalogue_image_number))

    n_sources = c.generate_catalogue(catalogue_image_path, solve=False)
    #c.generate_image_times()

    cataloguer_time = Utilities.finished_job("cataloguing stars", reducer_time)
    

    ## Moving object finder
    #mof = MovingObjectFinder.MovingObjectFinder(Constants.folder)
    #mof.find_moving_objects()
     
     
    ## ShiftFinder
    ## Gets the shift of each star for each image in the series
    sf = ShiftFinder.ShiftFinder(config, n_sources)
    #sf.generate_shifts()
    #reference_ids = sf.get_reference_ids()
     
    shift_finder_time = Utilities.finished_job("finding shifts", cataloguer_time)
     

    ## FluxFinder
    ff = FluxFinder.FluxFinder(config, n_sources)

    ## Find the flux of each star in each image then create a light curve
    ## Write the light curves to file
    #ff.make_light_curves()

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
    #ff.plot_adjusted_comparison(variable_ids, show=False)

    source_ids = da.get_source_ids()
    #ff.plot_all_light_curves(source_ids, adjusted=True)
    #ff.plot_adjusted_comparison(source_ids, show=False)

    reject_ids = np.setxor1d(source_ids, variable_ids)
    ff.plot_all_light_curves(reject_ids, plot_dir=config.rejects_dir, adjusted=True, show=False)

    print("[Main] Outputting results")
    da.output_results()
    da.plot_means_and_stds()
    da.create_thumbnails(ff)
    
    _ = Utilities.finished_job("everything", start_time)


## Rename all the RAW images to a more standard format
## Manually have to enter known fields which do not conform
## TODO: Copy flats for ones that need it?
def rename():
    ## Ideal format:
    ## lxxx_y_s_iii.fit
    ## xxx: llongitude integer part
    ## y: longitude decimal part, 0 or 5
    ## s: set number, 1..n_sets
    ## iii: image_number, 001..050 

    print("[TEST] Starting rename")

    ## l141 first set
    base_path = os.path.expanduser("~/mnt/jgt/2020/0206Trius")
    for i in range(1, 51):
        os.rename(
                os.path.join(base_path, "l141_{:03d}.fit".format(i)),
                os.path.join(base_path, "l141_1_{:03d}.fit".format(i))
        )

    ## l140.5 -> l140_5
    base_path = os.path.expanduser("~/mnt/jgt/2019/0221")
    for s in range(1, 8):
        for i in range(1, 51):
            os.rename(
                    os.path.join(base_path, "l140.5_{:1d}_{:03d}.fit".format(s, i)),
                    os.path.join(base_path, "l140_5_{:1d}_{:03d}.fit".format(s, i))
            )

    ## L140 -> l140_0
    base_path = os.path.expanduser("~/mnt/jgt/2019/0218")
    for s in range(1, 8):
        for i in range(1, 51):
            os.rename(
                    os.path.join(base_path, "L140_{:1d}_{:03d}.fit".format(s, i)),
                    os.path.join(base_path, "l140_0_{:1d}_{:03d}.fit".format(s, i))
            )

    ## Make l138_0 image numbers range from 1-50
    base_path = os.path.expanduser("~/mnt/uni/tmp_moving/0301")
    for i in range(40, 40+50):
        os.rename(
                os.path.join(base_path, "l138_0_1_{:03d}.fits".format(i)),
                os.path.join(base_path, "l138_0_1_{:03d}.fits".format(i-39)),
        )

    for i in range(51, 51+50):
        os.rename(
                os.path.join(base_path, "l138_0_3_{:03d}.fits".format(i)),
                os.path.join(base_path, "l138_0_3_{:03d}.fits".format(i-50)),
        )

    print("[TEST] Finished renaming")


def post_processing():
    config = Constants.Config(
        raw_image_dir = os.path.expanduser("~/mnt/uni/tmp_moving/0301"),
        image_prefix = "l138_0",
        n_sets = 9,
        bias_prefix = "bias",
        fits_extension = ".fits",
        fits_date_format = "%Y.%m.%dT%H:%M:%S.%f",
        has_filter_in_header = False,
    )

    pf = PeriodFinder.PeriodFinder(config)

    search_ids = [1194]

    for _i, path, source_id in Utilities.loop_variables(config, search_ids, adjusted=True):
        curve = np.genfromtxt(path, dtype=config.light_curve_dtype)
        P, P_err, A, A_err, phi, B = pf.period_search(source_id, path)
        main_sine = A * np.sin(2*np.pi/P * curve['time'] + phi) + B
        P, P_err, A, A_err, phi, B = pf.period_search(source_id, path, period_min = 9000, period_max = 10000)
        secondary_sine = A * np.sin(2*np.pi/P * curve['time'] + phi)

        plt.scatter(curve['time'], curve['counts'])
        plt.plot(curve['time'], main_sine + secondary_sine, color="red")
        plt.savefig(os.path.join(config.periods_dir, "TEST4.jpg"))
        plt.close()


if __name__ == "__main__":
    #rename()
    #setup()
    #main()
    post_processing()
