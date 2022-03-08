import pipeline
from pipeline import Pipeline, Config, Cataloguer, Utilities, PeriodFinder

from datetime import datetime
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Main function for testing things.
    Completely trashable.

    """

    config = Config(
        #raw_image_dir = os.path.expanduser("~/mnt/uni/tmp_moving/0301"),
        #raw_image_dir = os.path.expanduser("~/mnt/jgt/2022/0301"),
        image_prefix = "l138_0",
        n_sets = 9,
        bias_prefix = "bias",
        fits_extension = ".fits",
        fits_date_format = "%Y.%m.%dT%H:%M:%S.%f",
        has_filter_in_header = False,
        n_sample_periods = 100,
        amplitude_score_threshold = 0.85,
    )

    start_time = datetime.now()

    c = Cataloguer(config)
    catalogue_set_number = 1
    catalogue_image_number = 1
    catalogue_image_path = os.path.join(config.image_dir,
            config.image_format_str
            .format(catalogue_set_number, catalogue_image_number))

    n_sources = c.generate_catalogue(catalogue_image_path, solve=False)
    cataloguer_time = Utilities.finished_job("cataloguing stars", start_time)

    Pipeline.run_existing(config, n_sources, cataloguer_time)

def setup():
    """
    Delete certain results from the pipeline.
    Used for debugging
    """

    print("[Setup] Removing results and light curves")
    try:
        shutil.rmtree("./workspace/results/")
        shutil.rmtree("./workspace/rejects/")
        shutil.rmtree("./workspace/adjusted_light_curves/")
        shutil.rmtree("./workspace/periods/")
    except:
        pass

    

def post_processing():
    config = Config(
            raw_image_dir = os.path.expanduser("~/mnt/jgt/2022/0301"),
            #image_prefix = "l138_0",
            #image_prefix = "l140_0",
            image_prefix = "l141_5",
    )

    pf = PeriodFinder(config)

    #search_ids = [1194]
    #search_ids = [760]
    search_ids = [301, 448]

    for _i, path, source_id in Utilities.loop_variables(config, search_ids, adjusted=True):
        curve = np.genfromtxt(path, dtype=config.light_curve_dtype)

        print("[TEST] std before: {}".format(np.std(curve['counts'])))

        P, P_err, A, A_err, phi, offset = pf.period_search_curve(
                source_id, curve['time'], curve['counts'], curve['counts_err'])
        main_sine = A * np.sin(2*np.pi/P * curve['time'] + phi)
        curve['counts'] -= main_sine

        print("[POST] Period:    {:.6f} +/- {:.6f} hours".format(P/3600, P_err/3600))
        print("[POST] Amplitude: {:.6f} +/- {:.6f}".format(A*100, A_err*100))

        P, P_err, A, A_err, phi, offset = pf.period_search_curve(
                source_id, curve['time'], curve['counts'], curve['counts_err'])
        secondary_sine = A * np.sin(2*np.pi/P * curve['time'] + phi)
        curve['counts'] -= secondary_sine
        print("[POST] Period:    {:.6f} +/- {:.6f} hours".format(P/3600, P_err/3600))
        print("[POST] Amplitude: {:.6f} +/- {:.6f}".format(A*100, A_err*100))

        print("[TEST] std after: {}".format(np.std(curve['counts'])))

        plt.scatter(curve['time'], curve['counts'], marker="x")
        #plt.plot(curve['time'], main_sine + secondary_sine + offset, color="orange")
        #plt.plot(curve['time'], main_sine, color="red")
        plt.savefig(os.path.join(config.periods_dir, "TEST4.jpg"))
        plt.close()


if __name__ == "__main__":
    setup()
    main()
    #post_processing()
