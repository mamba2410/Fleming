import pipeline
from pipeline import Pipeline, Config, Cataloguer, Utilities, PeriodFinder
from pipeline import VariableDetector

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

    vt = 1.0
    at = 2.0

    config = Config(
        image_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images/"),
        fits_extension = ".fits",
        image_prefix = "l138_0",
        n_sets = 9,
        variability_threshold = vt,
        amplitude_score_threshold = at,
    )

    Pipeline.run_analysis(config, assume_already_adjusted=True)

    config = Config(
        image_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images/"),
        image_prefix = "l140",
        n_sets = 7,
        variability_threshold = vt,
        amplitude_score_threshold = at,
    )


    Pipeline.run_analysis(config, assume_already_adjusted=True)

    config = Config(
        image_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images/"),
        image_prefix = "l141_5",
        n_sets = 8,
        variability_threshold = vt0,
        amplitude_score_threshold = at,
    )

    Pipeline.run_analysis(config, assume_already_adjusted=True)


def main_amp_testing():
    config = Config(
        raw_image_dir = os.path.expanduser("~/mnt/jgt/2020/0212"),
        image_dir = os.path.expanduser("~/mnt/data/jgt_images/"),
        image_prefix = "l145_5",
        n_sets = 8,
    )

    source_ids = [301]
    n_sources = len(source_ids)
    means = np.ones(n_sources)
    stds = np.ones(n_sources)
    medians = np.ones(n_sources)
    n_positive = np.array([config.n_sets*config.set_size])

    vd = VariableDetector(config, source_ids, means, stds, medians, n_positive, adjusted=True)
    variable_ids = vd.amplitude_search(0.0)

    print(variable_ids)



def setup():
    """
    Delete certain results from the pipeline.
    Used for debugging
    """

    print("[Setup] Removing results and light curves")
    try:
        shutil.rmtree("./workspace/results/")
        #shutil.rmtree("./workspace/rejects/")
        #shutil.rmtree("./workspace/adjusted_light_curves/")
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
    #main_amp_testing()
    #post_processing()
