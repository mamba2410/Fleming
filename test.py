import pipeline
from pipeline import Pipeline, Config, Cataloguer, Utilities, PeriodFinder
from pipeline import VariableDetector, DataAnalyser

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

    wa = 1.5
    vt = 1.0
    at = 2.0

    config = Config(
        image_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images/"),
        fits_extension = ".fits",
        image_prefix = "l138_0",
        n_sets = 9,
        variability_threshold = vt,
        amplitude_score_threshold = at,
        period_width_adjustment = wa,
    )

    #Pipeline.run_analysis(config, assume_already_adjusted=True)

    config = Config(
        image_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images/"),
        image_prefix = "l140_0",
        n_sets = 7,
        variability_threshold = vt,
        amplitude_score_threshold = at,
        period_width_adjustment = wa,
    )


    #Pipeline.run_analysis(config, assume_already_adjusted=True)

    config = Config(
        image_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images/"),
        image_prefix = "l141_5",
        n_sets = 8,
        variability_threshold = vt,
        amplitude_score_threshold = at,
        period_width_adjustment = wa,
    )

    Pipeline.run_analysis(config, assume_already_adjusted=True)


def main_amp_testing():
    config = Config(
        image_dir = os.path.expanduser("~/mnt/data/jgt_images/"),
        image_prefix = "l141_5",
        n_sets = 8,
        amplitude_score_threshold = 1.0,
    )

    source_ids = [301, 448]
    n_sources = len(source_ids)

    #config = Config(
    #    image_dir = os.path.expanduser("~/mnt/data/jgt_images/"),
    #    image_prefix = "l140_0",
    #    n_sets = 8,
    #    amplitude_score_threshold = 1.0,
    #)

    #source_ids = [760]
    #n_sources = len(source_ids)

    #config = Config(
    #    image_dir = os.path.expanduser("~/mnt/data/jgt_images/"),
    #    image_prefix = "l138_0",
    #    n_sets = 9,
    #    amplitude_score_threshold = 1.0,
    #)

    source_ids = [47]
    n_sources = len(source_ids)


    da = DataAnalyser(config, adjusted=True)
    means, stds, meds, n_pos = da.get_means_and_stds()

    vd = VariableDetector(config, source_ids, means, stds, meds, n_pos, adjusted=True)
    _ = vd.std_dev_search(1e3, 1e3, 1e3, 1e3)
    variable_ids = vd.amplitude_search(-1.0)

    da.output_results(variable_ids, vd, out_dir=config.testing_dir)

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
        shutil.rmtree("./workspace/testing/")
    except:
        pass

    

def post_processing():
    config = Config(
            image_prefix = "l138_0",
            #image_prefix = "l140_0",
            #image_prefix = "l141_5",
    )

    pf = PeriodFinder(config)

    search_ids = [1194]
    #search_ids = [760]
    #search_ids = [301, 448]

    for _i, path, source_id in Utilities.loop_variables(config, search_ids, adjusted=True):
        curve = np.genfromtxt(path, dtype=config.light_curve_dtype)
        counts_change = np.copy(curve['counts'])

        print("[TEST] std before: {}".format(np.std(curve['counts'])))

        P, P_err, A, A_err, phi, offset = pf.period_search_curve(
                source_id, curve['time'], counts_change, curve['counts_err'])

        main_sine = A * np.sin(2*np.pi/P * curve['time'] + phi)
        counts_change -= main_sine

        #plt.scatter(curve['time'], curve['counts'], marker="x")
        #plt.plot(curve['time'], main_sine + offset, color="red")
        plt.scatter(curve['time'], counts_change, marker="x")
        plt.savefig(os.path.join(config.periods_dir, "TEST2.jpg"))
        plt.close()

        print("[POST] Period:    {:.6f} +/- {:.6f} hours".format(P/3600, P_err/3600))
        print("[POST] Amplitude: {:.6f} +/- {:.6f}".format(A*100, A_err*100))

        P, P_err, A, A_err, phi, offset = pf.period_search_curve(
                source_id, curve['time'], counts_change, curve['counts_err'])

        secondary_sine = A * np.sin(2*np.pi/P * curve['time'] + phi)
        counts_change -= secondary_sine

        print("[POST] Period:    {:.6f} +/- {:.6f} hours".format(P/3600, P_err/3600))
        print("[POST] Amplitude: {:.6f} +/- {:.6f}".format(A*100, A_err*100))

        print("[TEST] std after: {}".format(np.std(curve['counts'])))

        plt.scatter(curve['time'], curve['counts'], marker="x")
        plt.plot(curve['time'], secondary_sine + offset, color="red")
        plt.savefig(os.path.join(config.periods_dir, "TEST1.jpg"))
        plt.close()


if __name__ == "__main__":
    #setup()
    #main()
    #main_amp_testing()
    post_processing()
