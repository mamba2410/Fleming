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

    c = Cataloguer(config)

    catalogue_set_number = 1
    catalogue_image_number = 1

    catalogue_image_path = os.path.join(config.image_dir,
            config.image_format_str
            .format(catalogue_set_number, catalogue_image_number))

    n_sources = c.generate_catalogue(catalogue_image_path, solve=True, write_file=True)



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


if __name__ == "__main__":
    #setup()
    main()
    #main_amp_testing()

