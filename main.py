from pipeline import Utilities, Pipeline, Config, PeriodFinder

from datetime import datetime
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord 

def main():
    """
    Main entrypoint for the pipeline process.

    Here the end-user can write config objects for each specific field and 
    then run the pipeline.

    Use one `Config` object per field.

    Things you WILL need to change per field config:

    - `Config.raw_image_dir` - set this to the directory where the raw image files are
    located for a the field.
    - `Config.image_dir` - Set this to the directory where you want the reduced image
    files to be stored. You will need a lot of space (~3GB per field).
    - `Config.image_prefix` - name of the field, which is also the prefix for each image.


    Provided is a loop or two which allows for easy running of all fields in consecutive
    order.

    The `field_details` variable is an array containing fundamental details of 
    each field and things that are different across each field.
    Things like the image prefix, the number of sets the field has, 
    what the flat or bias prefix is and the base location of the raw images.

    The only things you'll need to change (to run on your computer) are where the raw
    images are stored on your drive.

    Pro tip: `~/` is a shortcut for your home folder, even on Windows.

    """

    start_time = datetime.now()
    print("[MAIN] Started everything at {}".format(start_time.strftime("%H:%M:%S")))

    ## ===================================================
    ## Nebulosity fields

    raw_base_dir = os.path.expanduser("~/mnt/jgt/")
    img_base_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images")

    ## Details for fields taken with Nebulosity
    ## SX10 configs are too different to be able to be put into the same array
    field_details = [
            # image_dir             image_prefix    n_sets  flat_prefix bias_prefix
            #[raw_base_dir+"/2019/1028", "l135",         6,      "noflat",   "bias"],
            #[raw_base_dir+"/2021/1210Trius", "l135_5",  4,      "dflat7",   "bias2"],
            #[raw_base_dir+"/2022/0105", "l136",         9,      "dflat",    "bias2"],
            #[raw_base_dir+"/2022/0117", "l136_5",       7,      "dflat",    "bias4"],
            #[raw_base_dir+"/2022/0121", "l137_0",       7,      "dflat",    "bias_shutter"],
            #[raw_base_dir+"/2022/0124", "l137_5",       7,      "dflat",    "bias"],
            #[raw_base_dir+"/2019/0218", "l140_0",       7,      "noflat",   "BiasLast"],
            #[raw_base_dir+"/2019/0221", "l140_5",       7,      "noflat",   "Bias_end_"],
            #[raw_base_dir+"/2020/0206Trius", "l141",    7,      "flat",     "bias2"],
            #[raw_base_dir+"/2020/0212", "l141_5",       8,      "flat",     "bias2"],
            #[raw_base_dir+"/2019/0226", "l196",         7,      "noflat",     "bias_end"],
            #[raw_base_dir+"/2019/0225", "l196_5",       6,      "noflat",     "bias_end"],
            #[raw_base_dir+"/2019/0204", "l197",         7,      "noflat",     "Bias"],
            #[raw_base_dir+"/2019/0131", "l197.5",       9,      "noflat",     "biasend"],
            [raw_base_dir+"/2019/0128", "l198",         7,      "noflat",     "bias_2"],
            #[raw_base_dir+"/2019/0129", "l198_5",       7,      "noflat",     "bias-1"],
            #[raw_base_dir+"/2019/0207", "l199",         7,      "noflat",     "Bias"],
    ]

    vt = 1.3
    at = 0.85

    ## Run for all Nebulosity fields
    for raw_image_dir, image_prefix, n_sets, flat_prefix, bias_prefix in field_details:

        config = Config(
            raw_image_dir = os.path.expanduser(raw_image_dir),
            image_dir = os.path.expanduser(img_base_dir),
            image_prefix = image_prefix,
            n_sets = n_sets,
            flat_prefix = flat_prefix,
            bias_prefix = bias_prefix,
            variability_threshold = vt,
            amplitude_score_threshold = at,
            catalogue_set_number = 1,
            catalogue_image_number = 1,
        )

        #Pipeline.run(config)
        Pipeline.run_analysis(config, assume_already_adjusted=True)


    ## ===================================================
    ## SX10 fields

    ## TODO: Add l138_0 and others
    field_details = [
            # image_dir             image_prefix    n_sets  flat_prefix bias_prefix
            #[raw_base_dir+"/2022/0301", "l138_0",         9,      "noflat",   "bias_end"],
    ]

    for raw_image_dir, image_prefix, n_sets, flat_prefix, bias_prefix in field_details:

        config = Config(
            raw_image_dir = os.path.expanduser(raw_image_dir),
            image_dir = os.path.expanduser(img_base_dir),
            image_prefix = image_prefix,
            n_sets = n_sets,
            flat_prefix = flat_prefix,
            bias_prefix = bias_prefix,
            fits_extension = ".fits",
            fits_date_format = "%Y.%m.%dT%H:%M:%S.%f",
            has_filter_in_header = False,
            variability_threshold = vt,
            amplitude_score_threshold = at,
        )

        #Pipeline.run(config)
        Pipeline.run_analysis(config, assume_already_adjusted=True)

    _ = Utilities.finished_job("all fields", start_time)


## TODO: Docs
def hand_pick():
    config = Config()

    hand_picked = [
        ["l135",   1005, 1], ## Sine
        ["l135_5",  148, 0], ## Transit
        ["l135_5",  590, 0], ## Long period?
        ["l135_5",  630, 0], ## Transit
        ["l136_5",  576, 1], ## Transit or sine?
        ["l136_5",  580, 1], ## Sine
        ["l136_5",  599, 2], ## Multi-period
        ["l137_0", 1084, 1], ## Sine
        ["l138_0", 1194, 3], ## Multi-period
        ["l138_0", 1384, 0], ## Eclipse?
        ["l140_0",  760, 1], ## Sine
        #["l140_0",  762, 1], ## Shadow of above?
        ["l140_5",  479, 0], ## Sharp transit
        ["l141",    125, 1], ## Sharp bottom
        ["l141_5",  276, 2], ## Multi-period
        ["l141_5",  301, 2], ## Multi-period
        ["l141_5",  448, 1], ## Sine
        ["l196",    507, 4], ## Sine, sharp bottom
        ["l196",   1287, 0], ## Eclipse
        ["l196",   1416, 0], ## Weird transit
        ["l196",   1637, 2], ## Transit
        #["l197.5", 1235, 0], ## Doesn't pick up period
        ["l197",     98, 4], ## Multi-period, appears irregular
        ["l197",    164, 0], ## Transit
        ["l197",    371, 3], ## Multi-period
        ["l197",    456, 2], ## Multi-period
        ["l197",    869, 1], ## Sine
        ["l197",   1077, 1], ## Sine
        ["l197",   1499, 2], ## Multi-period ?
        ["l197",   1844, 0], ## Transit
        ["l197",   2118, 2], ## Irregular?
        ["l197",   2374, 1], ## Sine, long period
        ["l198",    482, 1], ## Sine
        ["l198",   1577, 1], ## Sine
        #["l198",   1593, 1], ## Sine, shadow
        ["l198",   2133, 1], ## Flare? Irregular?
        ["l198",   2268, 0], ## Scoop?
        ["l198",   2640, 0], ## Transit
    ]

    hand_picked_unsure = [
        ["l135",    121, 1], # long?
        ["l135",   1080, 1],
        ["l135_5",  392, 1],
        ["l135_5",  398, 1],
        ["l135_5",  410, 1],
        ["l135_5",  590, 1], ## 60% 
        ["l135_5",  637, 1],
        ["l135_5",  839, 1],
        ["l135_5",  904, 0],
        ["l135_5",  955, 1],
        ["l135_5", 1000, 1],
        ["l135_5", 1004, 1],
        ["l137_0",  512, 1],
        ["l137_0",  686, 1],
        ["l196",   62, 1],
        ["l196",   96, 1],
        ["l196", 1080, 1],
        ["l196", 1497, 1],
        ["l196", 1635, 1], # weird transit?
        ["l197",  444, 1],
        ["l197",  479, 1], # weird transit?
        ["l197",  841, 1],
        ["l197",  957, 1],
        ["l197", 1118, 1],
        ["l197", 1293, 1],
        ["l197", 2140, 0], # transit?
        ["l197", 2462, 1],
        ["l198_5", 1340, 1], # weird
        ["l198_5", 1411, 1], # weird
        ["l198",  745, 1],
        ["l198", 1115, 0], # double transit? Shadow?
        ["l198", 1307, 1], # ??
        ["l198", 1589, 0], # double transit?
    ]

    n_variables = len(hand_picked)
    print("[HandPicked] Running for {} variables".format(n_variables))

    fname = "results_dr2.txt"
    final_results = open(os.path.join(config.hand_picked_dir, fname), "w")
    final_results.write(
            "RA [deg]\tDEC [deg]\tField id\tsource id\tperiod [h]\tperiod error [h]\tamplitude [%]\tamplitude error [%]\n")

    for field_id, source_id, n_periods in hand_picked:
    #for field_id, source_id, n_periods in hand_picked_unsure:
        config = Config(
            image_prefix = field_id,
        )
        pf = PeriodFinder(config)

        print("[Main] Hand picking for {}, {}".format(field_id, source_id))

        ## Grab the RA and DEC of the source
        cat = Utilities.read_catalogue(config)
        idx = np.where(cat['id'] == source_id)[0][0]
        ra  = cat['RA'][idx]
        dec = cat['DEC'][idx]

        ## Grab the adjusted light curve
        lc_path = os.path.join(config.adjusted_curve_dir,
                config.source_format_str.format(source_id))
        lc = np.genfromtxt(lc_path, dtype=config.light_curve_dtype)

        ## Plot the light curve
        plt.scatter(lc['time']/3600, lc['counts'], marker="x")

        ## Get the offset of the potential sine wave
        ps = pf.period_search_curve(
                source_id, lc['time'], lc['counts'], lc['counts_err'],
                n_samples=2000)
        offset = ps[-1]
        attempted_fit = np.ones(len(lc['counts'])) * offset

        if n_periods > 0:
            for _i in range(n_periods):
                _id, P, P_err, A, A_err, phi, offset = pf.period_search_curve(
                        source_id, lc['time'], lc['counts'], lc['counts_err'],
                        n_samples=2000)

                print("[Main] Fitting period {:5}s".format(P))

                ## Subtract the fit sine curve
                fit_sine = A*np.sin(2*np.pi/P * lc['time'] + phi)
                attempted_fit += fit_sine
                lc['counts'] -= fit_sine

                ## Write the period to the file
                final_results.write("{:.6f}\t{:.6f}\t{}\t{:04}\t{:1.5f}\t{:1.5f}\t{:02.4f}\t{:02.4f}\n"
                        .format(ra, dec, field_id, source_id, P/3600, P_err/3600, 100*A, 100*A_err))

            coord = SkyCoord(ra*u.degree, dec*u.degree)
            plt.title("Light curve and fitted plot for variable at\n{}"
                    .format(coord.to_string("hmsdms", precision=0)))
            plt.plot(lc['time']/3600, attempted_fit, color="red")

        else:
            ## If we don't have a period, just write coords
            final_results.write("{:.6f}\t{:.6f}\t{}\t{:04}\t{:.0f}\t{:4.0f}\t{:02.4f}\t{:02.4f}\n"
                    .format(ra, dec, field_id, source_id, 0, 0, 0, 0))

            plt.title("Light curve for variable at\n{}"
                    .format(coord.to_string("hmsdms", precision=0)))

        plt.xlabel("Time [hours]")
        plt.ylabel("Normalised flux [arb. u.]")

        ## Write the plot to disk
        fname = "final_{}_{}{:04}{}".format(config.image_prefix, config.identifier,
                source_id, config.plot_file_extension)
        plt.savefig(os.path.join(config.hand_picked_dir, fname))
        plt.close()


    final_results.close()



## TODO: Copy flats for ones that need it?
def rename():
    """
    Rename certain problematic images/fields.
    Should only need to run once per machine (unless you delete your
    JGT images folder).

    Assumes they are copied straight from the JGT backup.
    Change the base path for each field and they should rename to a format
    that will play nice with the pipeline.
    Also if you make any fields that need renaming, please add them here.
    It will make future observers' lives easier.

    Ideal format:
    lxxx_y_s_iii.fit
    xxx: llongitude integer part
    y: longitude decimal part, 0 or 5
    s: set number, 1..n_sets
    iii: image_number, 001..050 

    """

    print("[TEST] Starting rename")

    ## l141 first set
    base_path = os.path.expanduser("F:/2020/0206Trius")
    if os.path.exists(os.path.join(base_path, "l141_001.fit")):
        for i in range(1, 51):
            os.rename(
                    os.path.join(base_path, "l141_{:03d}.fit".format(i)),
                    os.path.join(base_path, "l141_1_{:03d}.fit".format(i))
            )

    ## l140.5 -> l140_5
    bsase_path = os.path.expanduser("F:/2019/0221")
    if os.path.exists(os.path.join(base_path, "l140.5_1_001.fit")):
        for s in range(1, 8):
            for i in range(1, 51):
                os.rename(
                        os.path.join(base_path, "l140.5_{:1d}_{:03d}.fit".format(s, i)),
                        os.path.join(base_path, "l140_5_{:1d}_{:03d}.fit".format(s, i))
                )

    ## L140 -> l140_0
    base_path = os.path.expanduser("F:/2019/0218")
    if os.path.exists(os.path.join(base_path, "L140_1_001.fit")):
        for s in range(1, 8):
            for i in range(1, 51):
                os.rename(
                        os.path.join(base_path, "L140_{:1d}_{:03d}.fit".format(s, i)),
                        os.path.join(base_path, "l140_0_{:1d}_{:03d}.fit".format(s, i))
                )

    ## Make l138_0 image numbers range from 1-50
    base_path = os.path.expanduser("F:/2022/0301")
    if not os.path.exists(os.path.join(base_path, "l138_0_1_001.fits")):
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

    ## l196
    base_path = os.path.expanduser("F:/2019/0226")
    if os.path.exists(os.path.join(base_path, "l196_1-1_001.fit")):
        os.rename(os.path.join(base_path, "l196_1_001.fit"), os.path.join(base_path, "test_l196_1_001.fit"))
        for i in range(1, 51):
            file_from = os.path.join(base_path, "l196_1-1_{:03d}.fit".format(i))
            file_to   = os.path.join(base_path, "l196_1_{:03d}.fit".format(i))
            os.rename(
                    file_from,
                    file_to
            )
            print("[RENAME] '{}' -> '{}'".format(file_from, file_to))

    ## l196.5
    base_path = os.path.expanduser("F:/2019/0225")
    if os.path.exists(os.path.join(base_path, "l196_5_001.fit")):
        for i in range(1, 51):
            file_from = os.path.join(base_path, "l196_5_{:03d}.fit".format(i))
            file_to   = os.path.join(base_path, "l196_5_1_{:03d}.fit".format(i))
            os.rename(
                    file_from,
                    file_to
            )
            print("[RENAME] '{}' -> '{}'".format(file_from, file_to))


    print("[Main] Finished renaming")

if __name__ == "__main__":
    #rename()
    #main()
    hand_pick()

