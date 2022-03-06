from pipeline import Utilities, Pipeline, Config

from datetime import datetime
import os

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

    ## Details for fields taken with Nebulosity
    ## SX10 configs are too different to be able to be put into the same array
    field_details = [
            # image_dir             image_prefix    n_sets  flat_prefix bias_prefix
            #["~/mnt/jgt/2019/1028", "l135",         6,      "noflat",   "bias"],
            #["~/mnt/jgt/2021/1210Trius", "l135_5",  4,      "dflat7",   "bias2"],
            #["~/mnt/jgt/2022/0105", "l136",         9,      "dflat",    "bias2"],
            #["~/mnt/jgt/2022/0117", "l136_5",       7,      "dflat",    "bias4"],
            #["~/mnt/jgt/2022/0121", "l137_0",       7,      "dflat",    "bias_shutter"],
            #["~/mnt/jgt/2022/0124", "l137_5",       7,      "dflat",    "bias"],
            #["~/mnt/jgt/2019/0218", "l140_0",       7,      "noflat",   "BiasLast"],
            #["~/mnt/jgt/2019/0221", "l140_5",       7,      "noflat",   "Bias_end_"],
            #["~/mnt/jgt/2020/0206Trius", "l141",    7,      "flat",     "bias2"],
            ["~/mnt/jgt/2020/0212", "l141_5",       8,      "flat",     "bias2"],
            ]

    ## Run for all Nebulosity fields
    for raw_image_dir, image_prefix, n_sets, flat_prefix, bias_prefix in field_details:

        config = Config(
            raw_image_dir = os.path.expanduser(raw_image_dir),
            image_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images"),
            image_prefix = image_prefix,
            n_sets = n_sets,
            flat_prefix = flat_prefix,
            bias_prefix = bias_prefix,
        )

        Pipeline.run(config)


    ## ===================================================
    ## SX10 fields

    ## TODO: Add l138_0 and others

    _ = Utilities.finished_job("all fields", start_time)


if __name__ == "__main__":
    main()
