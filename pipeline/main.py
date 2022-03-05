from pipeline import Pipeline

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
            #["~/mnt/jgt/2019/0218", "l140_0",       7,      "noflat",   "BiasLast"],
            #["~/mnt/jgt/2019/0221", "l140_5",       7,      "noflat",   "Bias_end_"],
            #["~/mnt/jgt/2020/0206Trius", "l141",    7,      "flat",     "bias2"],
            ["~/mnt/jgt/2020/0212", "l141_5",       8,      "flat",     "bias2"],
            ]

    ## Run for all fields
    for raw_image_dir, image_prefix, n_sets, flat_prefix, bias_prefix in field_details:

        config = Constants.Config(
            raw_image_dir = os.path.expanduser(raw_image_dir),
            image_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images"),
            image_prefix = image_prefix,
            n_sets = n_sets,
            flat_prefix = flat_prefix,
            bias_prefix = bias_prefix,
        )

        Pipeline.run(config)

    _ = Utilities.finished_job("everything", start_time)


if __name__ == "__main__":
    main()
