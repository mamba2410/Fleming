import os


class Config:
    
    ## Defaults for the config.
    ## These shouldn't need changing here
    ## Set these when creating the constructor in Main.py
    def __init__(self,
            pipeline_root       = ".",
            workspace_subdir    = "workspace",
            image_subdir        = "../images",
            flux_subdir         = "fluxes",
            light_curve_subdir  = "light_curves",
            adjusted_curve_subdir = "adjusted_light_curves",
            output_subdir       = "results",
            streak_subdir       = "streaks",
            moving_obj_subdir   = "moving_objects",

            time_fname          = "times",
            shift_fname         = "shift",
            id_map_fname        = "id_mapper",
            moving_obj_fname    = "moving_obj",
            astrometry_job_fname = "astrometryjob",

            reduced_prefix      = "",
            catalogue_prefix    = "catalogue_",
            flux_prefix         = "flux_",

            fits_extension      = ".fit",
            standard_file_extension = ".txt",
            plot_file_extension = ".jpg",

            line_ending         = "\n",
            identifier          = "id",
            table_format        = "ascii",

            variability_threshold   = 0.8,
            check_radius            = 10,
            cosmic_threshold        = 5,
            flux_cutoff             = 100,
            edge_limit              = 50,
            inner_radius            = 8,
            outer_radius            = 13,
            moving_obj_check_image  = 30,

            image_width         = 2432,
            image_height        = 1600,
            image_prefix        = "l137_0",
            has_sets            = True,
            set_size            = 50,
            n_sets              = 7,

            astrometry_timeout  = 1200,
            api_key             = "",
            api_key_file        = "astrometry_api_key.txt"
            ):

        self.pipeline_root          = pipeline_root
        self.workspace_dir          = os.path.join(pipeline_root, workspace_subdir)
        self.image_dir              = os.path.join(self.workspace_dir, image_subdir)
        self.flux_dir               = os.path.join(self.workspace_dir, flux_subdir)
        self.light_curve_dir        = os.path.join(self.workspace_dir, light_curve_subdir)
        self.adjusted_curve_dir     = os.path.join(self.workspace_dir, adjusted_curve_subdir)
        self.output_dir             = os.path.join(self.workspace_dir, output_subdir)
        self.streak_dir             = os.path.join(self.workspace_dir, streak_subdir)
        self.moving_obj_dir         = os.path.join(self.workspace_dir, moving_obj_subdir)

        fname = "{}{}".format(time_fname, standard_file_extension)
        self.time_path              = os.path.join(self.workspace_dir, fname)
        fname = "{}{}".format(shift_fname, standard_file_extension)
        self.shift_path             = os.path.join(self.workspace_dir, fname)
        fname = "{}{}".format(id_map_fname, standard_file_extension)
        self.id_map_path            = os.path.join(self.workspace_dir, fname)
        fname = "{}{}".format(moving_obj_fname, standard_file_extension)
        self.moving_obj_path        = os.path.join(self.workspace_dir, fname)
        fname = "{}{}".format(astrometry_job_fname, standard_file_extension)
        self.astrometry_job_path    = os.path.join(self.workspace_dir, fname)

        self.reduced_prefix         = reduced_prefix
        self.catalogue_prefix       = catalogue_prefix
        self.flux_prefix            = flux_prefix

        self.fits_extension             = fits_extension
        self.standard_file_extension    = standard_file_extension
        self.plot_file_extension        = plot_file_extension

        self.line_ending        = line_ending
        self.identifier         = identifier
        self.table_format       = table_format

        self.variability_threshold  = variability_threshold
        self.check_radius           = check_radius
        self.cosmic_threshold       = cosmic_threshold
        self.flux_cutoff            = flux_cutoff
        self.edge_limit             = edge_limit
        self.inner_radius           = inner_radius
        self.outer_radius           = outer_radius
        self.moving_obj_check_image = moving_obj_check_image

        self.image_width        = image_width
        self.image_height       = image_height
        self.image_prefix       = image_prefix
        self.has_sets           = has_sets
        self.set_size           = set_size
        self.n_sets             = n_sets

        self.astrometry_timeout = astrometry_timeout
        self.api_key_file       = api_key_file
        self.api_key            = api_key


        ## ===============
        ## Setting up other stuff

        with open(os.path.join(pipeline_root, api_key_file)) as f:
            self.api_key = f.read().replace(line_ending, "")


        fname = "{}{}{}".format(catalogue_prefix, image_prefix, standard_file_extension)
        self.catalogue_path = os.path.join(self.workspace_dir, fname)

        if self.has_sets:
            self.image_format_str = "{}{}_{}_{}{}".format(
                    reduced_prefix, image_prefix, "{:1}", "{:03}", fits_extension)
        else:
            self.image_format_str = "{}{}_{}{}".format(
                    reduced_prefix, image_prefix, "{:04}", fits_extension)

        self.source_format_str      = "{}_{}{}{}".format(
                image_prefix, identifier, "{:04}", standard_file_extension)


        ## ===============
        ## Create directories if they don't exist

        if not os.path.exists(self.workspace_dir):
            os.mkdir(self.workspace_dir)

        if not os.path.exists(self.image_dir):
            os.mkdir(self.image_dir)

        if not os.path.exists(self.flux_dir):
            os.mkdir(self.flux_dir)

        if not os.path.exists(self.light_curve_dir):
            os.mkdir(self.light_curve_dir)

        if not os.path.exists(self.adjusted_curve_dir):
            os.mkdir(self.adjusted_curve_dir)

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

        if not os.path.exists(self.streak_dir):
            os.mkdir(self.streak_dir)

        if not os.path.exists(self.moving_obj_dir):
            os.mkdir(self.moving_obj_dir)



