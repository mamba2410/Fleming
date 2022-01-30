import os
from os.path import join


## Per-user things
#folder = "C:\\Users\\callu\\uni-git\\masters-project\\fleming\\pipeline\\"  # Pipeline root
pipeline_dir = "/home/callum/repos/uni-git/masters-project/fleming/pipeline/"
api_key = ""                                                # Astrometry.net API key
with open(join(pipeline_dir, "astrometry_api_key.txt")) as f:
    api_key = f.read().replace("\n", "")


## General directories
workspace_subdir = "workspace/"         # Where to store everything
workspace_dir = join(pipeline_dir, workspace_subdir)
image_subdir = "../images/0121/"        # Where the JGT images are
fits_extension = ".fit"
#reduced_prefix = "r_"
reduced_prefix = ""                     # Not used
catalogue_prefix = "catalogue_"         #
standard_file_extension = ".txt"        #
table_format = "ascii"                  #
time_fname = "times" + standard_file_extension          #
shift_fname = "shifts" + standard_file_extension        #
catalogue_shift_fname = "catalogue_" + shift_fname      #
flux_subdir = "fluxes/"                                 #
light_curve_subdir = "light_curves/"                    #
flux_prefix = "fluxes_"                                 #
identifier = "id"                                       #
id_map_fname = "id_mapping" + standard_file_extension   #
adjusted_curves_subdir = "adjusted_light_curves/"       #
output_subdir = "results/"                              #
moving_obj_subdir = "moving_objects/"                   #
moving_obj_fname = "moving_objects.txt"                 #
streak_subdir = "streaks/"                              #
streak_fname = "streaks.txt"                            #
job_fname = "astrometryjob.txt"                         #
line_ending = "\n"                                      #

## Tweaking constants
variability_threshold = 0.8         # units?
check_radius = 10                   # Pixels?
cosmic_threshold = 5                # units?
flux_cutoff = 100                   # Counts?
edge_limit = 50                     # units?
inner_radius = 8                    # Pixels?
outer_radius = 13                   # Pixels?
moving_obj_check_image = 30         # units?

## Per-run constants
image_width = 2432          # Doesn't change unless not using
image_height = 1600         # SX Trius on the JGT
image_prefix = "l137_0"
has_sets = True             # Images split up into n*m blocks of exposures
set_size = 50               # Images per set
n_sets = 7                  # Number of sets total

