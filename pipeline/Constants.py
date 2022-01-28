
## General directories
working_directory = "workspace/"        # Where to store everything
image_directory = "reduced_images/jgt/0121/"     # Where the JGT images are
fits_extension = ".fit"
#reduced_prefix = "r_"
reduced_prefix = ""                     # Not used
catalogue_prefix = "catalogue_"         #
standard_file_extension = ".txt"        #
table_format = "ascii"                  #
time_file = "times" + standard_file_extension       #
shift_file = "shifts" + standard_file_extension     #
catalogue_shift_file = "catalogue_" + shift_file    #
flux_directory = "fluxes/"                          #
light_curve_directory = "light_curves/"             #
flux_prefix = "fluxes_"                             #
identifier = "id"                                   #
id_map_file = "id_mapping" + standard_file_extension    #
adjusted_curves_directory = "adjusted_light_curves/"    #
output_directory = "results/"                       #
moving_obj_folder = "moving_objects/"               #
moving_obj_file = "moving_objects.txt"              #
streak_folder = "streaks/"                          #
streak_file = "streaks.txt"                         #
job_file_name = "astrometryjob.txt"                 #

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
file_name = "l137_0"
has_sets = True             # Images split up into n*m blocks of exposures
set_size = 50               # Images per set
n_sets = 7                  # Number of sets total

## Per-user things
#folder = "C:\\Users\\callu\\uni-git\\masters-project\\fleming\\pipeline\\"  # Pipeline root
folder = "/home/callum/repos/uni-git/masters-project/fleming/pipeline/"
api_key = ""                                                # Astrometry.net API key
with open(folder+"astrometry_api_key.txt") as f:
    api_key = f.read().replace("\n", "")
api_key = "rozofdslppcdpqku"