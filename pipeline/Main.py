import Reducer
import ShiftFinder
import FluxFinder
import DataAnalyser
import Constants
import Utilities
import MovingObjectFinder
import StreakFinder
import Cataloguer

from datetime import datetime

import os.path

now = datetime.now()
wcs = None

current_time = now.strftime("%H:%M:%S")
print("Started at " + current_time)


# =============================================================================
#r = Reducer.Reducer(Constants.folder,"No filter",Constants.file_name, 9, 50)
#r.reduce(True)
#set_size, n_sets = r.get_set_info()
#Utilities.print_job("reducing images")


c = Cataloguer.Cataloguer(Constants.workspace_dir, Constants.image_prefix, Constants.has_sets, Constants.set_size, Constants.n_sets)
c.catalogue(os.path.join(Constants.workspace_dir, Constants.image_subdir, 
        "{}_{:1}_{:03}{}".format(Constants.image_prefix, 1, 1, Constants.fits_extension)))

#mof = MovingObjectFinder.MovingObjectFinder(Constants.folder)
#mof.find_moving_objects()
 
 
Utilities.print_job("cataloguing stars")
 
sf = ShiftFinder.ShiftFinder(Constants.workspace_dir, Constants.image_prefix, 
        Constants.has_sets, Constants.set_size, Constants.n_sets)
sf.get_all_shifts()
 
Utilities.print_job("finding shifts")
 
# =============================================================================
ff = FluxFinder.FluxFinder(Constants.workspace_dir, Constants.image_prefix,
        True, Constants.set_size, Constants.n_sets)
ff.find_all_fluxes()
ff.make_light_curves()

#strf = StreakFinder.StreakFinder(Constants.folder, c, ff)
#strf.find_all_streaks()

Utilities.print_job("making light curves")

r = None
c = None
sf = None

da = DataAnalyser.DataAnalyser(Constants.workspace_dir, Constants.image_prefix, True,
        Constants.n_sets, Constants.set_size)

da.get_means_and_stds(False)
da.get_variables(False)
da.plot_means_and_stds("_adjusted")
ids = da.get_ids_for_avg()


ff.make_avg_curve(ids)
ff.divide_by_average()
avg_fname = "{}_avg{}".format(Constants.image_prefix, Constants.standard_file_extension)
ff.plot_light_curve(None, os.path.join(Constants.workspace_dir, avg_fname), True)
Utilities.print_job("adjusting light curves")


da.get_means_and_stds(True)
da.get_variables(True)
da.plot_means_and_stds("_adjusted")
da.output_results()
da.create_thumbnails(ff)

Utilities.print_job("everything")

