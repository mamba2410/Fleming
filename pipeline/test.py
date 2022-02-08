import Constants
import Utilities
import Cataloguer
import ShiftFinder
import FluxFinder
import DataAnalyser

import os

def main():

    config = Constants.Config(
        image_dir = "/home/callum/mnt/data/jgtdata/l137_0",
        image_prefix = "l137_0",
        has_sets = True,
        set_size = 50,
        n_sets = 2,
    )
    
    
    c = Cataloguer.Cataloguer(config)
    catalogue_image = os.path.join(
            config.image_dir,
            config.image_format_str
                .format(1, 1))

    c.catalogue(catalogue_image, solve=False)
    Utilities.print_job("cataloguing stars")
     
    #sf = ShiftFinder.ShiftFinder(config)
    #sf.get_all_shifts()
    #Utilities.print_job("finding shifts")

    ff = FluxFinder.FluxFinder(config, c.get_n_sources())

    ff.find_all_fluxes()
    ff.make_light_curves()
    Utilities.print_job("making light curves")
    

    da = DataAnalyser.DataAnalyser(config)
    
    da.get_means_and_stds(adjusted=False)
    da.get_variables(ff, adjusted=False)
    da.plot_means_and_stds(adjusted=False)

    avg_ids = da.get_ids_for_avg()
    ff.make_avg_curve(avg_ids)

    ff.divide_by_average()

    avg_fname = "{}_avg{}".format(config.image_prefix, config.standard_file_extension)
    avg_path = os.path.join(config.workspace_dir, avg_fname)
    #ff.plot_light_curve(path=avg_path, adjusted=True, show=False)
    ff.plot_avg_light_curve(avg_path, adjusted=True, show=False)
    ff.plot_all_light_curves(adjusted=True, show=False)

    Utilities.print_job("adjusting light curves")
    
    
    da.get_means_and_stds(adjusted=True)
    da.get_variables(ff, adjusted=True)
    da.plot_means_and_stds(adjusted=True)
    da.output_results()
    da.create_thumbnails(ff)
    
    Utilities.print_job("everything")



if __name__ == "__main__":
    main()
