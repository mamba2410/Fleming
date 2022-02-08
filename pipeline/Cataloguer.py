from photutils import DAOStarFinder
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.table import Table
from astroquery.astrometry_net import AstrometryNet
import os
import matplotlib.pyplot as plt
import FluxFinder
import Utilities
import Constants


class Cataloguer:
    """
    Cataloguer class
    ----------
    """

    config = None           # Config object
    n_sources = 0           # Number of sources

    means = []
    stds = []
    var_means = []
    var_stds = []
    mean_bins = []
    avgs = []
    
    id_map = []
    solve = True        ## debug, solve astrometry_net plate
    wcs = None          ## TODO: Do something with this

    
    def __init__(self, config):
        """
        Parameters
        ----------

        config: Constants.Config
        """
        self.config = config


    def catalogue(self, image_path, solve=True):
        """
        Generates a catalogue of all the stars in the given image, as well as
        a list of the times at which each image was taken.

        Parameters
        ----------
        image_path: string
            path to the fits image to catalogue all the stars.
        solve: bool, optional
            should solve astrometry.net wcs


        """

        self.solve = solve

        #read in image data
        image_data = fits.getdata(image_path, ext=0)
        
        ## TODO: Magic number
        #build a catalogue of all stars in the image
        sources = find_stars(image_data, 5)
        self.remove_stars(sources, image_data)
        self.n_sources = len(sources['id'])

        ## Write all sources to a file
        Utilities.make_reg_file(self.config.workspace_dir, self.config.image_prefix, sources)
        
        print("[Cataloguer] Catalogued {} objects".format(self.n_sources))
        
        #add the RA and DEC for each star to the catalogue
        self.convert_to_ra_and_dec(image_path, sources)
        
        #write the catalogue to the catalogue file
        sources.write(self.config.catalogue_path, format=self.config.table_format, overwrite=True)
        
        ## Clear time file
        if(os.path.exists(self.config.time_path)):
            open(self.config.time_path, "w").close()

    
        #loop through all images within each set 
        for f, _s, _i in Utilities.loop_images(self.config):
                
                #build image filepath 
                if self.config.has_sets:
                    image_file = os.path.join(self.config.image_dir, f)
                else:
                    image_file = os.path.join(self.config.image_dir, f)
                
                #store the time which the current image was taken
                self.add_times(fits.getheader(image_file))
    


    def remove_stars(self, sources, image_data):
        """
        Removes stars from a source catalogue given cutoff parameters
        from the config object
        """
        
        to_remove = []
                
        for i in range(len(sources)):
            
            is_too_bright = sources['flux'][i] > self.config.flux_cutoff
            is_within_boundaries = Utilities.is_within_boundaries(
                    sources['xcentroid'][i],
                    sources['ycentroid'][i],
                    len(image_data[0]),
                    len(image_data),
                    self.config.edge_limit)

            if is_too_bright or not is_within_boundaries:
                to_remove.append(i)
        

        ## Remove rows from sources table
        for i in range(len(to_remove)):
            sources.remove_row(to_remove[i] - i)
    
            
        print("[Cataloguer] filtered out {} objects".format(len(to_remove)))
        



    def convert_to_ra_and_dec(self, image_file, sources=None):
        """
        Converts all the source positions in the image to RA and DEC.
        Uses World Coordinate System (wcs) from the fits image.
        """
        
        print("[Cataloguer] Getting coordinate system for image '{}'".format(image_file))

        # find the wcs assosiated with the fits image using astropy and the header
        self.wcs = self.get_wcs_header(image_file)
        
        ## make two new coloums of 0's
        #sources['RA'] = sources['xcentroid'] * 0
        #sources['DEC'] = sources['xcentroid'] * 0
    
        ## replace the 0's with ra and dec
        #for x in range(0,len(sources)):
        #    ra, dec = self.wcs.all_pix2world(sources['xcentroid'][x], sources['ycentroid'][x], 0) 
        #    sources['RA'][x] = ra
        #    sources['DEC'][x] = dec
    



    def add_times(self, image_header):
        """
        Add the time of the specified image file being taken to the times file

        Parameters
        ----------

        image_header: fits header for image to add
        """
        
        #open time file in appending mode 
        with open(self.config.time_path, "a+") as f:
            #write date of observation to file 
            f.write("{}{}".format(image_header['DATE-OBS'], self.config.line_ending))
        


    
    
    def get_means_and_stds(self, adjusted=False):
        """
        Find the means and standard deviations of all light curves generated.

        Parameters
        ==========

        adjusted: bool
            Has the light curves been adjusted (divided by average)

        """
        
        self.means = []
        self.stds = []
        
        cat = Table.read(self.config.catalogue_path, format=self.config.table_format)
        
        if adjusted:
            light_curve_dir = self.config.adjusted_curve_dir
            print("[Cataloguer] Using adjusted light curves")
        else:
            light_curve_dir = self.config.light_curve_dir
            print("[Cataloguer] Using non-adjusted light curves")
    

        #for each image file in the light curve directory 
        for path, source_id in Utilities.list_sources(self.config):
    
            #read light curve data from file
            t = Table.read(path, format=self.config.table_format)
            
                            
            ## TODO: Magic number
            #only plot data point if at least 5 non-zero counts are recorded
            if len(t['counts']) > 100:
                mean = Utilities.mean(t['counts'])
                std = Utilities.standard_deviation(t['counts'])
                value = std/mean
                
                # ??? Magic numbers
                if value > 0 and value < 2 and mean > 0.02 and mean < 80:
                    self.stds.append(value)
                    self.means.append(mean)
            
                    self.id_map.append(str(source_id))
            
                #if Utilities.is_above_line(std, mean, 17, 0, 0.05) and mean > 50:
                #if value < 0.01 and mean > 3:
                    #print(file.split("id")[1].split(".")[0],Utilities.mean(t['counts']))
        
        #compile results into single array
        a = [self.means, self.stds, self.id_map]
        #sort results by decreasing variability 
        Utilities.quicksort(a, True)

        ## TODO: Do something with this array
    



    def plot_means_and_stds(self, show=False):
        """
        Plot the means and standard deviations of all light curves generated.

        Parameters
        ----------

        show: bool, optional. Show plot to user.
        """

        plt.scatter(self.means, self.stds, marker='.')
        plt.scatter(self.var_means, self.var_stds, marker='.', color='red')
        plt.xlabel("Mean [counts]")
        plt.ylabel("Standard deviation [counts]")
        plt.title("Mean against standard deviation for each source in catalogue {}"
                .format(self.config.image_prefix))

        #ensure plot y axis starts from 0
        plt.xlim(80, 0)
        plt.gca().set_ylim(bottom=0)
        

        # Save image to file
        fname = "{}_mean_std{}".format(self.config.image_prefix, self.config.plot_file_extension)
        path = os.path.join(self.config.output_dir, fname)
        plt.savefig(path)

        if show:
            plt.show()

        plt.close()
    

    ## TODO: Magic numbers, make this more config friendly
    def get_ids_for_avg(self):
        """
        Find IDs of stars with high brightness to produce average light curve.
        The average light curve is subtracted from each printed light curve to reduce bias
        """
        print("[DEBUG] Calling `get_ids_for_avg` in Cataloguer")
        
        ids = []
        
        for i in range(len(self.means)):
            if not Utilities.is_above_line(self.means[i], self.stds[i], 2.2222*10**-9, 0.05777778, 0.001) and self.means[i] > 10^6:
       
                fname = self.config.source_format_str.format(self.id_map[i])
                light_curve_path = os.path.join(self.config.light_curve_dir, fname)
                t = Table.read(light_curve_path, format=self.config.table_format)
                
                # Add it if we have all time points
                # TODO: Add a >= 0.8* n?
                if len(t['time']) == self.config.set_size * self.config.n_sets:
                    ids.append(self.id_map[i])
    
        return ids
        
        


    ## TODO: Add threshold and FWHM and make config variables
    def get_wcs_header(self, file):
        """
        Get World Coordinate System header.
        [Astroquery docs.](https://astroquery.readthedocs.io/en/latest/astrometry_net/astrometry_net.html)

        Parameters
        ----------
        file: string
            File to send to astrometry.net

        Returns
        -------
        WCS
            FITS world coordinate system
        """

        ast = AstrometryNet()
        ast.TIMEOUT = self.config.astrometry_timeout
        ast.api_key = self.config.api_key
    
        wcs = None
        if not self.solve:
            return WCS(header=wcs)
    
        print("[Astrometry] Starting job")
        try:
            ## Solve on local machine
            wcs = ast.solve_from_image(
                    file, 
                    solve_timeout=self.config.astrometry_timeout,
                    force_image_upload=False)
        except Exception as e:
            print("\n[Astrometry] Error: WCS failed")
            print("............ Exception: {}".format(e))
        
        print("\n[Astrometry] Finished job")
            
        ## TODO: Gives warning with no image data
        return WCS(header=wcs)
            

    def get_n_sources(self):
        return self.n_sources

        
## End of class


## TODO: Magic numbers, make configs
def find_stars(image_data, threshold):
    """
    Catalogue all sources that meet the thresholds in the image.

    Parameters
    ----------
    image_data: fits image data
        Image data to find stars in.

    threshold: float
        Number of standard deviations above which is a star.

    """
    
    #get mean median and standard deviation of the image data
    mean, median, std = sigma_clipped_stats(image_data, sigma=3.0, maxiters=5)  
    
    #initiate finder object. Will find objects with a FWHM of 8 pixels
    # and 3-sigma times the background
    daofind = DAOStarFinder(fwhm=8, threshold=threshold*std) 
    
    #finds sources
    sources = daofind(image_data)

        
    for col in sources.colnames:    
        sources[col].info.format = '%.8g'  # for consistent table output
       
    return sources

                    

