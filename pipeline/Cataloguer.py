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


    """

    config = None           # Config object
    n_sources = 0           # Number of sources
    wcs = None          ## TODO: Do something with this

    
    def __init__(self, config):
        """
        Parameters
        ----------

        config: Constants.Config
        """
        self.config = config


    def generate_catalogue(self, image_path, solve=True):
        """
        Generates a catalogue of all the stars in the given image, as well as
        a list of the times at which each image was taken.

        Parameters
        ----------
        image_path: string
            path to the fits image to catalogue all the stars.
        solve: bool, optional
            should solve astrometry.net wcs?


        """

        #read in image data
        image_data = fits.getdata(image_path, ext=0)
        
        ## TODO: Magic number, threshold passed to find_stars
        #build a catalogue of all stars in the image
        sources = find_stars(image_data, 5)
        self.filter_stars(sources, image_data)
        self.n_sources = len(sources['id'])

        print("[Cataloguer] Found {} suitable sources".format(self.n_sources))

        ## Write all sources to a file
        ## TODO: Necessary?
        #Utilities.make_reg_file(self.config.workspace_dir, self.config.image_prefix, sources)
        
        #add the RA and DEC for each star to the catalogue
        self.convert_to_ra_and_dec(image_path, sources, solve=solve)
        
        #write the catalogue to the catalogue file
        sources.write(self.config.catalogue_path, format=self.config.table_format, overwrite=True)

        return self.n_sources

    


    def filter_stars(self, sources, image_data):
        """
        Removes stars from a source catalogue given cutoff parameters
        from the config object

        sources: Table
            
        image_data: 
            
        """
        
        to_remove = []
                
        for i in range(self.n_sources):
            
            is_too_bright = sources['flux'][i] > self.config.flux_cutoff
            is_within_boundaries = Utilities.is_within_boundaries(
                    sources['xcentroid'][i],
                    sources['ycentroid'][i],
                    self.config.image_width,
                    self.config.image_height,
                    self.config.edge_limit)

            if is_too_bright or not is_within_boundaries:
                to_remove.append(i)
        

        ## Remove rows from sources table
        for i in range(len(to_remove)):
            sources.remove_row(to_remove[i] - i)
    
            
        print("[Cataloguer] Filtered out {} objects".format(len(to_remove)))
        



    ## TODO: Add RA/DEC to sources table
    ##       and have error checking for astrometry wcs
    def convert_to_ra_and_dec(self, image_file, sources, solve=True):
        """
        Converts all the source positions in the image to RA and DEC.
        Uses World Coordinate System (wcs) from the fits image.

        """
        
        print("[Cataloguer] Getting coordinate system for image '{}'".format(image_file))

        # find the wcs assosiated with the fits image using astropy and the header
        self.wcs = self.get_wcs_header(image_file, solve=solve)
        
        ## make two new coloums of 0's
        #sources['RA'] = sources['xcentroid'] * 0
        #sources['DEC'] = sources['xcentroid'] * 0
    
        ## replace the 0's with ra and dec
        #for x in range(0,len(sources)):
        #    ra, dec = self.wcs.all_pix2world(sources['xcentroid'][x], sources['ycentroid'][x], 0) 
        #    sources['RA'][x] = ra
        #    sources['DEC'][x] = dec
    



    
    ## TODO: Add threshold and FWHM and make config variables
    def get_wcs_header(self, file, solve=True):
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
        if not solve:
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




    def generate_image_times(self):
        """
        Get a list of times of each measurement

        """

        print("[Cataloguer] Generating image time file")
        ## Clear time file if it exists
        if(os.path.exists(self.config.time_path)):
            open(self.config.time_path, "w").close()

        ## Open time file
        with open(self.config.time_path, "a+") as time_file:
    
            #loop through all images within each set 
            for image_name, _s, _i in Utilities.loop_images(self.config):
            
                image_file = os.path.join(self.config.image_dir, image_name)
                image_header = fits.getheader(image_file)

                ## Write time to file
                time_file.write("{}{}".format(
                            image_header['DATE-OBS'],
                            self.config.line_ending))

        
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

                    

