from photutils import DAOStarFinder
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
import os
import os.path
import Constants
import Utilities
from astropy.table import Table
import matplotlib.pyplot as plt
import FluxFinder
from astroquery.astrometry_net import AstrometryNet

#remove stars within 10 pixels on edge
#add xy dimensions of image as global variable

class Cataloguer:

    config = None           # Config object
    
    n_sources = 0           # Number of sources

    means = []
    stds = []
    var_means = []
    var_stds = []
    mean_bins = []
    avgs = []
    
    id_map = []
    wcs = None


    
    ## Constructor
    def __init__(self, config):
        self.config = config


    #generate a catalogue of all the stars in the given image, and a list of 
    #the times at which each image was taken 
    def catalogue(self, image_path):

        #read in first image data
        image_data = fits.getdata(image_path, ext=0)
        
        ## TODO: Magic number
        #build a catalogue of all stars in the image
        sources = find_stars(image_data, 5)
        self.remove_stars(sources, image_data)
        self.n_sources = len(sources['id'])


        ## ??
        Utilities.make_reg_file(
                os.path.join(self.workspace_dir), self.image_prefix, sources)
        
        print("[Cataloguer] Catalogued {} objects".format(self.n_sources))
        
        #add the RA and DEC for each star to the catalogue
        self.convert_to_ra_and_dec(image_path, sources)
        
        #write the catalogue to the catalogue file
        sources.write(self.config.catalogue_path, format=config.table_format, overwrite=True)
        
        ## Clear time tile
        if(os.path.exists(self.config.time_path)):
            open(self.config.time_path, "w").close()

    
        #loop through all images within each set 
        for s in range(1, self.config.n_sets + 1):
            for i in range(1, self.config.set_size + 1):
                
                #build image filepath 
                if self.config.has_sets:
                    image_file = os.path.join(self.config.workspace_dir, self.config.format_str.format(s, i))
                else:
                    image_file = os.path.join(self.config.workspace_dir, self.config.format_str.format(i))
                
                #store the time which the current image was taken
                self.add_times(self.config.time_path, fits.getheader(image_file))
    


    ## 
    def remove_stars(self, sources, image_data):
        
        to_remove = []
                
        for i in range(len(sources)):
            
            if sources['flux'][i] > Constants.flux_cutoff or not Utilities.is_within_boundaries(sources['xcentroid'][i], sources['ycentroid'][i], len(image_data[0]), len(image_data), Constants.edge_limit):
                to_remove.append(i)
        
        for i in range(len(to_remove)):
            sources.remove_row(to_remove[i] - i)
    
            
        print("[Cataloguer] removed {} objects".format(len(to_remove)))
        



    #investigate what this is doing 
    #convert all of the source x and y positions to RA and DEC
    def convert_to_ra_and_dec(self, image_file, sources=None):
        
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
    



    #add the time of the specified image file being taken to the times file
    def add_times(self, image_header):
        
        #open time file in appending mode 
        with open(self.config.time_path, "a+") as f:
            #write date of observation to file 
            f.write("{}{}".format(image_header['DATE-OBS'], self.config.line_ending))
        
    
    
    #plot the means and standard deviations of all light curves generated
    def get_means_and_stds(self, adjusted=False):
        
        #build path of the directory in which the light curves are stored
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
        for file in os.listdir(light_curve_dir):
            if file[:len(self.config.image_prefix)] == self.config.image_prefix:

                #read light curve data from file
                t = Table.read(os.path.join(light_curve_dir, file), format=self.config.table_format)

                                
                ## TODO: Magic number
                #only plot data point if at least 5 non-zero counts are recorded
                if len(t['counts']) > 100:
                    mean = Utilities.mean(t['counts'])
                    std = Utilities.standard_deviation(t['counts'])
                    value = std/mean
                    
                    image_id = file.split(self.config.identifier)[1].split(".")[0]
                    
                    # ??? Magic numbers
                    if value > 0 and value < 2 and mean > 0.02 and mean < 80:
                        self.stds.append(value)
                        self.means.append(mean)
                
                        self.id_map.append(str(image_id))
                
                    #if Utilities.is_above_line(std, mean, 17, 0, 0.05) and mean > 50:
                    #if value < 0.01 and mean > 3:
                        #print(file.split("id")[1].split(".")[0],Utilities.mean(t['counts']))
        
        #compile results into single array
        a = [self.means, self.stds, self.id_map]
        #sort results by decreasing variability 
        Utilities.quicksort(a, True)
    



    # TODO: Plot title
    def plot_means_and_stds(self):
        plt.scatter(self.means, self.stds, marker='.')
        plt.scatter(self.var_means, self.var_stds, marker='.', color='red')
        plt.xlabel("mean")
        plt.ylabel("standard deviation")
        plt.xlim(80, 0)
        
        #ensure plot y axis starts from 0
        plt.gca().set_ylim(bottom=0)
        plt.show()
        plt.close()
    


    # TODO:
    #finds which stars are variable by comparing their standard deviations
    #and means to those of the surrounding stars in the mean-std plot
    #I think a better method should be implemented
    def is_variable(self, index):
        
        #number of stars in either direction to compare to
        check_radius = 10
        
        #threshold for being defined as variable
        variability = 2
    
        total = 0
        
        #ensures only checks for surrounding stars where they exist
        #(no indexing errors)
        llim = index - check_radius
        
        if llim < 0:
            llim = 0
        
        ulim = index + check_radius + 1
        
        if ulim > len(self.means):
            ulim = len(self.means)
        
        #find average std for surrounding stars in plot
        for i in range(llim, ulim):
            if i != index:
                total += self.stds[i]
        
        avg = total / (ulim-llim)
        
        #if this star's std is much higher than the surrounding average
        #it is a variable star
        if self.stds[index] > avg * (1 + variability):
            
            self.var_means.append(self.means[index])
            self.var_stds.append(self.stds[index])
            
            return True
        return False
        
        
            
## TODO: Remove?
#    #print plots of all variable stars in dataset             
#    def get_variables(self):
#        
#        t = 0
#        
#        ff = FluxFinder.FluxFinder("/Users/Thomas/Documents/Thomas_test/", "l198", True, 7, 50)
#    
#        for i in range(len(self.means)):
#            if self.is_variable(i):
#                t+= 1
#                #print(self.id_map[i])
#                ff.plot_light_curve(self.id_map[i], None, True)
        


    ## TODO: Magic numbers
    #find ids of stars with high brightness to produce average light curve.
    #the average light curve is subtracted from each printed light curve to 
    #reduce bias
    def get_ids_for_avg(self):
        
        ids = []
        
        for i in range(len(self.means)):
            #if self.means[i] > 50 and self.stds[i] < 0.04:
            #if not Utilities.is_above_line(-0.0001, 0.03, self.means[i], self.stds[i], 0.01) and self.means[i] > 5:
            if not Utilities.is_above_line(self.means[i], self.stds[i], 2.2222*10**-9, 0.05777778, 0.001) and self.means[i] > 10^6:
       
                fname = self.config.source_format_str.format(self.id_map[i])
                light_curve_path = os.path.join(self.workspace_dir, Constants.light_curve_dir, fname)
                t = Table.read(light_curve_path, format = Constants.table_format)
                
                if len(t['time']) == self.set_size * self.n_sets:
                    ids.append(self.id_map[i])
    
        return ids
        
        


    # TODO: Docs
    def get_wcs_header(self, file):
        ast = AstrometryNet()
        ast.TIMEOUT = 1200
        ast.api_key = self.config.api_key
    
    
        print("[Astrometry] Starting job")
        try:
            wcs = ast.solve_from_image(file, solve_timeout = 1200)
        except:
            wcs = None
            print("[Astrometry] Error: No WCS found")
        
        print('[Astrometry] Finished job')
        if not wcs:
            print('[Astrometry] Error: Job failed')
            
        return WCS(header=wcs)
            

        
## End of class


## TODO: Magic numbers
#catalogue all sources that meet the thresholds in the image
def find_stars(image_data, threshold):
    
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

                    

