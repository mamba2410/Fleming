import os
from astropy.table import Table
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import CircularAnnulus
from astropy.io import fits
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import Constants

import Utilities

class FluxFinder:

    config  = None
    
    #list of all total shifts
    x_shifts = None
    y_shifts = None
    
    #catalogue of image 1
    catalogue = None
    n_sources = None
    n_measures = None
    
    
    def __init__(self, config, n_sources):
        """
        Initialises class and reads shifts and catalogue from files.
        """

        self.config = config
        self.n_sources = n_sources
        self.n_measures = self.config.set_size*self.config.n_sets
        
        self.catalogue = Utilities.read_catalogue(self.config)
        self.x_shifts, self.y_shifts = np.genfromtxt(self.config.shift_path).transpose()

        ## TODO: Sanity check if catalogue has same number of sources
        ## as arg passed



    def get_total_shift(self, set_number=0, image_number=0):
        """
        Get the total shift between the specified image and the first image.
        Requires `has_sets=True`

        Parameters
        ----------

        set_number: int
            Set of the image

        image_number: int
            Number of image in set

        """
    
        #file is just one long list of shifts, with a shift associated with 
        #each image - below is the expression required
        #to get the shift for the image belong to set number 'set' and image 
        #number 'image_number'
        index = (set_number-1)*self.config.set_size + image_number-1
        
        total_x = self.x_shifts[index]
        total_y = self.y_shifts[index]

        return total_x, total_y
    
    

    def find_fluxes(self, set_number=0, image_number=0):
        """
        Find the fluxes of all stars in a given image

        Parameters
        ----------

        set_number: int
            Set of the image

        image_number: int
            Number of image in set

        """

        print("[FluxFinder] Finding fluxes in image: set {:1}; image: {:03}"
                .format(set_number, image_number))

        #get total shift from first image to this one
        x_shift, y_shift = self.get_total_shift(set_number=set_number, image_number=image_number)
        
        #add the shift onto the positions of the stars in the first image
        #to find their positions in this image
        x = self.catalogue['xcentroid'] + x_shift
        y = self.catalogue['ycentroid'] + y_shift
        

        positions = np.array([x, y]).transpose()
        
        image_data = Utilities.get_image(self.config, set_number, image_number)[0].data
        x_max = self.config.image_width
        y_max = self.config.image_height

        # Local background subtraction

        # Define size of background aperture
        background_annuli = CircularAnnulus(
                positions,
                r_in=self.config.inner_radius,
                r_out=self.config.outer_radius)
        star_apertures = CircularAperture(positions, r=self.config.inner_radius-1)
        all_apertures = [star_apertures, background_annuli]

        # find counts in each aperture and annulus
        phot_table2 = aperture_photometry(image_data, all_apertures)

        n_sources_phot = len(phot_table2['id'])

        if n_sources_phot != self.n_sources:
            print("[FluxFinder] Warn: Aperture photometry found {} sources, self has {}"
                    .format(n_sources_phot, self.n_sources))

        for col in phot_table2.colnames:
            phot_table2[col].info.format = '%.8g'  # for consistent table output
    

        ## Set new columns to some unphysical value
        phot_table2['mean'] = phot_table2['aperture_sum_1']*0 - 1e9
        phot_table2['median'] = phot_table2['aperture_sum_1']*0 - 1e9
        phot_table2['residual_aperture_sum_mean'] = phot_table2['aperture_sum_1']*0 - 1e9
        phot_table2['residual_aperture_sum_med'] = phot_table2['aperture_sum_1']*0 - 1e9

        for i in range(n_sources_phot):
            
            x, y = positions[i]

            #check that largest aperture does not exceed the boundaries of the image
            if Utilities.is_within_boundaries(x, y, x_max, y_max, self.config.outer_radius):

                # Calc the mean background in the second aperture ring
                bkg_mean = phot_table2['aperture_sum_1'][i] / background_annuli[i].area
                phot_table2['mean'][i] = bkg_mean
                
                # Calc background level in each main aperture and subtract
                bkg_sum = bkg_mean * star_apertures[i].area
                final_sum = phot_table2['aperture_sum_0'][i] - bkg_sum
                phot_table2['residual_aperture_sum_mean'][i] = final_sum

                # Calc background median
                ann_mask = background_annuli[i].to_mask(method='center')
                weighted_data = ann_mask.multiply(image_data)
                phot_table2['median'][i] = np.median(weighted_data[weighted_data != 0])

                # Calc 
                bkg_med = phot_table2['median'][i]
                bkg_sum = bkg_med * star_apertures[i].area
                final_sum = phot_table2['aperture_sum_0'][i] - bkg_sum
                phot_table2['residual_aperture_sum_med'][i] = final_sum

            ## For debug purposes
            #else:
            #    print("[FluxFinder] Source {} is not within boundary {},{}; {},{}"
            #            .format(i, x_max, y_max, x, y))

        phot_table2['residual_aperture_sum_mean'].info.format = '%.8g'  # for consistent table output
        phot_table2['residual_aperture_sum_med'].info.format = '%.8g'  # for consistent table output
                
        return phot_table2
        


    ## TODO: Breaks if no sets
    def make_light_curves(self):
        if self.config.has_sets == False:
            print("[FluxFinder] Error: Cannot make light curves, not implemented for has_sets=False")
            exit()
        
        times = [line.rstrip(self.config.line_ending) for line in open(self.config.time_path)]
        n_times = len(times)

        dt_obs_start = datetime.strptime(times[0], "%Y-%m-%dT%H:%M:%S")

        ## Light curves
        ## Dim 0: each source
        ## Dim 1: time/counts
        ## Dim 2: measurement number
        light_curves = np.zeros((self.n_sources, 2, n_times))
          
        ## iterate over each image
        ## image_index iterates over each measurement
        ## j iterates over source index
        image_index = 0
        for _, s, i in Utilities.loop_images(self.config):
            #print("[FluxFinder] Making light curve for image: set {:1}; image: {:03}".format(s, i))
            
            source_table = self.find_fluxes(s, i)
                            
            #iterate through each source in the table
            for j in range(self.n_sources):
                    
                #get time matching the image
                date_and_time = times[image_index]

                ## TODO: Make time format a config param
                ## Use datetime objects to calc time difference
                dt_image = datetime.strptime(date_and_time, "%Y-%m-%dT%H:%M:%S")
                dt_elapsed = dt_image - dt_obs_start
                seconds_elapsed = dt_elapsed.seconds
                

                #if median has reasonable value
                #if source_table['median'][j] > 0: 
                light_curves[j][1][image_index] = float(source_table['residual_aperture_sum_med'][j])

                #else:
                #    print("[FluxFinder] Rejecting value for source {}, median is {}"
                #            .format(t['id'][j], t['median'][j]))

                ## Include time, regardless of value
                light_curves[j][0][image_index] = seconds_elapsed
            image_index += 1

        #loop through all light curves, writing them out to a file
        for j in range(self.n_sources):
            #build light curve path
            fname = self.config.source_format_str.format(self.catalogue['id'][j])
            path = os.path.join(self.config.light_curve_dir, fname)

            ## Write light curve
            ## Transpose to have columns
            ## Only one light curve per file
            np.savetxt(path, light_curves[j].transpose())
    
            
        
        
        
    def plot_light_curve(self, source_id=None, curve_path=None, adjusted=False, show=False, close=True):
        """
        Plot light curve of star with the given ID from catalogue

        Parameters
        ----------

        source_id: int
            ID of the source to plot light curve for.
            Plots average light curve if `None`

        curve_path: string, optional
            Path to save the image to

        adjusted: bool, optional
            Has the lighth curve been divided by the average flux?

        show: bool, optional
            Should we show the plot?

        close: bool, optional
            Should we close the plot? (useful for overplotting)

        """

        print("[FluxFinder] Plotting light curve for source {} (adjusted={})"
                .format(source_id, adjusted))

        
        if curve_path == None:
            fname = self.config.source_format_str.format(source_id)
            if adjusted:
                curve_path = os.path.join(self.config.adjusted_curve_dir, fname)
            else:
                curve_path = os.path.join(self.config.light_curve_dir, fname)

        #table = Table.read(curve_path, format=self.config.table_format)
        curve = np.genfromtxt(curve_path, dtype=[
            ('time', 'float64'),
            ('counts', 'float64'),
            ]).transpose()
        
        times = curve['time']
        fluxes = curve['counts']
        

        mean = np.mean(fluxes)
        median = np.median(fluxes)
        #normalised_fluxes = fluxes / mean
        normalised_fluxes = fluxes/median
        
        minimum = min(normalised_fluxes)
        maximum = max(normalised_fluxes)
    
        ## TODO: Multiple axes, seconds and minutes?
        #plt.figure(figsize = (12, 8))
        plt.scatter(times, normalised_fluxes, s=5, marker='x');
        plt.xlabel("Time [seconds]")
        plt.ylabel("Relative flux [counts/mean]")

        #plt.ylim(0.9 * minimum, maximum * 1.1)
        ## Debug
        #plt.ylim(0, 2)
        
        if source_id == None:
            fname = "{}_LC_avg.jpg".format(self.config.image_prefix)
            plt.title("Light curve for average of sources (adjusted={})"
                .format(adjusted))
        elif source_id >= 0:
            fname = "{}_LC_{}{:04}.jpg".format(
                self.config.image_prefix, self.config.identifier, source_id)
            plt.title("Light curve for source {:04} (adjusted={})"
                .format(source_id, adjusted))
        else:
            print("[FluxFinder] Error: Cannot plot light curve, source id '{}' invalid"
                    .format(source_id))
            plt.title("Light curve for unknown source (id {}) (adjusted={})"
                .format(source_id, adjusted))
            
        if show:
            plt.show()

        image_file = os.path.join(self.config.output_dir, fname)
        plt.savefig(image_file)

        if close:
            plt.close()


        
    def plot_all_light_curves(self, ids, adjusted=False, show=False):
        for _i, path, source_id in Utilities.loop_variables(self.config, ids, adjusted=adjusted):
            self.plot_light_curve(
                    source_id=source_id, curve_path=path, adjusted=adjusted, show=show, close=True)

    def plot_adjusted_comparison(self, ids, show=False):
        for _i, _path, source_id in Utilities.loop_variables(self.config, ids, adjusted=False):
            self.plot_light_curve(
                    source_id=source_id, curve_path=None, adjusted=False, show=False, close=False)
            self.plot_light_curve(
                    source_id=source_id, curve_path=None, adjusted=True, show=show, close=True)

    def plot_avg_light_curve(self, curve_path, adjusted=False, show=False):
        self.plot_light_curve(source_id=None, curve_path=curve_path, adjusted=adjusted, show=show)


    

    ## TODO: Duplicate somewhere else
    #use median
    def make_avg_curve(self, ids):
        """
        Makes an average light curve with the brightest stars to subtract noise
        from the resulting curves

        Parameters
        ----------
        ids: [int]
            IDs of sources to take average of
        """
        print("[DEBUG] Calling `make_avg_curve` in FluxFinder")
                
        #iterate through each star in the given list
        for i in range(len(ids)):
            source_id = ids[i]
            fname = self.config.source_format_str.format(source_id)
            path = os.path.join(self.config.light_curve_dir, fname)
            #print("source_id: {}".format(source_id))
            
            #get the light curve for star i
            t = Table.read(path, format=self.config.table_format)
            
            fluxes = t['counts']
            
            if i == 0:
                for j in range(len(fluxes)):
                    self.avg_fluxes.append(0)
                self.times = t['time']

            #find the mean flux of star i
            mean = Utilities.mean(fluxes)
            
            #add each flux measurement in the light curve at time index k 
            #of star i to a list which will store the average flux at each point
            #in time k. Note the flux added to the list is divided by the mean
            #to reduce any bias(?). The average contribution to each k is 1.
            for k in range(len(fluxes)):
                self.avg_fluxes[k] = self.avg_fluxes[k] + (fluxes[k] / mean)
            
        #calculate average flux at each time of measurement l
        for l in range(len(self.avg_fluxes)):
            self.avg_fluxes[l] = self.avg_fluxes[l] / len(ids)
                
        #generate table of average light curve
        light_curve = Table([self.times, self.avg_fluxes], names = ('time','counts') )

        #save average light curve
        fname = "{}_avg{}".format(self.config.image_prefix, self.config.standard_file_extension)
        path = os.path.join(self.config.workspace_dir, fname)

        light_curve.write(path, format=self.config.table_format, overwrite=True)



    def create_adjusted_light_curves(self):
        """
        Divides all light curves by the average light curve to remove noise.
        Also removes global effects of the field that vary over time.

        Creates an 'adjusted' light curve.

        """
        print("[DEBUG] Calling `divide_by_average` in FluxFinder")

        avg_curve = np.genfromtxt(self.config.avg_curve_path, dtype=[
            ('time', 'float64'),
            ('counts', 'float64'),
            ]).transpose()

        #for all files
        for path, source_id in Utilities.list_sources(self.config, adjusted=False):
            curve = np.genfromtxt(path, dtype=[
                ('time', 'float64'),
                ('counts', 'float64'),
                ]).transpose()

            
            curve['counts'] /= avg_curve['counts']
            med = np.median(curve['counts'])
            curve['counts'] /= med
            
            fname = self.config.source_format_str.format(source_id)
            out_path = os.path.join(self.config.adjusted_curve_dir, fname)

            np.savetxt(out_path, curve)

        

    ## TODO: Use set_number and image_number
    def get_thumbnail(self, image_n, x, y, size, add_shift=True):
        """
        Get a thumbnail of a source, returns pixel values.

        Parameters
        ----------

        image_n: int
            Total image number, accumulated across sets
        size: int
            Radius of box to select around center.

        """
        
        n = image_n % self.config.set_size
        
        if n == 0:
            n = self.config.set_size
            
        set_number = int((image_n-1) / self.config.set_size) + 1
        if add_shift:
            x_shift, y_shift = self.get_total_shift(set_number=set_number, image_number=n)
        else:
            x_shift = 0
            y_shift = 0
        
        self.has_sets = True
        image_file = Utilities.get_image(self.config, set_number=set_number, image_number=n)
        
        image = image_file[0].data
        
        
        x = x + x_shift
        y = y + y_shift
        
        
        ly = int(y-size)
        ry = int(y+size)
        
        lx = int(x-size)
        rx = int(x+size)
        
        if lx < 0:
            lx = 0
        
        if rx > self.config.image_width:
            rx = self.config.image_width-1
        
        if ly < 0:
            ly = 0
        
        if ry > self.config.image_height:
            ry = self.config.image_height - 1
            
            
        return image[ly:ry, lx:rx] # note your x,y coords need to be an int



    def map_id(self, id2, cat1, cat2, shifts, set_number):
        """
        Find ID of a star from catalogue 2 in catalogue 1
        Expects everything in Table objects

        """
        
        x_shift = 0
        y_shift = 0
        
        #find total shift between catalogues
        for i in range(set_number-1):
            x_shift += shifts['xshifts'][i]
            y_shift += shifts['yshifts'][i]
        
        #find expected position of star in catalogue 1 with id 'id2' in
        #the second catalogue
        expected_x = cat2['xcenter'][id2] - x_shift
        expected_y = cat2['ycenter'][id2] - y_shift
        
        #get x-y coords of stars in catalogue 1
        x1s = cat1['xcentroid']
        y1s = cat1['ycentroid']
        
        #find matching star id in catalogue 1 - error of 0.01 seems very large
        #perhaps implement iterative search? (increasing error if required)
        for i in range(len(x1s)):
            if Utilities.equals(expected_x, x1s[i], 0.01) and Utilities.equals(expected_y, y1s[i], 0.01):
                return i
            
        return None

