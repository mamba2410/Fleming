import os
from astropy.table import Table
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import CircularAnnulus
from astropy.io import fits
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
    
    #time at which first image was taken
    obs_start_time = []
    
    times = []
    
    avg_fluxes = []
    
    
    def __init__(self, config, n_sources):
        """
        Initialises class and reads shifts and catalogue from files.
        """

        self.config = config
        self.n_sources = n_sources
        
        self.catalogue = Table.read(self.config.catalogue_path, format=config.table_format)

        t = Table.read(self.config.shift_path, format=self.config.table_format)
        self.x_shifts = t['xshifts']
        self.y_shifts = t['yshifts']

        


    ## TODO: Remove
    def find_all_fluxes(self):
        """
        Find the fluxes of each star in each image
        Nothing happens with the output
        """
        
        #iterate over each image
        for _, s, i in Utilities.loop_images(self.config):
            print("[FluxFinder] Finding fluxes for image: set {:1}; image {:03}".format(s, i))
            self.find_fluxes(s, i)
                
    

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
    
    

    ## TODO: Magic numbers
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

        #get total shift from first image to this one
        x_shift, y_shift = self.get_total_shift(set_number=set_number, image_number=image_number)
        
        #add the shift onto the positions of the stars in the first image
        #to find their positions in this image
        x = self.catalogue['xcentroid'] + x_shift
        y = self.catalogue['ycentroid'] + y_shift
        
        positions = []
        
        for i in range(len(x)):
            positions.append((x[i], y[i]))
            
        
        image_data = Utilities.get_image(self.config, set_number, image_number)[0].data
        x_max = len(image_data[0])
        y_max = len(image_data)

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

        for col in phot_table2.colnames:
            phot_table2[col].info.format = '%.8g'  # for consistent table output
    

        ## Clear columns
        phot_table2['mean'] = phot_table2['aperture_sum_1']*0 - 1
        phot_table2['median'] = phot_table2['aperture_sum_1']*0 - 1
        phot_table2['residual_aperture_sum_mean'] = phot_table2['aperture_sum_1']*0 - 1
        phot_table2['residual_aperture_sum_med'] = phot_table2['aperture_sum_1']*0 - 1

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


            else:
                print("[FluxFinder] Source {} is not within boundary {},{}; {},{}"
                        .format(i, x_max, y_max, x, y))

        phot_table2['residual_aperture_sum_mean'].info.format = '%.8g'  # for consistent table output
        phot_table2['residual_aperture_sum_med'].info.format = '%.8g'  # for consistent table output
                
        # MEDIAN

        # Background sub using median
        #phot_table2['median'] = phot_table2['aperture_sum_1']*0 - 1
        #phot_table2['residual_aperture_sum_med'] = phot_table2['aperture_sum_1']*0 - 1


        ##for each source
        #for q in range(n_sources_phot):
        #    x, y = positions[q]
 
        #    xypos = (x, y)

        #    ## TODO: Magic number
        #    #check that largest aperture does not exceed the boundaries of the image
        #    if Utilities.is_within_boundaries(x, y, x_max, y_max, 15):
        #        annulus = CircularAnnulus(
        #                xypos,
        #                r_in=self.config.inner_radius, 
        #                r_out=self.config.outer_radius)
        #        ann_mask = annulus.to_mask(method='center')
        #        #print(ann_mask)
        #        weighted_data = ann_mask.multiply(image_data)
        #        print(weighted_data)
        #        phot_table2['median'][q] = np.median(weighted_data[weighted_data != 0])

        #        
        #        # Calc the median background in the second aperture ring
        #        bkg_med = phot_table2['median'][q]
        #    
        #        # Calc background level in each main aperture and subtract
        #        bkg_sum = bkg_med * apertures.area
        #        final_sum = phot_table2['aperture_sum_0'][q] - bkg_sum
        #    
        #    
        #        phot_table2['residual_aperture_sum_med'][q] = final_sum
        #    else:
        #        print("[FluxFinder] Source {} is not within boundary {},{}; {},{}"
        #                .format(q, x_max, y_max, x, y))
        #
        #phot_table2['residual_aperture_sum_med'].info.format = '%.8g'  # for consistent table output
        
        return phot_table2
        


    

    ## TODO: Closer look
    ## TODO: Breaks if no sets
    def make_light_curves(self):
        if self.config.has_sets == False:
            print("[FluxFinder] Error: Cannot make light curves, not implemented for has_sets=False")
            exit()
        
        times = [line.rstrip(self.config.line_ending) for line in open(self.config.time_path)]
        
        light_curves = []
          
        #iterate over each image
                
        for _, s, i in Utilities.loop_images(self.config):
                #print("Finding fluxes in image {}".format(str((set-1)*self.set_size + i)))
                print("[FluxFinder] Making light curve for image: set {:1}; image: {:03}".format(s, i))
                
                t = self.find_fluxes(s, i)
                                
                #iterate through each source in the table
                for j in range(len(t['id'])):
                    
                    #add table objects to the light curves array when it is
                    #empty
                    if i == 1 and s == 1:
                        light_curves.append(Table(names = ('time','counts')))
                        
                    #get time matching the image
                    #date_and_time = times[(set-1) * self.set_size * 2 + (2*i-2)]
                    date_and_time = times[(s-1)*self.config.set_size + i-1]
                    

                    ## TODO: use datetime here
                    #just interested in time part of date and time for now
                    time = date_and_time.split("T")[1]
                    
                    #split time into hours, minutes and seconds
                    time_split = time.split(":")
                    time_split[0] = int(time_split[0])
                    time_split[1] = int(time_split[1])
                    time_split[2] = int(time_split[2])
                    
                    #if the observation start time is empty, add the current
                    #time in hours minutes and seconds 
                    if len(self.obs_start_time) == 0:
                        self.obs_start_time = time_split
                    
                    #get difference between time of start and time of 
                    #observation of current image, in hours minutes and seconds
                    time_elapsed_hms = Utilities.diff(time_split, self.obs_start_time)
                    
                    #if the hours value is positive this observation was taken before midnight
                    #the time elapsed from the start would be negative as the hour value
                    #woudl reset to 0
                    if time_elapsed_hms[0] >= 0:
                        #get time elapsed in seconds
                        time_elapsed = time_elapsed_hms[0] * 3600 + time_elapsed_hms[1] * 60 + time_elapsed_hms[2]
                    else:
                        time_elapsed_hms_till_midnight = Utilities.diff([23, 59, 60], self.obs_start_time)
                        time_elapsed = time_elapsed_hms_till_midnight[0] * 3600 + time_elapsed_hms_till_midnight[1] * 60 + time_elapsed_hms_till_midnight[2]
                        time_elapsed += time_split[0] * 3600 + time_split[1] * 60 + time_split[0]
                        
                    
                    #if median has reasonable value
                    if t['median'][j] > 0: 
                        light_curves[j].add_row((time_elapsed, str(t['residual_aperture_sum_med'][j])))
                    else:
                        print("[FluxFinder] Rejecting value for source {}, median is {}"
                                .format(t['id'][j], t['median'][j]))
                    #light_curves[j].add_row((time_elapsed, str(t['residual_aperture_sum_med'][j])))

        #loop through all light curves, writing them out to a file
        for j in range(len(light_curves)):
            #build light curve path
            fname = self.config.source_format_str.format(self.catalogue['id'][j])
            path = os.path.join(self.config.light_curve_dir, fname)

            table = light_curves[j]
            table.write(path, format=self.config.table_format, overwrite=True)
            
    
    
    
    
    def map_id(self, id2, cat1, cat2, shifts, set_number):
        """
        Find ID of a star from catalogue 2 in catalogue 1
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
            
        
        
        
    def plot_light_curve(self, source_id=None, curve_path=None, adjusted=False, show=False):
        """
        Plot light curve of star with the given ID from catalogue

        Parameters
        ----------

        source_id: int
            ID of the source to plot light curve for.
            Plots average light curve if `None`

        curve_path: string, optional
            Path to save the image to

        adjusted: bool
            Has the lighth curve been divided by the average flux?

        """

        print("[FluxFinder] Plotting light curve for source {} (adjusted={})"
                .format(source_id, adjusted))

        
        if curve_path == None:
            fname = self.config.source_format_str.format(source_id)
            if adjusted:
                curve_path = os.path.join(self.config.adjusted_curve_dir, fname)
            else:
                curve_path = os.path.join(self.config.light_curve_dir, fname)

        table = Table.read(curve_path, format=self.config.table_format)
        
        times = table['time']
        fluxes = table['counts']
        

        mean = Utilities.mean(fluxes)
        
        for i in range(len(fluxes)):
            fluxes[i] = fluxes[i]/mean
        
        minimum = min(fluxes)
        maximum = max(fluxes)
    
        ## TODO: Multiple axes, seconds and minutes?
        #plt.figure(figsize = (12, 8))
        plt.scatter(times, fluxes, s=5, marker = 'x');
        plt.xlabel("Time [seconds]")
        plt.ylabel("Relative flux [counts/mean]")
        plt.ylim(0.9 * minimum, maximum * 1.1)
        
        if source_id == None:
            fname = "LC_{}_avg.jpg".format(self.config.image_prefix)
            plt.title("Light curve for average of sources (adjusted={})"
                .format(adjusted))
        elif source_id >= 0:
            fname = "LC_{}_{}{:04}.jpg".format(
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
        plt.close()


        
    def plot_all_light_curves(self, ids, adjusted=False, show=False):
        for _i, path, source_id in Utilities.loop_variables(self.config, ids, adjusted=adjusted):
            self.plot_light_curve(source_id=source_id, curve_path=path, adjusted=adjusted, show=show)


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



    def divide_by_average(self):
        """
        Divides all light curves by the average light curve to remove noise.
        Also removes global effects of the field that vary over time.

        Creates an 'adjusted' light curve.

        """
        print("[DEBUG] Calling `divide_by_average` in FluxFinder")
                    
        #for all files
        for path, source_id in Utilities.list_sources(self.config, adjusted=False):
            t = Table.read(path, format=self.config.table_format)
                                            
            this_fluxes = t['counts']
            this_times = t['time']
            
            for i in range(len(this_fluxes)):
                time = this_times[i]
                
                #divide each flux measurement by the average flux measurement
                #at its respective point in time
                for j in range(len(self.avg_fluxes)):
                    if time == self.times[j]:
                        this_fluxes[i] = this_fluxes[i] / self.avg_fluxes[j]
                        #print(this_fluxes[i])
            
            #export adjusted light curve
            light_curve = Table([this_times, this_fluxes], names = ('time','counts'))

            fname = self.config.source_format_str.format(source_id)
            out_path = os.path.join(self.config.adjusted_curve_dir, fname)

            light_curve.write(out_path, format=self.config.table_format, overwrite=True)

        



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
            #x_shift, y_shift = self.get_total_shift(n, set_number)
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


# =============================================================================
#     #find the fluxes of all stars ina given image
#     def find_fluxes(self, image_number, set):
# 
#         #get total shift from first image to this one
#         x_shift, y_shift = self.get_total_shift(image_number, set)
#         
#         #add the shift onto the positions of the stars in the first image
#         #to find their positions in this image
#         x = self.catalogue['xcentroid'] + x_shift
#         y = self.catalogue['ycentroid'] + y_shift
#         
#         positions = (x, y) 
#         
#         # Shape and size of aperture object
#         apertures = CircularAperture(positions, r=9) 
#         
#         #build path to image 
#         data_path = self.directory + self.config.working_directory + self.config.image_directory + self.config.reduced_prefix + self.image_names
#         if not self.has_sets:
#             data_path += "000" + str(image_number)
#         else:
#             data_path += "_" + str(set) + "_" + Utilities.format_index(image_number)
#         
#         data_path += self.config.fits_extension
#             
#         #read in image data
#         image_data = fits.getdata(data_path, ext = 0)
#                 
#         # Local background subtraction
# 
#         # Define size of background aperture
#         annulus_apertures = CircularAnnulus(positions, r_in=10, r_out=15)
#         apertures = CircularAperture(positions, r=9)
#         apers = [apertures, annulus_apertures]
# 
#         # find counts in each aperture and annulus
#         phot_table2 = aperture_photometry(image_data, apers)
# 
#         for col in phot_table2.colnames:
#             phot_table2[col].info.format = '%.8g'  # for consistent table output
#     
#         # MEAN
#         # Background sub using mean
#             
#         phot_table2['residual_aperture_sum_mean'] = phot_table2['aperture_sum_1']*0 - 1
#         phot_table2['mean'] = phot_table2['aperture_sum_1']*0 - 1
# 
# 
#         for i in range(len(phot_table2)):
#             
#             x = positions[0][i]
#             y = positions[1][i]
#             
#             
#             #check that largest aperture does not exceed the boundaries of the image
#             if Utilities.is_within_boundaries(x, y, len(image_data[0]), len(image_data), 15):
#             
#                 # Calc the mean background in the second aperture ring
#                 bkg_mean = phot_table2['aperture_sum_1'][i] / annulus_apertures.area()
#                 phot_table2['mean'][i] = bkg_mean
#                 
#                 # Calc background level in each main aperture and subtract
#                 bkg_sum = bkg_mean * apertures.area()
#                 final_sum = phot_table2['aperture_sum_0'][i] - bkg_sum
#                 phot_table2['residual_aperture_sum_mean'][i] = final_sum
#         
#         phot_table2['residual_aperture_sum_mean'].info.format = '%.8g'  # for consistent table output
#                 
#         # MEDIAN
# 
#         # Background sub using median
#         phot_table2['median'] = phot_table2['aperture_sum_1']*0 - 1
#         phot_table2['residual_aperture_sum_med'] = phot_table2['aperture_sum_1']*0 - 1
# 
# 
#         #for each source
#         for q in range(0,len(phot_table2)):
#             x = positions[0][q]
#             y = positions[1][q]   
#             xypos = (x, y)
#             
#             #check that largest aperture does not exceed the boundaries of the image
#             if Utilities.is_within_boundaries(x, y, len(image_data[0]), len(image_data), 15):
#                 annulus = CircularAnnulus(xypos, r_in=10, r_out=15)
#                 ann_mask = annulus.to_mask(method='center')[0]
#                 weighted_data = ann_mask.multiply(image_data)
#                 phot_table2['median'][q] = np.median(weighted_data[weighted_data != 0])
#                 
#                 # Calc the median background in the second aperture ring
#                 bkg_med = phot_table2['median'][q]
#             
#                 # Calc background level in each main aperture and subtract
#                 bkg_sum = bkg_med * apertures.area()
#                 final_sum = phot_table2['aperture_sum_0'][q] - bkg_sum
#             
#             
#                 phot_table2['residual_aperture_sum_med'][q] = final_sum
#         
#         phot_table2['residual_aperture_sum_med'].info.format = '%.8g'  # for consistent table output
#         
#         #build path for flux table output
#         out_path = self.flux_dir + self.config.flux_prefix + self.image_names
#         
#         if not self.has_sets:
#             out_path += "000" + str(image_number)
#         else:
#             out_path += "_" + str(set) + "_" + Utilities.format_index(image_number)
#             
#         out_path += self.config.standard_file_extension
#         phot_table2.write(out_path,  format = self.config.table_format, overwrite = True)
#         
#     def make_light_curves(self):
#         
#         #get time of observation for each image
#         times_file = self.directory + self.config.working_directory + self.config.time_file
#         
#         times = [line.rstrip('\n') for line in open(times_file)]
#         
#         light_curves = []
#           
#         #iterate over each flux table
#         for set in range(1, self.n_sets+1):
#             for i in range(1, self.set_size + 1):
#                 print(i)
#                 
#                 #build flux path
#                 file = self.flux_dir + self.config.flux_prefix + self.image_names
#                 if not self.has_sets:
#                     file += "000" + str(i) 
#                 else:
#                     file += "_" + str(set) + "_" + Utilities.format_index(i)
#                 file += self.config.standard_file_extension
#                 
#                 #read flux data
#                 t = Table.read(file, format = self.config.table_format)
#                 
#                 #iterate through each source in the table
#                 for j in range(len(t['id'])):
#                     
#                     #add table objects to the light curves array when it is
#                     #empty
#                     if i == 1 and set == 1:
#                         light_curves.append(Table(names = ('time','counts')))
#                         
#                     #get time matching the image
#                     date_and_time = times[(set-1) * self.set_size + (i-1)]
#                     
#                     #just interested in time part of date and time for now
#                     time = date_and_time.split("T")[1]
#                     
#                     #split time into hours, minutes and seconds
#                     time_split = time.split(":")
#                     time_split[0] = int(time_split[0])
#                     time_split[1] = int(time_split[1])
#                     time_split[2] = int(time_split[2])
#                     
#                     #if the observation start time is empty, add the current
#                     #time in hours minutes and seconds 
#                     if len(self.obs_start_time) == 0:
#                         self.obs_start_time = time_split
#                     
#                     #get difference between time of start and time of 
#                     #observation of current image, in hours minutes and seconds
#                     time_elapsed_hms = Utilities.diff(time_split, self.obs_start_time)
#                     
#                     #if the hours value is positive this observation was taken before midnight
#                     #the time elapsed from the start would be negative as the hour value
#                     #woudl reset to 0
#                     if time_elapsed_hms[0] >= 0:
#                         #get time elapsed in seconds
#                         time_elapsed = time_elapsed_hms[0] * 3600 + time_elapsed_hms[1] * 60 + time_elapsed_hms[2]
#                     else:
#                         time_elapsed_hms_till_midnight = Utilities.diff([23, 59, 60], self.obs_start_time)
#                         time_elapsed = time_elapsed_hms_till_midnight[0] * 3600 + time_elapsed_hms_till_midnight[1] * 60 + time_elapsed_hms_till_midnight[2]
#                         time_elapsed += time_split[0] * 3600 + time_split[1] * 60 + time_split[0]
#                         
#                     
#                     #if median has reasonable value
#                     if t['median'][j] > 0: 
#                         light_curves[j].add_row((time_elapsed, str(t['residual_aperture_sum_med'][j])))
#                         
# 
#         #loop through all light curves, writing them out to a file
#         for j in range(len(light_curves)):
#             #build light curve path
#             file = self.light_curve_dir + self.image_names + self.config.identifier + str(j+1) + self.config.standard_file_extension
#            
#             table = light_curves[j]
#             table.write(file, format = self.config.table_format, overwrite=True)
#             
# =============================================================================
