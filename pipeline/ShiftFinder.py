from astropy.table import Table
from astropy.io import fits
from photutils import centroid_2dg
import os
import Constants
import Utilities
import numpy as np

class ShiftFinder:
    
    #directory containing raw images
    workspace_dir = None
    
    #image name regex
    image_prefix = None
    
    #are the images stored in sets?
    has_sets = None
    
    #number of sets 
    n_sets = None
    
    #number of images in each set
    set_size = None
    
    def __init__(self, workspace_dir, image_prefix, has_sets, set_size, n_sets):
        self.workspace_dir = workspace_dir
        self.image_prefix = image_prefix
        self.has_sets = has_sets
        self.n_sets = n_sets
        self.set_size = set_size


    
    # TODO: Quicksort, magic numbers
    #finds the coordinates of the brightest star in the image
    def get_reference_coordinates(self):
        
        #reads in catalogue
        table_fname = "{}{}{}".format(
                Constants.catalogue_prefix, self.image_prefix, Constants.standard_file_extension)
        table_path = os.path.join(self.workspace_dir, table_fname)
        t = Table.read(table_path, format=Constants.table_format)

        xs = t["xcentroid"]
        ys = t["ycentroid"]
        
        fluxes = t["flux"]
        max = 0
        
        Utilities.quicksort([fluxes, xs, ys], True)

        # Magic numbers
        i = int(len(fluxes)*0.8)
        n = 10
        
        if i + n >= len(fluxes):
            n = (len(fluxes) - i)
        
        x = xs[i:i+10]
        y = ys[i:i+10]
        
# =============================================================================
#         #finds index of brightest star
#         #DO THIS ON MULTIPLE STARS
#         for i in range(len(fluxes)):
#             if float(fluxes[i]) > fluxes[max]:
#                 max = i
#             
#         x = xs[max]
#         y = ys[max]
# =============================================================================
        
        return x, y
    

    ## TODO: Fit warnings are probably in here
    def find_shift(self, previous_x, previous_y, image_path, size=10):

        # choose how big a box to draw around the star
        #size = 10  # 10 will make a 20x20 pixel box
        
        # open the next image and load it as an array
        image = fits.open(image_path)
        data_array = image[0].data
        
        # Select just the box around the star
        data_square = data_array[int(previous_y)-size:int(previous_y)+size, int(previous_x)-size:int(previous_x)+size] # note your x,y coords need to be an int
        
        #get coordinates of the centre of the star in the data square
        x, y = centroid_2dg(data_square)

        #calculate shift between previous position and newly calculated
        #positions        
        x_shift = ((x - size) - (previous_x - int(previous_x)))
        # also correcting for the shift due to using int
        y_shift = ((y - size) - (previous_y - int(previous_y)))
        
        return x_shift, y_shift
    
        


    ## TODO: Fix fit warnings
    def get_all_shifts(self):
        
        #build filepath to file in which shifts will be stored       
        shift_path = os.path.join(Constants.workspace_dir, Constants.shift_fname)
        
        #empty shift file if it exists
        if(os.path.exists(shift_path)):
            open(shift_path, "w").close()
        
        x_shifts = []
        y_shifts = []
        
        #get coordinates of 10 bright stars in image to use as a reference
        xs, ys = self.get_reference_coordinates()
        
        j = 0
        
        
        #iterate through each image in each set
        for s in range(1, self.n_sets + 1):
            for i in range(1, self.set_size+1):    

                print("Finding shifts in image: set {:1}; image {:03}".format(s, i))

                ## TODO: Make this a utils function
                #build image file path
                if self.has_sets:
                    format_str = "{}{}_{}_{}{}".format(Constants.reduced_prefix, 
                            self.image_prefix, "{:1}", "{:03}", Constants.fits_extension)
                else:
                    format_str = "{}{}_{}{}".format(Constants.reduced_prefix, 
                            self.image_prefix, "{:04}", Constants.fits_extension)

                if self.has_sets:
                    image_path = os.path.join(self.workspace_dir, 
                            Constants.image_subdir, format_str.format(s, i))
                else:
                    image_path = os.path.join(self.workspace_dir, 
                            Constants.image_subdir, format_str.format(i))

                
                avg_x = []
                avg_y = []
                
                for k in range(len(xs)):
                    #find shift between the x and y of the reference star in the 
                    #previous image and the current image
                    x_shift, y_shift = self.find_shift(xs[k], ys[k], image_path)
                    
                    #append the total shift so far to the arrays containing
                    #the shifts
                    avg_x.append(x_shift)
                    avg_y.append(y_shift)
                
                if s == 1 and i == 1:
                    prev_x_shift = 0
                    prev_y_shift = 0
                else:
                    prev_x_shift = x_shifts[j-1]
                    prev_y_shift = y_shifts[j-1]
                    
                
                med_x = np.median(avg_x)
                med_y = np.median(avg_y)
                for l in range(len(xs)):
                    #update previous x and y
                    xs[l] += med_x
                    ys[l] += med_y
                
                x_shifts.append(med_x+prev_x_shift)
                y_shifts.append(med_y+prev_y_shift)
                print(x_shifts[j], y_shifts[j])

                j+= 1
                
        
        #make table of x and y shifts for output
        table = Table([x_shifts, y_shifts], names = ('xshifts','yshifts'))
        
        #export shifts as table 
        table.write(shift_path, format = Constants.table_format, overwrite=True)
    
## TODO: Remove
#    #I believe this is redundant
#    def find_shift_between_all_catalogues(self, image_size):  
#        
#        previous_cat = Table.read(self.directory + Constants.working_directory + Constants.catalogue_prefix + self.image_names + "_1" + Constants.standard_file_extension, format = Constants.table_format)
#        x_shifts = []
#        y_shifts = []
#        
#        for i in range(2, self.n_sets+1):
#            cat = Table.read(self.directory + Constants.working_directory + Constants.catalogue_prefix + self.image_names + "_" + str(i) + Constants.standard_file_extension, format = Constants.table_format)
#            x_shift, y_shift = self.find_shift_between_catalogues(previous_cat, cat, image_size)
#            
#            x_shifts.append(x_shift)
#            y_shifts.append(y_shift)
#            
#            previous_cat = cat
#            
#        shift_file = self.directory + Constants.working_directory + Constants.catalogue_shift_file
#
#        table = Table([x_shifts, y_shifts], names = ('xshifts','yshifts'))
#        
#        table.write(shift_file, format = Constants.table_format, overwrite=True)
            
# =============================================================================
# #I believe this is redundant 
# def find_shift_between_catalogues(catalogue1, catalogue2):
# 
#     x1s = catalogue1['xcentroid']
#     y1s = catalogue1['ycentroid']
#     
#     x2s = catalogue2['xcentroid']
#     y2s = catalogue2['ycentroid']
#     
#     fluxes = catalogue1["flux"]
#     
#     max = 0
#     
#     for i in range(len(fluxes)):
#         if float(fluxes[i]) > fluxes[max] and x1s[i] > Constants.image_width - 200 and x1s[i] < Constants.image_width + 200 and y1s[i] > Constants.image_height - 200 and y1s[i] < Constants.image_height + 200:
#             max = i
#     
#     x = x1s[max]
#     y = y1s[max]
#     
#     distances = set()
#     
#     for i in range(len(x1s)):
#                 
#         if i != max:
#             shifts = [0, 0]
#         
#             shifts[0] = round(x1s[i] - x)
#             shifts[1] = round(y1s[i] - y)
#             
#             distances.add(str(shifts[0]) + " " +  str(shifts[1]))
#             
#     for i in range(len(x2s)):
#         matches = 0
#         for j in range(len(x2s)):
#             if i != j:
#                 shift = [round(x2s[j] - x2s[i]), round(y2s[j] - y2s[i])]
#                 
#                 for k in range(len(distances)):
#                     if distances[k][0] - 3 < shift[0] and distances[k][0] + 3 > shift[0] and distances[k][1] - 3 < shift[1] and distances[k][1] + 3 > shift[1] and not k in matched:
#                         matches += 1
#                 
#                 if j > 0.01 * len(x2s) and matches < 0.1 * j:
#                     #print(matches, j)
# 
#                     break
#         
#         print(matches)
#         if matches > 0.3*len(distances):
# 
#             return x2s[i] - x1s[max], y2s[i] - y1s[max]
#                 
# 
#     
# =============================================================================
    


        
        
    

