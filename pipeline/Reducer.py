from __future__ import print_function
from astropy.io.fits import getdata
from astropy.io.fits import getheader
from astropy.io.fits import HDUList
from astropy.io.fits import PrimaryHDU
from astropy.io.fits import getval
import Utilities
import numpy as np
import PositionFinder

import os
import Constants
    
class Reducer:
    
    config = None
    
    #filter to use
    image_filter = None
        
    #bias file
    bias_frames = []
    
    #flatfield file
    flatfield_frames = []
    
    master_bias = None
    
    master_flat = None
    
    flat_median = None
    
    
    def __init__(self, config, image_filter):
        """
        Takes Config object and image filter, creates a Reducer object
        Also gets bias and flat frames

        """

        self.config = config
        self.image_filter = image_filter
        self.get_bias_and_flatfield()
        
        
    #get the bias and flatfield data files 
    def get_bias_and_flatfield(self):
        """
        Finds bias and flats in image_dir and populates an internal 
        array with the image data

        """

        for file in os.listdir(self.config.image_dir):
            
            if file[:len(self.config.bias_prefix)] == self.config.bias_prefix:
                self.bias_frames.append(getdata(
                    os.path.join(self.config.image_dir, file)))

            if file[:len(self.config.flat_prefix)] == self.config.flat_prefix:
                self.flatfield_frames.append(getdata(
                    os.path.join(self.config.image_dir, file)))
      
        print("[Reducer] Found {} bias frames".format(len(self.bias_frames)))
        print("[Reducer] Found {} flatfield frames".format(len(self.flatfield_frames)))
        
    def create_master_bias(self):
        """
        Creates a master bias from median of each pixel in bias frames
        """

        self.master_bias = np.median(self.bias_frames)
        print("[Reducer] Created master bias")
    
    def create_master_flat(self):
        """
        Creates a master flat from median of each pixel in flat frames.
        Also finds the median pixel value across all flat frames.
        """
        self.master_flat = np.median(self.flatfield_frames)
        #get median of flatfield
        self.flat_median = np.median(self.master_flat)
                
        print("[Reducer] Created master flatfield")
        
            

    def reduce(self, skip_existing=False):
        """
        Subtract master bias and divide by master flat for all images in 
        original image directory

        """

        self.create_master_bias()
        self.create_master_flat()
        
        #loop througheach file in directory 
        for file in os.listdir(self.config.image_dir):
            #if the first strlen characters are the image name regex then
            #this file is an image to be processed
            #if file[:strlen]==self.image_names and Constants.fits_extension in file:
            if file[:len(self.config.image_prefix)] == self.config.image_prefix:

                #if raw images are stored in sets
                if(self.config.has_sets):
                    
                    #if has sets, then files are stored with the following 
                    #suffix format 'name_1_001', 'name_3_020' etc. The 
                    #following code splits the name string to find the
                    #set number and image number
                    
                    array = file.split("_")
                    set_number   = int(array[len(array)-2])
                    image_number = int(array[len(array)-1].split(".")[0])
                    
                    ## Throws warning if we have a set size larger than specified in config.
                    if image_number > self.config.set_size:
                        print("[Reducer] Warning: Config has set size {}, but found image number {}"
                                .format(self.config.set_size, image_number))

                    
                    if set_number > self.config.n_sets:
                        print("[Reducer] Warning: Config has n_sets {}, but found set number {}"
                                .format(self.config.n_sets, set_number))
                                            
                    fname = self.config.image_format_str.format(set_number, image_number)
                    reduced_path = os.path.join(self.config.image_dir, fname)
                
                    print("[Reducer] Creating reduced image '{}'".format(reduced_path))
                else:

                    ## If no sets, then assumed file is named
                    ## <image_prefix>_<image_number>.<extension>
                    tmp = file.split("_")[1]
                    tmp = file.split(".")[0]
                    image_numer = int(tmp)

                    fname = self.config.image_format_str.format(image_number)
                    reduced_path = os.path.join(self.config.image_dir, fname)

                ## If reduced file already exists and skip_existing=True then skip this 
                if skip_existing and os.exists(reduced_path):
                    continue

                path = os.path.join(self.config.image_dir, file)
                #get image filter value
                image_filter=getval(path, 'FILTER',ignore_missing_end=True)

                #if the filter of the image matches the required filter then
                if image_filter[:1]==self.image_filter or image_filter == self.image_filter:
                    #get image data and header
                    data = getdata(path, ignore_missing_end=True) 
                    
                    image_height = len(data)
                    image_width  = len(data[0])
                    if image_height != self.config.image_height or image_width != self.config.image_width:
                        print("[Reducer] Error: Image dimensions different to what is given in config:")
                        print("......... (w,h): Config: ({},{}); image: ({},{})"
                                .format(self.config.image_width, self.config.image_height,
                                    image_width, image_height))
                    

                    #subtract bias from image
                    data = data - self.master_bias
                    
                    #divide image by flatfield divided by median of flatfield
                    data = data / (self.master_flat / self.flat_median)
                    
                    head=getheader(path, ignore_missing_end=True)
                    
                    
                        
                    
                    ## ?? Not sure why this is here
                    ## TODO: Magic numbers
                    if False and self.config.has_sets:
                        
                        x = int((image_width*0.1) + ((set_number*self.config.set_size + image_number)/400)*(0.7*width))
                        y = int(image_height/2)
                        
                        print(x, y, set_number, image_number)
                        for k in range(-3, 3, 1):
                            for l in range(-3, 3, 1):
                                
                                n = (k**2 + l**2)**0.5
                                if n == 0:
                                    n = 1
                                    
                                data[x + k][y + l] = 4.2/n * 8000
                        
                    #export processed image to file 
                    hdu = PrimaryHDU(data, head)
                    hdul = HDUList([hdu], None)
                    hdul.writeto(reduced_path, overwrite=True)

                    
    

