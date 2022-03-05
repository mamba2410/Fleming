from astropy.table import Table
from PIL import Image
from astropy.visualization import (ZScaleInterval, LinearStretch, ImageNormalize)

import os 
import numpy as np
import matplotlib.pyplot as plt

from . import Utilities, Config

class DataAnalyser:
    
    config = None

    n_sources = 0
    source_ids = None

    source_means = None
    source_medians = None
    source_stds = None

    adjusted_source_means = None
    adjusted_source_medians = None
    adjusted_source_stds = None

    n_variables = 0
    variable_scores = None
    variable_ids = None
    variable_mask = None

    results_table = None
    

    def __init__(self, config):
        """
        Does all the number-crunching for the pipeline.
        Further refining the data from the light curves and 
        transforming it.

        Parameters
        ----------

        config: Config
            Config object

        """

        self.config = config


        ## Set our internal number of sources to how many we have
        ## in the catalogue
        cat = Utilities.read_catalogue(self.config)

        self.source_ids = cat['id']
        self.n_sources = len(self.source_ids)

        
    def get_means_and_stds(self, source_ids=None, adjusted=False):
        """
        Plot the mean against standard deviation for the light curves
        of all sources with ids as given.
        If none are given, do it for all sources.

        Note: No sigma-clipping is done (yet).

        Parameters
        ----------

        source_ids: numpy array (int), optional
            Array of all ids to use.

        adjusted: bool, optional
            Use adjusted light curves?

        Returns
        -------

        source_means: numpy array (float64)
            Means of all sources given.

        source_stds: numpy array (float64)
            Standard deviations of all sources given.

        source_medians: numpy array (float64)
            Medians of all sources given.

        """
        
        ## Make empty arrays
        source_medians = np.empty(self.n_sources)
        source_means   = np.empty(self.n_sources)
        source_stds    = np.empty(self.n_sources)
        
        if source_ids == None:
            source_ids = self.source_ids

        ## TODO: loop better
        ## For each source
        for i, path, source_id in Utilities.loop_variables(self.config, source_ids, adjusted=adjusted):
            
            ## Read light curve data from file
            lc = np.genfromtxt(path, dtype=self.config.light_curve_dtype)

            ## Get number of measurements
            n_measures = len(lc['time'])
                            
            ## TODO: Magic number, fix this logic
            if n_measures > self.config.set_size*self.config.n_sets / 3:
                
                source_means[i]   = np.mean(lc['counts'])
                source_stds[i]    = np.std(lc['counts'])
                source_medians[i] = np.median(lc['counts'])
                                                   
                ### TODO: Value checking
                #value = std/mean
                #if value > 0 and value < 2: #and mean > 0.02 and mean < 80:
                #    self.stds.append(value)
                #    self.means.append(mean)
                #    self.id_map.append(int(source_id))

        #Utilities.quicksort([self.means, self.stds, self.id_map], True)

        if adjusted:
            self.adjusted_source_means   = source_means
            self.adjusted_source_stds    = source_stds
            self.adjusted_source_medians = source_medians
        else:
            self.source_means   = source_means
            self.source_stds    = source_stds
            self.source_medians = source_medians
        
        return source_means, source_stds, source_medians


    def plot_means_and_stds(self, adjusted=True):
        """
        Plot means against standard deviations of each source.
        Overplotted in red are the stars deemed variable.

        Parameters
        ----------

        adjusted: bool, optional
            Has light curve been divided by average?

        """

        #plt.figure(figsize=(16, 10))

        plt.scatter(
                self.source_means,
                self.source_stds,
                marker='.');
        plt.scatter(
                self.source_means[self.variable_mask],
                self.source_stds[self.variable_mask],
                marker='.',
                color = 'red');

        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("Mean [counts]")
        plt.ylabel("Standard deviation [counts]")
        plt.title("Mean against standard deviation for all sources in field {}"
                .format(self.config.image_prefix))
        
        #ensure plot y axis starts from 0
        #plt.gca().set_ylim(bottom=1e-9)
        #plt.xlim(self.means[len(self.means)-1] * 1.3, self.means[0] * 0.7)
        
        if adjusted:
            fname = "std_dev_mean_adjusted{}".format(self.config.plot_file_extension)
        else:
            fname = "std_dev_mean{}".format(self.config.plot_file_extension)

        path = os.path.join(self.config.output_dir, fname)

        plt.savefig(path)
        plt.close()
    



    ## TODO: Why using index range?
    ## TODO: Better way of determining score
            
    
    def get_source_ids(self):
        """
        Give all the IDs of each source known by the object.

        Returns
        -------

        self.source_ids: numpy array
            array of all source ids known

        """

        return self.source_ids
            
        
        
    def get_ids_for_avg(self):
        """
        Find IDs of stars to produce average light curve.
        The average light curve is subtracted from each printed light curve to reduce bias
        Eligable stars must have a value for all images and must not be variable

        Assumes already have variable candidates
        Always adjusted=False

        Returns
        -------

        avg_ids: numpy array
            IDs of sources suitable for finding the average light curve

        """

        print("[DEBUG] Calling `get_ids_for_avg` in DataAnalyser")
        

        ## Get ids which aren'y deemed variable
        flat_ids = np.setxor1d(self.source_ids, self.variable_ids)

        avg_ids = np.zeros(len(flat_ids))
    
        accepted_index= 0
        for _i, path, source_id in Utilities.loop_variables(self.config, flat_ids, adjusted=False):
            curve = np.genfromtxt(path).transpose()
            n_measures = len(curve[0])

            ## Find indices where counts are positive
            positive_indices = np.where(curve[1] > 0)[0]
            
            ## If all our counts are positive
            ## (This means we didn't discard it in the flux finding phase)
            if len(positive_indices) == n_measures:
                avg_ids[accepted_index] = source_id
                accepted_index += 1

        avg_ids = np.trim_zeros(avg_ids)
        return avg_ids
    


    ## Note: Callum changed this algorithm
    def make_avg_curve(self, ids):
        """
        Creates an average light curve from the given ID array.

        Divides each light curve by the median counts (to account for different star 
        brightnesses).
        Takes the mean value at each time and considers this the average curve.

        (Numpy does the loop addition/division for us)

        Parameters
        ----------

        ids: numpy array
            IDs of the sources to use for the average curve.

        Returns
        -------

        avg_curve: numpy array 2d
            Light curve of times and counts of the average of each light curve.

        """

        print("[DEBUG] Calling `make_avg_curve` in DataAnalyser")

        ## TODO: Make this config-able
        n_measures = self.config.n_sets*self.config.set_size
        total_counts = np.zeros(n_measures)
        total_vars = np.zeros(n_measures)

        for i, path, source_id in Utilities.loop_variables(self.config, ids, adjusted=False):

            curve = np.genfromtxt(path, dtype=self.config.light_curve_dtype).transpose()
            
            fluxes = curve['counts']
            errs = curve['counts_err']

            #print("Source id: {}; {}, {}".format(source_id, np.min(fluxes), np.max(fluxes)))

            ## Normalise each curve to ~1
            ## Median may be better if we have huge outliers/clouds?
            #source_mean = np.mean(fluxes)
            #fluxes /= source_mean
            source_median = np.median(fluxes)

            ## TODO: Temporary bugfix
            ## Sometimes light curves write all zeros
            ## so we end up with nan here
            if source_median > 0:
                fluxes /= source_median
                errs /= source_median
                total_counts += fluxes
                total_vars += errs*errs

            
        
        ## Create average value at each measurement
        mean_counts = total_counts/n_measures
        mean_err = np.sqrt(total_vars)/n_measures

        ## TODO: Proper sigma clipping
        delete_idx = np.where(mean_err>1)[0]
        mean_counts[delete_idx] = 0

        write_me = np.array([curve['time'], mean_counts, mean_err]).transpose()

        np.savetxt(self.config.avg_curve_path, write_me)

        return write_me


    def create_avg_curve(self):
        """
        Wrapper function user should call to create the average light curve.
        Prepares the object to be in the state it expects to be 
        when creating the average light curve.

        """

        mean, std, med = self.get_means_and_stds(source_ids=None, adjusted=False)
        sorted_ids = self.sort_brightness()
        
        vd = VariableDetector(self.config, self.n_sources, mean, std, med, sorted_ids)
        exclude_ids = vd.std_dev_search(self.config.avg_exclude_threshold, adjusted=False)

        avg_ids = self.get_ids_for_avg(exclude_ids)
        self.make_avg_curve(avg_ids)


        
    def output_results(self):
        """
        Save the table of variable stars, in order of decreasing variability

        """
        
        #a = [self.results_table['variability'], self.results_table['id'], self.results_table['xcentroid'], self.results_table['ycentroid'], self.results_table['RA'], self.results_table['DEC']]
        
        cat = Utilities.read_catalogue(self.config)
        #indices = np.where(cat['id'] == self.variable_ids)
        _intersect, indices, _indices2 = np.intersect1d(cat['id'], self.variable_ids,
                return_indices=True, assume_unique=True)

        t = [('id', 'int64'),
            ('variability', 'float64'),
            ('xcentroid', 'float64'),
            ('ycentroid', 'float64'),
            ('RA', 'float64'),
            ('DEC', 'float64')]

        n_variables = len(self.variable_ids)

        ## Build the results table
        ## Wanted to do it better with numpy but it's hard to name 
        ## the columns
        results = np.empty(n_variables, dtype=t)
        for j, idx in enumerate(indices):
            results['id'][j] = cat['id'][idx]
            results['variability'][j] = self.variable_scores[idx]
            results['xcentroid'][j] = cat['xcentroid'][idx]
            results['ycentroid'][j] = cat['ycentroid'][idx]
            results['RA'][j] = cat['RA'][idx]
            results['DEC'][j] = cat['DEC'][idx]


        results_fname = "{}_results{}".format(self.config.image_prefix, self.config.standard_file_extension)
        results_path = os.path.join(self.config.output_dir, results_fname)

        #np.savetxt(results_path, results)
        np.savetxt(results_path, results, fmt='%04d %.10f %.10f %.10f %.10f %.10f')

        #Utilities.make_reg_file(self.config.output_dir,
        #        self.config.image_prefix + "_variables", self.results_table)

        self.results_table = results


       
    def remove_cosmics(self, lc, std):
        """
        Remove cosmic rays in light curve
        I don't know the algorithm but it seems to work

        Parameters
        ----------

        lc: numpy array 2d
            Light curve of source
            
        """
        
        counts = lc['counts']
        cosmic_index = -1
        
        n_measures = len(counts)
        
        for i in range(n_measures):
            
            m = counts[i]
            
            if i == 0:
                l = counts[i+1]
            else:
                l = counts[i-1]
                
            if i == len(counts) - 1:
                r = counts[i-1]
            else:
                r = counts[i+1]
        
            if m - r > self.config.cosmic_threshold * std and m - l > self.config.cosmic_threshold * std:
                
                if cosmic_index != -1:
                    return False
              
                cosmic_index = i
            
        if cosmic_index == -1:
            return False
        
        if i == 0:
            replacement = counts[1]
        else:
            replacement = counts[i-1]
                
        lc['counts'][cosmic_index] = replacement
            
        return True


          
    ## TODO: Magic numbers
    def create_thumbnails(self, ff, adjusted=False):
        """
        Creates images of the brightest and dimmest frames of each source deemed variable

        Parameters
        ----------

        ff: FluxFinder
            FluxFinder object used to get a thumbnail slice

            
        """

        print("[DataAnalyser] Making thumbnails")
        
        #for i, (path, source_id) in enumerate(Utilities.list_sources(self.config, adjusted=adjusted)):
        for i, path, source_id in Utilities.loop_variables(self.config, self.results_table['id']):
            curve = np.genfromtxt(path, dtype=self.config.light_curve_dtype).transpose()
            
            c = curve['counts']

            i_dim = np.argmin(c)
            i_bright = np.argmax(c)
            
            i_x = self.results_table['xcentroid'][i]
            i_y = self.results_table['ycentroid'][i]
            
            print("[DataAnalyser] Creating thumbnail for source id {:04}, centroid {},{}"
                    .format(source_id, i_x, i_y))
            
            ## Magic numbers
            dim = ff.get_thumbnail(i_dim+1, i_x, i_y, 20, True)
            bright = ff.get_thumbnail(i_bright+1, i_x, i_y, 20, True)

        
            fig = plt.figure()
            fig.add_subplot(1, 2, 1)
            plt.axis('off')

            dim_norm = ImageNormalize(dim, interval=ZScaleInterval(), stretch=LinearStretch())
            plt.imshow(dim, origin='upper', cmap='gray', norm = dim_norm)

            
            fig.add_subplot(1, 2, 2)
            plt.axis('off')

            bright_norm = ImageNormalize(bright, interval=ZScaleInterval(), stretch=LinearStretch())
            plt.imshow(bright, origin='upper', cmap='gray', norm = bright_norm)
            
            fname= "thumb_{}_{}{:04}{}".format(self.config.image_prefix, self.config.identifier,
                    source_id, self.config.plot_file_extension)
            path = os.path.join(self.config.output_dir, fname)

            plt.savefig(path)
            plt.close()

