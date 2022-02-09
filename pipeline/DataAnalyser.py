from astropy.table import Table
import Constants
import os 
import Utilities
import matplotlib.pyplot as plt
import FluxFinder
import numpy as np
from PIL import Image
from astropy.visualization import (ZScaleInterval, LinearStretch,
                                   ImageNormalize)

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
    

    def __init__(self, config):
        self.config = config

        cat = Utilities.read_catalogue(self.config)

        self.source_ids = cat['id']
        self.n_sources = len(self.source_ids)

        
    #plot the means and standard deviations of all light curves generated
    ## TODO: Globals are overwritten on adjusted/non-adjusted
    def get_means_and_stds(self, source_ids=None, adjusted=False):
        """
        Plot the means and standard deviations of all light curves generated

        Loops over all sources.

        Parameters
        ----------

        adjusted: bool, optional
            use adjusted light curve?
        """
        
        source_medians = np.empty(self.n_sources)
        source_means   = np.empty(self.n_sources)
        source_stds    = np.empty(self.n_sources)
        
        if source_ids == None:
            source_ids = self.source_ids

        ## TODO: loop better
        #for each file in the light curve directory 
        for i, path, source_id in Utilities.loop_variables(self.config, source_ids, adjusted=adjusted):
            
            #read light curve data from file
            lc = np.genfromtxt(path, dtype=[
                ('time', 'float64'),
                ('counts', 'float64')
                ])
            n_measures = len(lc['time'])
            
                
                            
            ## TODO: Magic number
            #only plot data point if at least 5 non-zero counts are recorded
            if n_measures > self.config.set_size*self.config.n_sets / 3:
                
                source_means[i]   = np.mean(lc['counts'])
                source_stds[i]    = np.std(lc['counts'])
                source_medians[i] = np.median(lc['counts'])
                                                   
                ### TODO: Value checking
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


    def plot_means_and_stds(self, adjusted=False):
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
    



    ## TODO: Merge with `is_variable` in Cataloguer
    ## TODO: Why using index range?
    def get_variable_score(self, source_index, adjusted=False):
        """
        Returns a score determining the variability of the star.
        
        Score is the standard deviation of the source compared
        to the median of the surrounding stars (in index-space)

        Parameters
        ----------

        index: int
            Index of variable star to check?

        """
        
       
        llim = source_index - self.config.check_radius
        ulim = source_index + self.config.check_radius + 1
        
        if llim < 0:
            llim = 0
        
        if ulim > self.n_sources:
            ulim = self.n_sources
        
        ## Debug
        #llim = 0
        #ulim = self.n_sources

        ## Median of standard deviations
        ## TODO: Why index range instead of flux range?
        ## TODO: Move this out of here
        if adjusted:
            median_std = np.median(self.adjusted_source_stds[llim:ulim])
            variable_score = self.adjusted_source_stds[source_index]/median_std - 1
        else:
            median_std = np.median(self.source_stds[llim:ulim])
            variable_score = self.source_stds[source_index]/median_std - 1

        #std_threshold = median_std * (1 + self.config.variability_threshold)

        #if self.source_stds[index] > std_threshold:
        #    self.variable_mask = np.append(self.variable_mask, source_index)
        #    #self.variable_means.append(self.means[index])
        #    #self.variable_stds.append(self.stds[index])

        return variable_score
        
        
    def get_source_ids(self):
        return self.source_ids

    def get_variable_ids(self, adjusted=True):
        """
        Calculates the variability score of all stars in the catalogue
        then builds a list of source indexes which are deemed variable

        TODO: Bug of id0000 never deemed variable

        """

        cat = Utilities.read_catalogue(self.config)
        n_sources_cat = len(cat['id'])

        if n_sources_cat != self.n_sources:
            print("[DEBUG] Catalogue has {} sources, DA has {}"
                    .format(n_sources_cat, self.n_sources))


        self.variable_scores = np.zeros(self.n_sources)

        ## Loop over sources in the catalogue
        for i in range(self.n_sources):
            self.variable_scores[i] = self.get_variable_score(i, adjusted=adjusted)

        self.variable_mask = np.where(self.variable_scores > self.config.variability_threshold)[0]
        self.variable_ids = np.copy(self.source_ids[self.variable_mask])

            #else:
            #    print("[DEBUG] Rejecting, score of {}".format(self.variable_scores[i]))
                
                ## TODO: Check if catalogue is always ordered in ids
                #row_index = np.where(cat['id'].data==variable_id)

                ## TODO: Build table? Join with existing data
                #if 'RA' in cat.colnames:
                #    self.results_table.add_row([
                #        int(variable_id),
                #        cat['xcentroid'][row_index],
                #        cat['ycentroid'][row_index],
                #        variability,
                #        cat['RA'][row_index],
                #        cat['DEC'][row_index]])
                #else:
                #    self.results_table.add_row([
                #        int(variable_id),
                #        cat['xcentroid'][row_index],
                #        cat['ycentroid'][row_index], 
                #        variability,
                #        0,
                #        0])
        

        ## TODO: write results_table to file?

        print("[DataAnalyser] Found {} variables out of {} sources"
                .format(len(self.variable_ids), self.n_sources))

        return self.variable_ids
            
    
    ## TODO: Magic numbers
    def create_thumbnails(self, ff, adjusted=False):
        """
        Creates images of the brightest and dimmest frames of each source deemed variable

        Parameters
        ----------

        ff: FluxFinder
            FluxFinder object used to get a thumbnail slice
        """
        
        #for i, (path, source_id) in enumerate(Utilities.list_sources(self.config, adjusted=adjusted)):
        for i, path, source_id in Utilities.loop_variables(self.config, self.results_table['id']):
            t = Table.read(path, format=self.config.table_format)
            
            i_dim = 0
            i_bright = 0
            c = t['counts']
            
            for j in range(len(c)):
                if c[j] > c[i_bright]:
                    i_bright = j
                if c[j] < c[i_dim]:
                    i_dim = j
            
            i_x = self.results_table['xcentroid'][i]
            i_y = self.results_table['ycentroid'][i]
            
            print("[DataAnalyser] Creating thumbnail for source id {:04}, centroid {},{}".format(
                source_id, i_x, i_y))
            
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
            
            

            
        
        
    def get_ids_for_avg(self):
        """
        Find IDs of stars to produce average light curve.
        The average light curve is subtracted from each printed light curve to reduce bias
        Eligable stars must have a value for all images and must not be variable

        Assumes already have variable candidates
        Always adjusted=False

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

        (Numpy does the loop addition/division for us)

        """

        print("[DEBUG] Calling `make_avg_curve` in DataAnalyser")

        ## TODO: Make this config-able
        n_measures = self.config.n_sets*self.config.set_size
        total_counts = np.zeros(n_measures)

        for i, path, source_id in Utilities.loop_variables(self.config, ids, adjusted=False):

            curve = np.genfromtxt(path, dtype=[
                ('time', 'float64'),
                ('counts', 'float64'),
                ]).transpose()
            
            fluxes = curve['counts']

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
                total_counts += fluxes

            
        
        ## Create average value at each measurement
        mean_counts = total_counts/n_measures

        ## TODO: Proper sigma clipping
        delete_idx = np.where(mean_counts>10)[0]
        mean_counts[delete_idx] = 1

        write_me = np.array([curve['time'], mean_counts]).transpose()

        np.savetxt(self.config.avg_curve_path, write_me)

        return write_me


    def create_avg_curve(self):
        _ = self.get_means_and_stds(source_ids=None, adjusted=False)
        variable_ids = self.get_variable_ids(adjusted=False)
        avg_ids = self.get_ids_for_avg()
        self.make_avg_curve(avg_ids)


    def get_variables(self):
        _means, stds, _meds = self.get_means_and_stds(source_ids=None, adjusted=True)

        ## Remove any cosmic rays in the light curve
        for i, path, source_id in Utilities.loop_variables(self.config, self.source_ids, adjusted=True):
            curve = np.genfromtxt(path, dtype=[
                ('time', 'float64'),
                ('counts', 'float64'),
                ]).transpose()

            if self.remove_cosmics(curve, stds[i]):
                #print("[DataAnalyser] Removed cosmics on id {}".format(source_id))
                np.savetxt(path, curve)

        variable_ids = self.get_variable_ids(adjusted=True)

        return variable_ids


        
    ## TODO: Quicksort
    #save table of variable stars (in order of decreasing variability)
    def output_results(self):
        
        a = [self.results_table['variability'], self.results_table['id'], self.results_table['xcentroid'], self.results_table['ycentroid'], self.results_table['RA'], self.results_table['DEC']]
        
        
        Utilities.quicksort(a, False)

        results_fname = "{}_results{}".format(self.config.image_prefix, self.config.standard_file_extension)
        results_path = os.path.join(self.config.output_dir, results_fname)

        self.results_table.write(results_fname, format=self.config.table_format, overwrite=True)
        Utilities.make_reg_file(self.config.output_dir,
                self.config.image_prefix + "_variables", self.results_table)


       
    def remove_cosmics(self, lc, std):
        """
        Remove cosmic rays in light curve
        I don't know the algorithm but it seems to work

        Parameters
        ----------

        lc: numpy array
            
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
          

