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

    n_variables = 0
    variable_scores = None
    variable_ids = None
    variable_mask = None
    

    def __init__(self, config, source_ids):
        self.config = config
        self.source_ids = source_ids
        self.n_sources = len(source_ids)

        
    #plot the means and standard deviations of all light curves generated
    def get_means_and_stds(self, adjusted=False):
        """
        Plot the means and standard deviations of all light curves generated

        Loops over all sources, removes cosmics if adjusted.

        Parameters
        ----------

        adjusted: bool, optional
            use adjusted light curve?
        """
        
        source_medians = np.empty(n_sources)
        source_means = np.empty(n_sources)
        source_stds  = np.empty(n_sources)
        

        ## TODO: loop better
        #for each file in the light curve directory 
        for source_index, path, source_id in Utilities.loop_variables(self.config, self.source_ids, adjusted=adjusted):
            
            #read light curve data from file
            lc = np.genfromtxt(path, dtype=[
                ('time', 'float64'),
                ('counts', 'float64')
                ])
            n_measures = len(lc['time'])
            
            ## TODO: Move this
            if adjusted:
                if self.remove_cosmics(t):
                    print("[DataAnalyser] Removed cosmics on id {}", source_id)
                
                            
            ## TODO: Magic number
            #only plot data point if at least 5 non-zero counts are recorded
            if n_measures > self.config.set_size*self.config.n_sets / 3:
                
                self.source_means[i]   = np.mean(lc['counts'])
                self.source_stds[i]    = np.std(lc['counts'])
                self.source_medians[i] = np.median(lc['counts'])
                                                   
                ## TODO: Value checking
                if value > 0 and value < 2: #and mean > 0.02 and mean < 80:
                    self.stds.append(value)
                    self.means.append(mean)
                    self.id_map.append(int(source_id))

        #Utilities.quicksort([self.means, self.stds, self.id_map], True)
        


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
    def get_variable_score(self, source_index):
        """
        Returns a score determining the variability of the star

        the next part is just the same as the is_variable method
        in Cataloguer - need to investigate if there is a reason for this
        Merge both together?

        Parameters
        ----------

        index: int
            Index of variable star to check?

        """
        
       
        llim = source_index - self.config.check_radius
        ulim = source_index + self.config.check_radius + 1
        
        if llim < 0:
            llim = 0
        
        if ulim > len(self.means):
            ulim = len(self.means)
        

        ## Median of standard deviations
        ## TODO: Why index range instead of flux range?
        ## TODO: Move this out of here
        median_std = np.median(self.source_stds[llim:ulim])
        #std_threshold = median_std * (1 + self.config.variability_threshold)

        #if self.source_stds[index] > std_threshold:
        #    self.variable_mask = np.append(self.variable_mask, source_index)
        #    #self.variable_means.append(self.means[index])
        #    #self.variable_stds.append(self.stds[index])

        variable_score = (self.source_stds[source_index] - median_std)/ median_std
        return variable_score
        
        

    def get_variable_ids(self):
        """
        Calculates the variability score of all stars in the catalogue
        then builds a list of source indexes which are deemed variable

        """

        #cat = np.genfromtxt(self.config.catalogue_path, dtype=[
        #    ('id', 'int64'),
        #    ('xcentroid', 'float64'),
        #    ('ycentroid', 'float64'),
        #    ('sharpness', 'float64'),
        #    ('roundness1', 'float64'),
        #    ('roundness2', 'float64'),
        #    ('npix', 'float64'),
        #    ('sky', 'float64'),
        #    ('peak', 'float64'),
        #    ('flux', 'float64'),
        #    ('mag', 'float64'),
        #    ])
        #self.results_table = Table(names = ('id', 'xcentroid', 'ycentroid', 'variability', 'RA', 'DEC'))
        
        print("[DEBUG] Catalogue has {} sources, DA has {}"
                .format(len(cat['id']), n_sources))

        self.variable_scores = np.zeros(self.n_sources)

        ## Loop over sources in the catalogue
        for i in range(self.n_sources):
            self.variable_scores[i] = self.get_variable_score(i)
            
            if self.variable_scores[i] > self.config.variability_threshold:

                ## TODO: Do this after the loop and join with no.where
                self.variable_ids = np.append(self.variable_ids, self.source_ids[i])
                self.variable_mask = np.append(self.variable_mask, i)
                
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
            
            

            
        
        
    ## TODO: Magic numbers
    #same function in cataloguer - remove one
    #probably from cataloguer
    #comments in other
    def get_ids_for_avg(self):
        """
        Find IDs of stars with high brightness to produce average light curve.
        The average light curve is subtracted from each printed light curve to reduce bias

        """

        print("[DEBUG] Calling `get_ids_for_avg` in DataAnalyser")
        
        ids = []
    
    
        for i in range(len(self.means)-1, int((len(self.means)-1) * 0.9), -1):
            #if self.means[i] > 50 and self.stds[i] < 0.04:
            #if not Utilities.is_above_line(-0.0001, 0.03, self.means[i], self.stds[i], 0.01) and self.means[i] > 5:
            #if not Utilities.is_above_line(self.means[i], self.stds[i], 2.2222*10**-9, 0.05777778, 0.001) and self.means[i] > 10^6:
            if not self.id_map[i] in self.variable_ids:
                fname = self.config.source_format_str.format(self.id_map[i])
                light_curve_path = os.path.join(self.config.light_curve_dir, fname)

                t = Table.read(light_curve_path, format=self.config.table_format)
                
                if len(t['time']) == self.config.set_size * self.config.n_sets:
                    ids.append(self.id_map[i])
        return ids
    


    #duplicate method in ff - remove this?
    #comments in other
    def make_avg_curve(self, ids):
        print("[DEBUG] Calling `make_avg_curve` in DataAnalyser")

        for i in range(len(ids)):
            fname = self.config.source_format_str.format(self.id_map[i])
            light_curve_path = os.path.join(self.config.light_curve_dir, fname)
            
            t = Table.read(light_curve_path, format=self.config.table_format)
            
            fluxes = t['counts']
            
            if i == 0:
                for j in range(len(fluxes)):
                    self.avg_fluxes.append(0)
                self.times = t['time']

            
            mean = Utilities.mean(fluxes)
                        
            for k in range(len(fluxes)):
                self.avg_fluxes[k] = self.avg_fluxes[k] + (fluxes[k] / mean)
            
        for l in range(len(self.avg_fluxes)):
            self.avg_fluxes[l] = self.avg_fluxes[l] / len(ids)
        
        light_curve = Table([self.times, self.avg_fluxes], names = ('time','counts') )

        fname= "{}_avg{}".format(self.image_prefix, self.config.standard_file_extension)
        path = os.path.join(self.config.light_curve_dir, fname)
        light_curve.write(path, format=self.config.table_format, overwrite=True)



    ## TODO: 
    #duplicate method in ff - remove this?
    #comments in other
    def divide_by_average(self):
        print("[DEBUG] Calling `divide_by_average` in DataAnalyser")
                    
        #for all source files
        for path, source_id in Utilities.list_sources(self.config):
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


       
    def remove_cosmics(self, t):
        """
        Remove cosmic rays in image

        Parameters
        ----------

        t: Table
            Table contaning source information. Counts will be modified
            
        """
        
        counts = t['counts']
        cosmic_index = -1
        
        if len(counts) == 0:
            return
        
        std = Utilities.standard_deviation(counts)

        
        for i in range(len(counts)):
            
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
            
        print("[DataAnalyser] cosmic detected at approx {}s in id:"
            .format(35*cosmic_index))
        ## TODO: Change the 35 (image cadence)
        
        if i == 0:
            replacement = counts[1]
        else:
            replacement = counts[i-1]
                
        t['counts'][cosmic_index] = replacement
            
        return True
          

