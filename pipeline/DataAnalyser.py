from astropy.table import Table

import os 
import numpy as np
import matplotlib.pyplot as plt

from . import Utilities, Config, VariableDetector

class DataAnalyser:
    """
    DataAnalyser class.

    Does all the number-crunching for the pipeline.
    Further refining the data from the light curves and 
    transforming it.

    All arrays are sorted in order of descending brightness
    as per the catalogue.

    One per set of light curves (adjusted or not).

    Attributes
    ----------

    config: Config
        Config object for the field

    n_sources: int
        Number of sources the object knows about

    source_ids: numpy array
        IDs of all sources we know about.

    source_means: numpy array
        Means of the light curves of each source.

    source_medians: numpy array
        Medians of the light curves of each source.

    source_stds: numpy array
        Standard deviations of the light curves of each source.

    n_positive: numpy array
        Number of positive values in each light curve

    adjusted: bool
        Are we an object for adjusted or non-adjusted set?

    """
    
    config = None

    n_sources = 0
    source_ids = None

    source_means = None
    source_medians = None
    source_stds = None
    n_positive = None

    adjusted = False
    

    def __init__(self, config, adjusted=False):
        self.config = config
        self.adjusted = adjusted


        ## Set our internal number of sources to how many we have
        ## in the catalogue
        cat = Utilities.read_catalogue(self.config)

        ## Sort our sources in order of brightness
        brightness_indices = np.flip(np.argsort(cat['flux']))
        self.source_ids = cat['id'][brightness_indices]

        self.n_sources = len(self.source_ids)

        
    def get_means_and_stds(self):
        """
        Calculate the statistics of each light curve and store
        in the object.

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

        n_positive: numpy array
            Number of positive counts in the light curve

        """
        
        ## Make empty arrays
        source_medians = np.zeros(self.n_sources)
        source_means   = np.zeros(self.n_sources)
        source_stds    = np.zeros(self.n_sources)
        n_positives    = np.zeros(self.n_sources)

        ## For each source
        for i, path, source_id in Utilities.loop_variables(self.config, \
                self.source_ids, adjusted=self.adjusted):
            
            ## Read light curve data from file
            lc = np.genfromtxt(path, dtype=self.config.light_curve_dtype)

            ## Get number of measurements
            n_measures = len(lc['time'])
                            
            ## TODO: Magic number, fix this logic
            if n_measures > self.config.set_size*self.config.n_sets / 3:
                
                source_means[i]   = np.mean(lc['counts'])
                source_stds[i]    = np.std(lc['counts'])
                source_medians[i] = np.median(lc['counts'])
                n_positives[i]    = len(np.where(lc['counts'] > 0)[0])
                                                   
                ### TODO: Value checking, do somewhere else
                #value = std/mean
                #if value > 0 and value < 2: #and mean > 0.02 and mean < 80:
                #    self.stds.append(value)
                #    self.means.append(mean)
                #    self.id_map.append(int(source_id))

        self.source_means   = source_means
        self.source_stds    = source_stds
        self.source_medians = source_medians
        self.n_positives    = n_positives
        
        return source_means, source_stds, source_medians, n_positives



    ## TODO: Allow overplotting of variables
    def plot_means_and_stds(self):
        """
        Plot the mean against standard deviation for the light curves
        of all sources with ids as given.
        If none are given, do it for all sources.

        Parameters
        ----------

        adjusted: bool, optional
            Has light curve been divided by average?

        """

        plt.scatter(
                self.source_means,
                self.source_stds,
                marker='.');

        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("Mean [counts]")
        plt.ylabel("Standard deviation [counts]")
        plt.title("Mean against standard deviation for all sources in field {}"
                .format(self.config.image_prefix))
        

        if self.adjusted:
            fname = "std_dev_mean_adjusted{}".format(self.config.plot_file_extension)
        else:
            fname = "std_dev_mean{}".format(self.config.plot_file_extension)

        path = os.path.join(self.config.output_dir, fname)

        plt.savefig(path)
        plt.close()
    
    
    def get_source_ids(self):
        """
        Give all the IDs of each source known by the object.

        Returns
        -------

        self.source_ids: numpy array
            array of all source ids known

        """

        return self.source_ids
            
        
        
    def get_ids_for_avg(self, exclude_ids):
        """
        Find IDs of stars to produce average light curve.
        The average light curve is subtracted from each printed light curve to reduce bias
        Eligable stars must have a value for all images and must not be variable

        Assumes already have variable candidates
        Always adjusted=False

        Parameters
        ----------

        exclude_ids: numpy array
            IDs of all the sources which we should also avoid.
            Usually given by the prototype variable finder.

        Returns
        -------

        avg_ids: numpy array
            IDs of sources suitable for finding the average light curve

        """

        n_measures = self.config.set_size * self.config.n_sets
        avg_ids = np.zeros(self.n_sources)
    
        accepted_index= 0
        for i, _path, source_id in Utilities.loop_variables( \
                self.config, self.source_ids, adjusted=self.adjusted):
            
            ## If we run across a source we should avoid
            ## skip this loop iteration
            if source_id in exclude_ids:
                continue

            ## If all our counts are positive
            ## (This means we didn't discard it in the flux finding phase)
            if self.n_positives[i] == n_measures:
                avg_ids[accepted_index] = source_id
                accepted_index += 1

        ## Get rid of trailing zeros
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
            else:
                print("[DataAnalyser] Error: Encountered non-positive median when creating \
                        an average light curve")

            
        
        ## Create average value at each measurement
        mean_counts = total_counts/n_measures
        mean_err = np.sqrt(total_vars)/n_measures

        ## TODO: Error checking, this shouldn't trigger
        delete_idx = np.where(mean_err>1)[0]
        mean_counts[delete_idx] = 0

        write_me = np.array([curve['time'], mean_counts, mean_err]).transpose()

        np.savetxt(self.config.avg_curve_path, write_me)

        #return write_me

        
    ## TODO: Move me
    def output_results(self, variable_ids, vd, out_dir=None):
        """
        Save the table of variable stars.

        Parameters
        ----------

        variable_ids: numpy array
            IDs of all the sources deemed variable.

        vd: VariableDetector
            Variable detector object associated with this data analyser object

        Results
        -------

        results_table: numpy array
            Numpy array containing all details of the sources.
            Same object that is written to disk.

        """
        
        #a = [self.results_table['variability'], self.results_table['id'], self.results_table['xcentroid'], self.results_table['ycentroid'], self.results_table['RA'], self.results_table['DEC']]
        
        cat = Utilities.read_catalogue(self.config)

        ## Find out where the variables are in the catalogue
        _intersect, indices, _indices2 = np.intersect1d(cat['id'], variable_ids,
                return_indices=True, assume_unique=True)

        ## Columns of the results table
        t = [('id', 'int64'),
            ('variability_score', 'float64'),
            ('amplitude_score', 'float64'),
            ('flux', 'float64'),
            ('xcentroid', 'float64'),
            ('ycentroid', 'float64'),
            ('RA', 'float64'),
            ('DEC', 'float64')]

        n_variables = len(variable_ids)

        ## Retrieve the scores of each variable
        ## TODO: breaks if you don't do both searches. Only one should be valid enough
        variability_score, amplitude_score = vd.get_scores(variable_ids)

        ## Build the results table
        ## j is index of source in results table
        ## idx is index of source in catalogue
        results = np.empty(n_variables, dtype=t)
        for j, idx in enumerate(indices):
            results['id'][j] = cat['id'][idx]
            results['variability_score'][j] = variability_score[j]
            results['amplitude_score'][j] = amplitude_score[j]
            results['flux'][j] = cat['flux'][idx]
            results['xcentroid'][j] = cat['xcentroid'][idx]
            results['ycentroid'][j] = cat['ycentroid'][idx]
            results['RA'][j] = cat['RA'][idx]
            results['DEC'][j] = cat['DEC'][idx]


        results_fname = "{}_results{}".format(self.config.image_prefix, self.config.standard_file_extension)
        if out_dir == None:
            results_path = os.path.join(self.config.output_dir, results_fname)
        else:
            results_path = os.path.join(out_dir, results_fname)

        np.savetxt(results_path, results, fmt='%04d %.10f %.10f %.10f %.10f %.10f %.10f %.10f')

        #Utilities.make_reg_file(self.config.output_dir,
        #        self.config.image_prefix + "_variables", self.results_table)

        return results



