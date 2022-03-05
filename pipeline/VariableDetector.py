from . import Utilities, Config, PeriodFinder

import numpy as np

class VariableDetector:

    config = None

    n_sources = 0
    source_ids = None

    means = None
    stds  = None
    medians = None

    def __init__(self, config, source_ids, means, stds, medians, n_positive, adjusted=True):
        self.config = config

        self.n_sources = len(source_ids)
        self.source_ids = source_ids

        self.means = means
        self.stds = stds
        self.medians = medians
        self.n_positive = n_positive

        self.adjusted = adjusted


    def get_variable_score(self, source_id):
        """
        Returns a score determining the variability of the star.

        Standard deviation search.
        
        Score is the standard deviation of the source compared
        to the median of the surrounding stars (in index-space)

        Parameters
        ----------

        index: int
            Index of variable star to check?

        adjusted: bool, optional
            Has the source light curve been adjusted?

        """
        
        source_index = np.where(self.source_ids == source_id)[0][0]
       
        llim = source_index - self.config.check_radius
        ulim = source_index + self.config.check_radius + 1
        
        if llim < 0:
            llim = 0
        
        if ulim > self.n_sources:
            ulim = self.n_sources


        ## Median of standard deviations
        median_std = np.median(self.source_ids[llim:ulim])
        variable_score = self.source_ids[source_index]/median_std - 1

        return variable_score
        
        
    ## TODO: move signal to noise out of here
    def std_dev_search(self, threshold):
        """
        Calculates the variability score of all stars in the catalogue
        then builds a list of source indexes which are deemed variable.

        Parameters
        ----------
        
        adjusted: bool, optional
            Has the source light curve been adjusted?


        Returns
        -------

        self.variable_ids: numpy array
            Numpy array of all ids which are deemed to be variable candidates.

        """

        self.variable_scores = np.zeros(self.n_sources)

        ## Loop over sources in the catalogue
        for i, path, source_id in Utilities.loop_variables( \
                self.config, self.source_ids, adjusted=self.adjusted):

            self.variable_scores[i] = self.get_variable_score(source_id)
            lc = np.genfromtxt(path, dtype=self.config.light_curve_dtype)

        signal_to_noise = self.medians/self.stds
        self.variable_mask = np.where(
                (self.variable_scores > self.config.variability_threshold)
                & (self.variable_scores < self.config.variability_max)
                & (signal_to_noise > self.config.min_signal_to_noise)
                # TODO: something with self.n_positives
                )[0]


        variable_ids = np.copy(self.source_ids[self.variable_mask])
        self.variable_ids_s = variable_ids

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

        print("[VariableDetector] std search: Found {} variables out of {} sources"
                .format(len(variable_ids), self.n_sources))

        return self.variable_ids_s


    def amplitude_search(self):
        """

        """

        pf = PeriodFinder.PeriodFinder(self.config)

        amplitude_score = np.zeros(self.n_sources)

        for i, path, source_id in Utilities.loop_variables( \
                self.config, self.source_ids, adjusted=self.adjusted):

            print("[VariableDetector] Finding period for source {}/{}".format(i, self.n_sources))
            ## TODO: Make period range config variable
            _P, _P_err, A, A_err, _phi, _offset = pf.period_search(
                    source_id, path, n_samples=self.config.n_sample_periods)
            amplitude_score[i] = A/A_err

        np.savetxt(self.config.workspace_dir + "/amplitude_test.txt", amplitude_score)
        variable_mask = np.where(amplitude_score > self.config.amplitude_score_threshold)[0]

        variable_ids = np.copy(self.source_ids[variable_mask])
        self.variable_ids_a = variable_ids

        print("[VariableDetector] amplitude search: Found {} variables out of {} sources"
                .format(len(variable_ids), self.n_sources))

        return variable_ids


    def get_variables(self):
        """
        Function the user should call to get final IDs of sources deemed variable
        candiates.

        Returns
        -------

        variable_ids: numpy array
            IDs of each source considered variable

        """
        ## Get variables judged by different methods
        variable_ids_s = self.std_dev_search(self.config.variability_threshold)
        variable_ids_a = self.amplitude_search()

        ## Stars which apppear in both are returned
        #variable_ids = np.intersect1d(variable_ids_s, variable_ids_a)
        variable_ids = np.concatenate(variable_ids_s, variable_ids_a)

        return variable_ids


