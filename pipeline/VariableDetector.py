from . import PeriodFinder, Utilities, Config

import numpy as np

class VariableDetector:
    """
    VariableDetector class.

    Does everything to do with determining variability in a source.

    Usually associated with a DataAnalyser object.

    Attributes
    ----------


    config: Config
        Config object for the field

    n_sources: int
        Number of sources the object knows about

    source_ids: numpy array
        IDs of all sources we know about.

    means: numpy array
        Means of the light curves of each source.

    medians: numpy array
        Medians of the light curves of each source.

    stds: numpy array
        Standard deviations of the light curves of each source.

    n_positive: numpy array
        Number of positive values in each light curve

    adjusted: bool
        Are we an object for adjusted or non-adjusted set?

    variable_scores: numpy array
        Scores for standard deviation search
    
    amplitude_scores: numpy array
        Scores for period search

    period_stats: numpy array
        Table of primary period, amplitude, phase, offset and errors for each source.

    variable_ids_s: numpy array
        IDs of all sources deemed variable through the standard deviation check

    variable_ids_a: numpy array
        IDs of all sources deemed variable through the amplitude search

    """

    config = None

    n_sources = 0
    source_ids = None

    means = None
    stds  = None
    medians = None
    n_positive = None

    variable_scores = None
    amplitude_scores = None
    period_stats = None

    variable_ids_s = None
    variable_ids_a = None

    adjusted = False

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

        source_id: int
            ID of the source to get the score of

        Returns
        -------

        variable_score: int
            Score of the variable

        """
        
        source_index = np.where(self.source_ids == source_id)[0][0]
       
        llim = source_index - self.config.check_radius
        ulim = source_index + self.config.check_radius + 1
        
        if llim < 0:
            llim = 0
        
        if ulim > self.n_sources:
            ulim = self.n_sources


        ## Median of standard deviations
        median_std = np.median(self.stds[llim:ulim])
        variable_score = self.stds[source_index]/median_std - 1

        return variable_score
        
        
    ## TODO: move signal to noise calculation out of here
    ## TODO: Make quality checks optional (for average )
    def std_dev_search(self, threshold, threshold_upper, min_snr):
        """
        Calculates the variability score of all stars in the catalogue
        then builds a list of source indexes which are deemed variable.

        Parameters
        ----------
        
        threshold: float64
            Threshold score, above which is deemed variable

        Returns
        -------

        self.variable_ids_s: numpy array
            Numpy array of all ids which are deemed to be variable candidates.

        """

        self.variable_scores = np.zeros(self.n_sources)

        ## Loop over sources in the catalogue
        for i, path, source_id in Utilities.loop_variables( \
                self.config, self.source_ids, adjusted=self.adjusted):

            self.variable_scores[i] = self.get_variable_score(source_id)
            #lc = np.genfromtxt(path, dtype=self.config.light_curve_dtype)

        signal_to_noise = self.medians/self.stds
        self.variable_mask = np.where(
                (self.variable_scores > threshold)
                & (self.variable_scores < threshold_upper)
                & (signal_to_noise > min_snr)
                # TODO: something with self.n_positives
                )[0]


        variable_ids = np.copy(self.source_ids[self.variable_mask])
        self.variable_ids_s = variable_ids
        

        print("[VariableDetector] std search: Found {} variables out of {} sources"
                .format(len(variable_ids), self.n_sources))

        return variable_ids


    def amplitude_search(self, amplitude_score_threshold):
        """
        Finds the most prominent period for each source.
        Reports the period, amplitude, phase, offset and uncertainties.

        Parameters
        ----------

        amplitude_score_threshold: float64
            Threshold score, above which is deemed variable

        Returns
        -------

        period_stats: numpy array
            Table of the statistics for each source

        """

        pf = PeriodFinder(self.config)

        amplitude_score = np.zeros(self.n_sources)

        period_stats = np.zeros(self.n_sources, dtype=[
            ('period', 'float64'),
            ('period_err', 'float64'),
            ('amplitude', 'float64'),
            ('amplitude_err', 'float64'),
            ('phi', 'float64'),
            ('offset', 'float64'),
            ])

        for i, path, source_id in Utilities.loop_variables( \
                self.config, self.source_ids, adjusted=self.adjusted):

            print("[VariableDetector] Finding period for source {}/{}".format(i, self.n_sources))

            lc = np.genfromtxt(path, dtype=self.config.light_curve_dtype).transpose()

            ## TODO: Make period range config variable
            #period_stats[i] = pf.period_search(
            #        source_id, path, n_samples=self.config.n_sample_periods)
            period_stats[i] = pf.period_search_curve(source_id, 
                    lc['time'], lc['counts'], lc['counts_err'], n_samples=self.config.n_sample_periods)
            A = period_stats[i]['amplitude']
            A_err = period_stats[i]['amplitude_err']

            ## TODO: Use signal to noise or amplitude per std?
            #amplitude_score[i] = A/A_err
            #amplitude_score[i] = A/self.stds[i]
            if period_stats['period'][i] > 0:
                main_period = A*np.sin(2*np.pi/period_stats['period'][i] * lc['time'] + period_stats['phi'][i]) + period_stats['offset'][i]
                subtracted_curve = lc['counts'] - main_period
                new_std = np.std(subtracted_curve)
                amplitude_score[i] = self.stds[i]/new_std
            else:
                amplitude_score[i] = 0

        self.period_stats = period_stats
        self.amplitude_scores = amplitude_score

        variable_mask = np.where(amplitude_score > amplitude_score_threshold)[0]
        variable_ids = np.copy(self.source_ids[variable_mask])
        self.variable_ids_a = variable_ids

        print("[VariableDetector] amplitude search: Found {} variables out of {} sources"
                .format(len(variable_ids), self.n_sources))

        return variable_ids


    def get_variables(self):
        """
        Do all searches to find sources which are variable.
        Recommended to call each search yourself, this is just
        here for convenience.

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

        ## Stars which were in one or the other are returned
        ## TODO: May contain duplicate IDs
        variable_ids = np.concatenate(variable_ids_s, variable_ids_a)

        return variable_ids

    def get_period_stats(self, ids=[]):
        """
        Retrieve the period stats for the given IDs

        Parameters
        ----------

        ids: array, optional
            IDs to get stats of. If empty, gives all known stats

        Returns
        -------

        period_stats: numpy array
            Table of the statistics for each source

        """
        if len(ids) == 0:
            return self.period_stats
        else:
            _intersect, indices, _indices2 = np.intersect1d(self.source_ids, ids,
                return_indices=True, assume_unique=True)
            return self.period_stats[indices]


    def get_scores(self, ids=[]):
        """
        Retrieve the scores of each given source.

        Parameters
        ----------

        ids: array, optional
            IDs to get stats of. If empty, gives all known scores

        Returns
        -------

        period_stats: numpy array
            Table of the statistics for each source

        """

        if len(ids) == 0:
            return self.variable_scores, self.amplitude_scores
        else:
            _intersect, indices, _indices2 = np.intersect1d(self.source_ids, ids,
                return_indices=True, assume_unique=True)
            return self.variable_scores[indices], self.amplitude_scores[indices]


