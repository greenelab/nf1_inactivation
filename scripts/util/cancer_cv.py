'''
(C) Gregory Way 2016
NF1 Inactivation Classifier for Glioblastoma
cancer_cv.py

Description:
Classes to perform cross validation on pancan datasets

Usage:
Import only
'''


class cancer_cv(object):
    """
    A class to determine cross validation intervals for pancancer tissues with
    careful consideration to split cv folds with equal number of mutations and
    tissues if applicable.
    """

    def __init__(self, y_dict, fold, hold, min_samples, seed, tissues='auto'):
        """
        :param y_dict: keys are tissue types, values are pandas DataFrame
                       holding mutation status for each sample
        :param fold: integer of how many cross validation folds to extract
        :param hold: boolean of whether or not to extract holdout set
        :param seed: set a random seed for reproducibility
        :param tissues: a list of tissues to use in analysis
                        (defaults to auto-detect tissues)
        :return: a dictionary holding y matrices for each fold
        """

        self.y_dictionary = y_dict
        self.num_folds = fold
        self.seed = seed
        self.holdout = hold
        self.percent_held = None
        self.cv_folds = None
        self.x_folds = None

        if tissues == 'auto':
            # Determine the appropriate tissues to use since not all tissues
            # will have enough mutations
            self.tissues = []
            for tissue in y_dict.keys():
                tissue_DF = y_dict[tissue]
                mut_samples = tissue_DF[tissue_DF['Status'] == 1]
                if mut_samples.shape[0] > min_samples:
                    self.tissues.append(tissue)
        else:
            self.tissues = tissues

    def holdout_samples(self, fraction_held):
        self.fraction_held = fraction_held

    def assign_folds(self):
        """
        Subsets X and Y matrices to implement cross validation logic.
        self.holdout_samples must be called prior if self.holdout = T

        Output:
        a dictionary holding evenly distributed sample names
        """

        import random
        import numpy as np
        import pandas as pd

        random.seed(self.seed)

        cv_folds = dict.fromkeys(range(0, self.num_folds))
        if self.holdout:
            cv_folds['holdout'] = pd.DataFrame(columns=['Status', 'Samples',
                                                        'Tissue'])
        for fold in cv_folds.keys():
            cv_folds[fold] = pd.DataFrame(columns=['Status', 'Samples',
                                                   'Tissue'])
        tissue_idx = 0
        for tissue in self.tissues:
            tissue_DF = self.y_dictionary[tissue]

            # Determine how many samples
            mut_samples = tissue_DF[tissue_DF['Status'] == 1]
            mut_samples = mut_samples['Samples'].tolist()
            not_mut_samples = tissue_DF[tissue_DF['Status'] == 0]
            not_mut_samples = not_mut_samples['Samples'].tolist()
            num_mut = len(mut_samples)
            num_not_mut = len(not_mut_samples)

            if self.holdout:
                # For hold out 3/4 are not mutated
                num_held = int(num_mut * self.fraction_held)
                hold_ix_mut = random.sample(range(0, num_mut), num_held)

                num_held = int(num_not_mut * self.fraction_held)
                hold_ix_not_mut = random.sample(range(0, num_not_mut),
                                                num_held)

                # Get the holdout samples
                mut_hold = [mut_samples[i] for i in hold_ix_mut]
                not_mut_hold = [not_mut_samples[i] for i in hold_ix_not_mut]
                accept = mut_hold + not_mut_hold
                hold_sam = tissue_DF[tissue_DF['Samples'].isin(accept)]

                # Randomly shuffle the rows
                hold_sam = hold_sam.iloc[np.random.permutation(len(hold_sam))]

                # Append them to the holdout key
                if tissue_idx == 0:
                    cv_folds['holdout'] = hold_sam
                else:
                    cv_folds['holdout'] = pd.concat([cv_folds['holdout'],
                                                     hold_sam])

                # Remove from tissue_DF
                tissue_DF = tissue_DF[~tissue_DF['Samples'].isin(accept)]

                # Now determine folds for the remaining samples
                mut_samples = tissue_DF[tissue_DF['Status'] == 1]
                mut_samples = mut_samples['Samples'].tolist()
                not_mut_samples = tissue_DF[tissue_DF['Status'] == 0]
                not_mut_samples = not_mut_samples['Samples'].tolist()
                num_mut = len(mut_samples)
                num_not_mut = len(not_mut_samples)

            # Because there are much less mutated samples than not, shuffle the
            # samples and assign manually to each fold
            range_folds = list(range(0, self.num_folds))
            random.shuffle(mut_samples)
            random.shuffle(range_folds)
            cv_fold_ix = 0
            samp_iter = 0
            cv_fold_samples = dict.fromkeys(range_folds)
            for key in cv_fold_samples.keys():
                cv_fold_samples[key] = []

            for samp in mut_samples:
                cur_fold = range_folds[cv_fold_ix]
                cv_fold_samples[cur_fold].append(samp)
                samp_iter += 1
                if samp_iter == int(num_mut / self.num_folds):
                    samp_iter = 0
                    cv_fold_ix += 1
                if cv_fold_ix > max(range_folds):
                    cv_fold_ix = range_folds[random.sample(range_folds, 1)[0]]

            # Non mutated samples are more prevalent, so assign randomly
            # without worrying too much about uneven numbers between folds
            fold_index = np.random.choice(range_folds, num_not_mut)
            for samp_idx in range(0, num_not_mut):
                fold = fold_index[samp_idx]
                samp = not_mut_samples[samp_idx]
                cv_fold_samples[fold].append(samp)

            # Reassign each of these samples to the actual sample information
            for key in cv_fold_samples.keys():
                cur_samples = cv_fold_samples[key]
                y_subset = tissue_DF[tissue_DF['Samples'].isin(cur_samples)]
                if tissue_idx == 0:
                    cv_folds[key] = y_subset
                else:
                    cv_folds[key] = pd.concat([cv_folds[key], y_subset])

            tissue_idx += 1

        self.cv_folds = cv_folds

    def assign_x_matrix_folds(self, full_x):
        """
        Uses previously established cross validation intervals to split full x
        """

        # Initialize the x_folds constructor
        self.x_folds = dict.fromkeys(self.cv_folds.keys())

        # Assign the x_matrix to each of the folds
        for cv_key in self.x_folds.keys():
            # Get the current fold
            cur_fold = self.cv_folds[cv_key]

            # Subset x matrix
            y_samples = cur_fold['Samples'].tolist()
            x_subset = full_x[y_samples]

            # Assign to dictionary
            self.x_folds[cv_key] = x_subset

    def split_train_test(self, fold):
        """
        Will partition x and y into training and testing sets with equal parts

        :param fold: the given fold to extract test set

        Output:
        training and testing x's and y's according to the subset fold
        """

        import pandas as pd

        # The test set is the holdout set if the fold == 'all'
        if fold == 'all':
            num_fold_iter = 0
            for f in range(self.num_folds):
                x_fold = self.x_folds[f]
                y_fold = self.cv_folds[f]
                if f != 'holdout':
                    if num_fold_iter == 0:
                        x_train = x_fold
                        y_train = y_fold['Status']
                        num_fold_iter += 1
                    else:
                        x_train = pd.concat([x_train, x_fold], axis=1)
                        y_train = pd.concat([y_train, y_fold['Status']])

            x_test = self.x_folds['holdout']
            y_test = self.cv_folds['holdout']['Status']

        # Remove a fold for the test set if fold not 'all'
        else:
            remain_folds = list(range(self.num_folds))
            remain_folds.remove(fold)

            # Concatenate the remaining x's and y's
            num_fold_iter = 0
            for f in remain_folds:
                x_fold = self.x_folds[f]
                y_fold = self.cv_folds[f]
                if num_fold_iter == 0:
                    x_train = x_fold
                    y_train = y_fold['Status']
                    num_fold_iter += 1
                else:
                    x_train = pd.concat([x_train, x_fold], axis=1)
                    y_train = pd.concat([y_train, y_fold['Status']])

            x_test = self.x_folds[fold]
            y_test = self.cv_folds[fold]['Status']

        return x_train.T, x_test.T, y_train, y_test
