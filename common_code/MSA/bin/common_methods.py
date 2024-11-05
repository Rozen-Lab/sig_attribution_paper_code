"""
Common methods module used across different scripts
"""
import os
import math
import json
import copy
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import distance
from statsmodels.stats.proportion import proportion_confint
import matplotlib.pyplot as plt

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which # python 2 option
    from shutil import which

    return which(name) is not None

def make_folder_if_not_exists(folder):
    if not os.path.exists(folder):
        try:
            os.makedirs(folder)
        except:
            warnings.warn("Could not create a folder %s" % folder)

def merge_plots(input_plots_filenames, merged_plot_filename, merge_tool = 'pdfunite'):
    if is_tool(merge_tool):
        if merge_tool=='pdfunite':
            merge_command = merge_tool + ' $(ls -v ' + input_plots_filenames + ') ' + merged_plot_filename
        elif merge_tool=='pdftk':
            merge_command = merge_tool + ' $(ls -v ' + input_plots_filenames + ') output ' + merged_plot_filename
        else:
            raise ValueError('Output plots merge failed: %s tool unknown (please use pdftk or pdfunite)' % merge_tool)
    else:
        warnings.warn('Output plots merge failed: %s tool not available' % merge_tool)
        return

    print('Merging output plots with command:', merge_command)
    try:
        os.system(merge_command)
    except:
        warnings.warn("Could not create a merged file %s. Please check if input files/folders exist." % merged_plot_filename)
    else:
        print("Resulting merged pdf: %s" % merged_plot_filename)

def write_data_to_JSON(data, JSON_output_file, json_orient = 'columns', indent = True):
    output_folder = JSON_output_file.rsplit('/',1)[0]
    make_folder_if_not_exists(output_folder)

    json_ready = False
    # convert pandas objects to JSON if needed
    if isinstance(data, pd.DataFrame) or isinstance(data, pd.Series):
        data_to_save = data.to_json(orient = json_orient)
        json_ready = True
    # account for pandas objects in dictionary
    elif isinstance(data, dict):
        data_to_save = {}
        for key, value in data.items():
            if isinstance(value, pd.DataFrame) or isinstance(value, pd.Series):
                data_to_save[key] = value.to_json(orient = json_orient)
            else:
                data_to_save[key] = value
    # assume other objects
    else:
        data_to_save = copy.deepcopy(data)

    output_file = open(JSON_output_file, 'w')

    if json_ready:
        output_file.write(data_to_save)
    else:
        output_file.write(json.dumps(data_to_save, indent=4 if indent else None, sort_keys = True))

    output_file.close()

def read_data_from_JSON(JSON_input_file, pandas = True, json_orient = 'columns'):
    input_file = open(JSON_input_file, 'r')
    input_JSON = ''.join(input_file.readlines())
    data = json.loads(input_JSON)
    if pandas:
        for key, value in data.items():
            # extract dataframe and convert index/columns to str
            df = pd.read_json(value, orient = json_orient, convert_axes = False)
            df.columns = df.columns.astype(str)
            df.index = df.index.astype(str)
            data[key] = df
    input_file.close()
    return data

def plot_array_as_histogram(arrays, labels, title, bins=[0.01*k for k in range(101)], savepath='./hist.pdf'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # bins = [0.025*k for k in range(42)]
    for array, label in zip(arrays, labels):
        plt.hist(array, label = label, histtype='step', bins = bins)

    legend = plt.legend()
    plt.gca().add_artist(legend)

    ax.set_title(title, fontsize=14, pad=10)
    # ax.set_xlabel('Attribution ($\%$)', fontsize=12)
    # ax.set_ylabel('Number of attributions', fontsize=12)
    # ax.set_xticks([0.1*k for k in range(11)])

    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

def calculate_confidence_interval(array, confidence=95):
    # return np.percentile(array, [(100-confidence)/2, (100+confidence)/2])
    return [np.percentile(array, (100-confidence)/2), np.percentile(array, (100+confidence)/2)]

def calculate_similarity(first_sample, second_sample, metric='Cosine', normalise=False):
    """
    Calculate various similarity metrics as 1 minus distance functions provided
    by scipy.spatial package (see https://docs.scipy.org/doc/scipy/reference/spatial.distance.html).
    If input samples are zeroes, the similarity can be undefined (nan) for some metrics.

    Parameters:
    ----------
    first_sample: array_like
        First input array (e.g. reconstructed sample)
    sedond_sample: array_like
        Second input array (e.g. true sample)
    metric: string
        One of the supported metrics: Jensen-Shannon, Cosine, Correlation, Chebyshev, L1, L2, L3
        (see https://docs.scipy.org/doc/scipy/reference/spatial.distance.html)
    normalise: boolean
        If true, normalise input vectors by the norm corresponding to chosen metric.
        L2 norm is considered for non-Minkowski metrics.

    Returns:
    -------
    similarity: double.
        Smaller distance gives similarity closer to one.
        Note: distance is not always normalised to 1 (e.g. JS vs minkowski metrics).
    -------
    """
    assert len(first_sample) == len(second_sample)
    if normalise:
        if 'L1' in metric or 'Manhattan' in metric or 'manhattan' in metric:
            if np.linalg.norm(first_sample, ord=1)!=0:
                first_sample = first_sample/np.linalg.norm(first_sample, ord=1)
            if np.linalg.norm(second_sample, ord=1)!=0:
                second_sample = second_sample/np.linalg.norm(second_sample, ord=1)
        elif 'L3' in metric:
            if np.linalg.norm(first_sample, ord=3)!=0:
                first_sample = first_sample/np.linalg.norm(first_sample, ord=3)
            if np.linalg.norm(second_sample, ord=3)!=0:
                second_sample = second_sample/np.linalg.norm(second_sample, ord=3)
        elif 'Chebyshev' in metric or 'chebyshev' in metric:
            if np.linalg.norm(first_sample, ord=np.inf)!=0:
                first_sample = first_sample/np.linalg.norm(first_sample, ord=np.inf)
            if np.linalg.norm(second_sample, ord=np.inf)!=0:
                second_sample = second_sample/np.linalg.norm(second_sample, ord=np.inf)
        else:
            # consider L2 norm
            if np.linalg.norm(first_sample, ord=2)!=0:
                first_sample = first_sample/np.linalg.norm(first_sample, ord=2)
            if np.linalg.norm(second_sample, ord=2)!=0:
                second_sample = second_sample/np.linalg.norm(second_sample, ord=2)
    if 'Cosine' in metric or 'cosine' in metric:
        similarity = 1 - distance.cosine(first_sample, second_sample)
    elif 'Correlation' in metric or 'correlation' in metric:
        similarity = 1 - distance.correlation(first_sample, second_sample)
    elif 'Chebyshev' in metric or 'chebyshev' in metric:
        similarity = 1 - distance.chebyshev(first_sample, second_sample)
    elif 'L1' in metric or 'Manhattan' in metric or 'manhattan' in metric:
        similarity = 1 - distance.minkowski(first_sample, second_sample, p=1)
    elif 'L2_normalised_by_first' in metric:
        if normalise:
            warnings.warn('Double normalisation using normalise flag with ', metric)
        similarity = 1 - distance.euclidean(first_sample, second_sample)/np.linalg.norm(first_sample, ord=2)
    elif 'L2' in metric or 'Euclidean' in metric or 'euclidean' in metric:
        similarity = 1 - distance.euclidean(first_sample, second_sample)
    elif 'L3' in metric:
        similarity = 1 - distance.minkowski(first_sample, second_sample, p=3)
    elif 'Jensen-Shannon' in metric or 'jensen-shannon' in metric:
        similarity = 1 - distance.jensenshannon(first_sample, second_sample)
    elif 'Canberra' in metric or 'canberra' in metric:
        similarity = 1 - distance.canberra(first_sample, second_sample)
    elif 'eucl_norm_squared' in metric:
        similarity = 1 - 0.5*np.var(first_sample-second_sample)/(np.var(first_sample)+np.var(second_sample))
    else:
        raise ValueError("Unknown metric: ", metric)

    return similarity

def calculate_confusion_matrix(signatures, reco_table, truth_table, number_of_samples = None):
    """
    Calculate confusion matrix

    Parameters:
    ----------
    signatures: list
        List of strings with signature names.
        Should be available in input reco and truth tables as columns.
    reco_table: pandas dataframe
        Dataframe containing the signature attribution table for a set of samples
    truth_table: pandas dataframe
        Dataframe containing the true signature activities (known from simulations)
        for the same set of samples (order has to be the same)
    number_of_samples: int or None
        If provided as an int, cap the number of samples by this number.
        Otherwise, all provided samples are used

    Returns:
    -------
    true_positives, true_negatives, false_positives, false_negatives: list of doubles
    -------
    """
    true_positives = true_negatives = false_positives = false_negatives = 0
    assert len(reco_table) == len(truth_table)

    true_positives_signatures = []
    false_negatives_signatures = []
    false_positives_signatures = []

    if not number_of_samples:
        number_of_samples = len(reco_table)

    for i in range(number_of_samples):
        reco_sample = reco_table.iloc[i]
        truth_sample = truth_table.iloc[i]
        for sig in signatures:
            if truth_sample[sig] > 0 and reco_sample[sig] > 0:
                true_positives += 1
                true_positives_signatures.append(truth_sample[sig])
            elif truth_sample[sig] == 0 and reco_sample[sig] == 0:
                true_negatives += 1
            elif truth_sample[sig] > 0 and reco_sample[sig] == 0:
                false_negatives += 1
                false_negatives_signatures.append(truth_sample[sig])
            elif truth_sample[sig] == 0 and reco_sample[sig] > 0:
                false_positives += 1
                false_positives_signatures.append(reco_sample[sig])
            else:
                raise ValueError(
                    "Invalid signature attribution values -- please check the tables")
    
    return true_positives, true_negatives, false_positives, false_negatives

def calculate_stat_scores(signatures, reco_table, truth_table, number_of_samples = None):
    """
    Calculate various statistical scores.

    Parameters:
    ----------
    signatures: list
        List of strings with signature names.
        Should be available in input reco and truth tables as columns.
    reco_table: pandas dataframe
        Dataframe containing the signature attribution table for a set of samples
    truth_table: pandas dataframe
        Dataframe containing the true signature activities (known from simulations)
        for the same set of samples (order has to be the same)
    number_of_samples: int or None
        If provided as an int, cap the number of samples by this number.
        Otherwise, all provided samples are used

    Returns:
    -------
    sensitivity, specificity, precision, accuracy, F1, MCC: list of doubles
    List of various scores (see e.g. https://en.wikipedia.org/wiki/Confusion_matrix for info)
    MCC stands for Matthews correlation coefficient.
    -------
    """
    
    true_positives, true_negatives, false_positives, false_negatives = calculate_confusion_matrix(signatures, reco_table, truth_table, number_of_samples)
    
    sensitivity = np.float64(true_positives) / (true_positives + false_negatives)
    specificity = np.float64(true_negatives) / (true_negatives + false_positives)
    precision = np.float64(true_positives) / (true_positives + false_positives)

    # # verbose info
    # print("Missed sigs mean/median/stdev/max:",np.mean(false_negatives_signatures),np.median(false_negatives_signatures),np.std(false_negatives_signatures),np.max(false_negatives_signatures))
    # plot_array_as_histogram([true_positives_signatures, false_positives_signatures, false_negatives_signatures], ['True positives', 'False positives', 'False negatives'], title = 'Signature attribution distributions', savepath="distributions.pdf")

    accuracy = np.float64(true_positives + true_negatives) / (true_positives + true_negatives + false_positives + false_negatives)
    F1 = 2 * np.float64(true_positives) / (2 * true_positives + false_positives + false_negatives)
    MCC = np.float64(true_positives * true_negatives - false_positives * false_negatives) / math.sqrt((true_positives + false_positives)
            * (true_positives + false_negatives) * (true_negatives + false_positives) * (true_negatives + false_negatives))

    return sensitivity, specificity, precision, accuracy, F1, MCC

def calculate_sensitivity_CI(signatures, reco_table, truth_table, number_of_samples = None, CI_method = 'beta'):
    """
    Calculate confidence interval on sensitivity as a binomial proportion

    Parameters:
    ----------
    signatures: list
        List of strings with signature names.
        Should be available in input reco and truth tables as columns.
    reco_table: pandas dataframe
        Dataframe containing the signature attribution table for a set of samples
    truth_table: pandas dataframe
        Dataframe containing the true signature activities (known from simulations)
        for the same set of samples (order has to be the same)
    number_of_samples: int or None
        If provided as an int, cap the number of samples by this number.
        Otherwise, all provided samples are used
    CI_method: string
        default: 'beta' (Clopper-Pearson) method to use for confidence interval in proportion_confint funciton.
        Supported methods: https://www.statsmodels.org/dev/generated/statsmodels.stats.proportion.proportion_confint.html

    Returns:
    -------
    sensitivity_CI: confidence interval as a list of doubles
    -------
    """
    true_positives, _, _, false_negatives = calculate_confusion_matrix(signatures, reco_table, truth_table, number_of_samples)

    # calculate CIs on sensitivity and specificity
    sensitivity_CI = proportion_confint(true_positives, true_positives + false_negatives, alpha=0.05, method=CI_method)

    return sensitivity_CI

def calculate_specificity_CI(signatures, reco_table, truth_table, number_of_samples = None, CI_method = 'beta'):
    """
    Calculate confidence interval on specificity as a binomial proportion

    Parameters:
    ----------
    signatures: list
        List of strings with signature names.
        Should be available in input reco and truth tables as columns.
    reco_table: pandas dataframe
        Dataframe containing the signature attribution table for a set of samples
    truth_table: pandas dataframe
        Dataframe containing the true signature activities (known from simulations)
        for the same set of samples (order has to be the same)
    number_of_samples: int or None
        If provided as an int, cap the number of samples by this number.
        Otherwise, all provided samples are used
    CI_method: string
        default: 'beta' (Clopper-Pearson) method to use for confidence interval in proportion_confint funciton.
        Supported methods: https://www.statsmodels.org/dev/generated/statsmodels.stats.proportion.proportion_confint.html

    Returns:
    -------
    specificity_CI: confidence interval as a list of doubles
    -------
    """
    _, true_negatives, false_positives, _ = calculate_confusion_matrix(signatures, reco_table, truth_table, number_of_samples)

    # calculate CIs on specificity
    specificity_CI = proportion_confint(true_negatives, true_negatives + false_positives, alpha=0.05, method=CI_method)

    return specificity_CI

def clean_up_labels(signature, category, mutation_type='SBS'):
    # clean up labels for nicer plots
    sig_clone = copy.deepcopy(signature)

    if type(sig_clone.index) != pd.MultiIndex:
        if mutation_type=='SBS':
            sig_clone.index = sig_clone.index.str.replace('\[' + category + '\]', category[0])
        elif mutation_type=='DBS':
            sig_clone.index = sig_clone.index.str.replace(category, '')
        elif mutation_type=='ID':
            sig_clone.index = sig_clone.index.str.replace(category, '')
            sig_clone.index = sig_clone.index.str.replace('_', '')

            # special treatment for some indel categories - removes '+' sign in front of integers
            if '+' in category:
                sig_clone.index = pd.Series(sig_clone.index).apply(lambda x : x[1:] if x.startswith("+") else x)

            # for certain indel categories, increase the last integer by 1
            if category in ['DEL_C_1', 'DEL_T_1', 'DEL_repeats_2', 'DEL_repeats_3', 'DEL_repeats_4', 'DEL_repeats_5+']:
                sig_clone.index = pd.Series(sig_clone.index).apply(lambda x : str(int(x[0])+1)+'+' if x.endswith("+") else str(int(x)+1))
    else:
        sig_clone = sig_clone.droplevel(0)

    sig_clone.index.name = ''
    return sig_clone

def categorise_mutations(input_data, categories, mutation_type='SBS', condensed=False, strand_bias=False, add_non_transcribed=False):
    """
    Categorise the dataframe by a given a list of categories

    Parameters:
    ----------
    input_data: pandas dataframe
        Input dataframe. If strand_bias flag is True, a single-column dataframe is expected with
        strand information (see rearrange_stranded_data function)

    mutation_type: string ('SBS' by default)
        Mutation type used in clean_up_labels function.

    categories: list of strings
        Mutation categories to condense mutations to, e.g.: ['C>A','C>G','C>T','T>A','T>C','T>G']

    condensed: boolean
        If True, the function returns a dataframe with mutations summed up for given categories,
        otherwise, a dictionary of dataframes is returned with keys for each category

    strand_bias: boolean
        If True, rearrange_stranded_data function is called on a dataframe to obtain columns on strands.

    add_non_transcribed: boolean flag
        Add the non-transcribed region (flag is passed to rearrange_stranded_data function)

    Returns:
    -------
    categorised_mutations: pandas dataframe
        Dataframe with categorised mutations if condensed flag is True, otherwise
        a dictionary of dataframes with keys for each category listed in input categories
    -------
    """
    data = copy.deepcopy(input_data)
    if strand_bias:
        data = rearrange_stranded_data(data, add_non_transcribed = add_non_transcribed)

    # dataframe in case of condensed format, dictionary of dataframes otherwise
    if condensed:
        if strand_bias:
            if add_non_transcribed:
                categorised_mutations = pd.DataFrame(columns=['transcribed strand','untranscribed strand','non-transcribed'])
            else:
                categorised_mutations = pd.DataFrame(columns=['transcribed strand','untranscribed strand'])
        else:
            if isinstance(data, pd.Series):
                categorised_mutations = pd.DataFrame(columns=['N'])
            else:
                categorised_mutations = pd.DataFrame(columns=data.columns)
    else:
        categorised_mutations = {}

    for category in categories:
        mutations_per_category = clean_up_labels(data.filter(like=category, axis=0), category, mutation_type=mutation_type)

        # Sum up identical categories (if MultiIndex)
        if mutations_per_category.index.nlevels>1:
            mutations_per_category = mutations_per_category.groupby(level=1).sum()

        # sum up categories for condensed format, otherwise leave it in dictionary
        if condensed:
            categorised_mutations.loc[category] = mutations_per_category.sum()
        else:
            categorised_mutations[category] = mutations_per_category

    return categorised_mutations

def rearrange_stranded_data(input_data, add_non_transcribed = False):
    """
    Rearrange stranded data of a single sample into two different columns of the dataframe

    Parameters:
    ----------
    input_data: pandas dataframe
        Supported types:
        1. Single-column dataframe with a 3-level MultiIndex of 192-context PCAWG format (e.g. WGS_PCAWG.192.csv on Synapse).
        Essential to have a 0-level of the MultiIndex, with 'T' and 'U' values.
        2. Single-column dataframe with regular index.
        Essential to have categories ending with "-transcribed" and "-untranscribed" to filter upon,
        or starting with "T:" or "U:" (SigProfiler output).
    add_non_transcribed: boolean flag
        Add the non-transcribed region (useful in contexts higher than 192)
        Essential to have a 0-level of the MultiIndex with 'N' value, or categories
        ending with "-nontranscribed" to filter upon, or starting with "N:" (SigProfiler output).
    Returns:
    -------
    rearranged_data: pandas dataframe
        Dataframe with removed strand info in index, but with two columns
        (transcribed and untranscribed strands), or three columns if non-transcribed region is required.
    -------
    """
    data = copy.deepcopy(input_data)
    T_strand, U_strand = pd.DataFrame(), pd.DataFrame()
    if type(input_data.index) == pd.MultiIndex:
        T_strand = data.iloc[data.index.get_level_values(0) == 'T']
        U_strand = data.iloc[data.index.get_level_values(0) == 'U']
        T_strand = T_strand.droplevel(0)
        U_strand = U_strand.droplevel(0)
        if add_non_transcribed:
            N_strand = data.iloc[data.index.get_level_values(0) == 'N']
            N_strand = N_strand.droplevel(0)
    else:
        if not data.filter(like='T:').empty and not data.filter(like='U:').empty:
            T_strand = data.filter(like='T:')
            U_strand = data.filter(like='U:')
            T_strand.index = T_strand.index.str.replace('T:', '')
            U_strand.index = U_strand.index.str.replace('U:', '')
            if add_non_transcribed:
                N_strand = data.filter(like='N:')
                N_strand.index = U_strand.index.str.replace('N:', '')
        elif not data.filter(like='-transcribed').empty and not data.filter(like='-untranscribed').empty:
            T_strand = data.filter(like='-transcribed')
            U_strand = data.filter(like='-untranscribed')
            T_strand.index = T_strand.index.str.replace('-transcribed', '')
            U_strand.index = U_strand.index.str.replace('-untranscribed', '')
            if add_non_transcribed:
                N_strand = data.filter(like='-nontranscribed')
                N_strand.index = U_strand.index.str.replace('-nontranscribed', '')
        else:
            raise ValueError("Can't find transcribed or untranscribed strands, please check input data:", data.index)
    if add_non_transcribed:
        rearranged_data = pd.concat([T_strand, U_strand, N_strand], names=['transcribed strand','untranscribed strand', 'non-transcribed'], axis=1, sort=True)
        rearranged_data.columns = ['transcribed strand', 'untranscribed strand', 'non-transcribed']
    else:
        rearranged_data = pd.concat([T_strand, U_strand], names=['transcribed strand','untranscribed strand'], axis=1, sort=True)
        rearranged_data.columns = ['transcribed strand', 'untranscribed strand']
    return rearranged_data

def mutational_burden(data, absolute=True):
    # calculate the absolute total number of mutations in the dataset
    if isinstance(data, dict):
        burden = 0
        for key, dataframe in data.items():
            burden += mutational_burden(dataframe, absolute=absolute)
        return burden
    elif isinstance(data, pd.DataFrame):
        if absolute:
            return abs(sum(data.sum()))
        else:
            return sum(data.sum())
    elif isinstance(data, pd.Series):
        if absolute:
            return abs(data.sum())
        else:
            return data.sum()
    else:
        raise ValueError("Unknown data type:",type(data))

def normalise_mutations(input_data):
    # normalise the number of mutations by the total number of mutations
    data = copy.deepcopy(input_data)
    total_mutational_burden = mutational_burden(data)
    if isinstance(data, pd.DataFrame):
        data = 100*data/total_mutational_burden
    elif isinstance(data, dict):
        for category, mutations in data.items():
            data[category] = 100*mutations/total_mutational_burden
    else:
        raise ValueError("Unknown data type:",type(data))

    return data

def calculate_poisson_CL(x, confidence=0.95):
    # calculate Poisson confidence limits
    confidence_limit = [stats.chi2.ppf((1-confidence)/2, 2*x)/2, stats.chi2.ppf((1+confidence)/2, 2*(x+1))/2]
    return confidence_limit

def calculate_errors(input_data, normalise = False, normalisation_factor = None):
    # calculate confidence intervals for error bars (assuming Poisson)
    if isinstance(input_data, dict):
        errors = {}
        for key, dataframe in input_data.items():
            errors[key] = calculate_errors(dataframe, normalise = normalise, normalisation_factor = mutational_burden(input_data))
        return errors
    elif isinstance(input_data, pd.DataFrame) or isinstance(input_data, pd.Series):
        # ensure numeric format
        if isinstance(input_data, pd.DataFrame):
            data = input_data.astype('float')
        else:
            data = pd.to_numeric(input_data)
        # calculate Poisson confidence limits
        confidence_limit = calculate_poisson_CL(data.T)
        # treat NaNs
        confidence_limit = np.nan_to_num(confidence_limit)
        # calculate errors (offsets) for plotting
        errors = np.abs(data.T.to_numpy()-confidence_limit)
        # swap axes for matplotlib format
        if not isinstance(input_data, pd.Series):
            errors = np.swapaxes(errors, 0, 1)
        if normalise:
            if not normalisation_factor:
                normalisation_factor = mutational_burden(input_data)
            errors = 100*errors/normalisation_factor
        return errors
    else:
        raise ValueError("Unknown data type:",type(input_data))
