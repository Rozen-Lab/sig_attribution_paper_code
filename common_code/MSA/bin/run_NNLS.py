import argparse
import os
import warnings
import pandas as pd
import numpy as np
import copy
import time
from scipy.optimize import nnls
from scipy.spatial import distance
from common_methods import make_folder_if_not_exists, calculate_similarity

def bootstrap_mutation_table(input_dataframe, method="classic", fitted=None, residuals=None):
    input_mutations = copy.deepcopy(input_dataframe)
    if method=="classic":
        # classic bootstrap (resampling with replacement, i.i.d. assumption)
        original_index = input_mutations.index
        # bootstrap using pd.sample
        bootstrap_mutations = input_mutations.sample(n=len(input_mutations.index), replace=True)
        # keep the old index (mutation categories)
        bootstrap_mutations.index = original_index
    elif method=="bootstrap_residuals":
        # semi-parametric model of bootstrapping residuals
        # input mutations dataframe is not actually used, as provided fitted
        # dataframe is added with resampled residuals dataframe
        if residuals is None:
            raise ValueError("Residuals are not provided for %s method, please use the residuals flag" % method)
        elif fitted is None:
            raise ValueError("Fitted values are not provided for %s method, please use the fitted flag" % method)
        else:
            # bootstrap residuals and add them to input mutations
            return fitted + bootstrap_mutation_table(residuals, method="classic")
    elif method=="poisson":
        # Poisson bootstrap (n>100 and i.i.d. assumption)
        bootstrap_mutations = input_mutations.applymap(lambda x: x*np.random.poisson(1))
    elif method=="binomial" or "multinomial" in method:
        # mutational burdens for each sample (column) for normalisation
        sums = input_mutations.sum(axis=0)
        normalised_input_mutations = input_mutations.div(sums, axis=1)
        # in case if mutational burden is zero for some samples, replace nans with zeros
        normalised_input_mutations.fillna(0, inplace=True)
        bootstrap_mutations = pd.DataFrame(index = input_mutations.index, columns = input_mutations.columns)
        for i in range(len(bootstrap_mutations.columns)):
            if method=="binomial":
                # simple independent binomials for each category
                bootstrap_mutations.iloc[:,i] = np.random.binomial(sums.iloc[i], normalised_input_mutations.iloc[:,i])
            elif method=="multinomial":
                # Weighted multinomial with fixed number of draws equalling mutational burden
                bootstrap_mutations.iloc[:,i] = np.random.multinomial(sums.iloc[i], normalised_input_mutations.iloc[:,i])
            elif method=="multinomial_weight":
                # Mulinomial counts with equal probabilities, equivalent to classic bootstrap
                n = len(normalised_input_mutations.iloc[:,i])
                counts = np.random.multinomial(n, [1/n for i in range(n)])
                input_data = input_mutations.iloc[:,i]
                bootstrap = np.concatenate([[input_data[j]]*counts[j] for j in range(len(input_data))]).ravel().astype(int)
                np.random.shuffle(bootstrap)
                bootstrap_mutations.iloc[:,i] = bootstrap
            else:
                raise ValueError("Unknown bootstrap method: %s" % method)
    else:
        raise ValueError("Unknown bootstrap method: %s" % method)
    return bootstrap_mutations

def remove_weak_signatures(selected_mutations, initial_signatures, weak_threshold = 0.01, similarity_index=-3, verbose = False):
    """
    Perform a loop to remove all weak signatures with given penalty (weak_threshold).
    The NNLS is intially applied on a given input set of signatures and mutations.
    The similarity of the reconstructed profile (linear combination of input signatures)
    to the initial sample is calculated, called base similarity (column with index similarity_index in stat_info dataframe).
    The actual loop to remove weak signatures is executed, in which all signatures that decrease the
    similarity by less than a given threshold (weak_threshold) are removed.
    The output is the final list of remaining signatures.

    Parameters:
    ----------
    selected_mutations: list of ints
        List of mutation counts for a single sample. Example: list of 96 integers,
        for a 96-context SBS mutations.

    initial_signatures: pandas dataframe
        Dataframe of input signatures with the index of mutation categories, with the order
        that has to be the same as for the input mutations.

    weak_threshold: float
        Similarity decrease threshold to exclude weakest signatures (default: 0.01)

    similarity_index: into
        Index of the similarity metric (e.g. 3rd from the end of stat_info column list: -3)

    verbose: boolean
        Verbosity flag: lots of output for debugging if set to True

    Returns:
    -------
    significant_signatures: list of strings
        The list of remaining signatures passing the optimisation loop.
    -------
    """

    significant_signatures = copy.deepcopy(initial_signatures)
    # calculate the base similarity with an input set of signatures
    _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    base_similarity = stat_info[similarity_index]

    if verbose:
        print('Starting the removal loop.')
        print('Current significant signatures:')
        print(significant_signatures.columns.tolist())
        print('Current base similarity:', base_similarity)

    all_weak_signatures_removed = False
    while all_weak_signatures_removed == False:
        signatures_contribution = {}
        if verbose:
            print('Current significant signatures:')
            print(significant_signatures.columns.tolist())
        for signature in significant_signatures.columns.tolist():
            signatures_to_run = copy.deepcopy(significant_signatures)
            signatures_to_run = signatures_to_run.drop(signature, axis=1)
            if verbose:
                print('Running without signature:', signature)
            _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, signatures_to_run)
            new_similarity = stat_info[similarity_index]
            contribution = base_similarity - new_similarity
            signatures_contribution[signature] = contribution

        if verbose:
            print('Signatures contribution for sample', sample)
            print(signatures_contribution)

        least_contributing_signature = min(signatures_contribution, key=signatures_contribution.get)
        if signatures_contribution[least_contributing_signature]<weak_threshold and len(significant_signatures.columns.tolist())>1:
            significant_signatures = significant_signatures.drop(least_contributing_signature, axis=1)
            if verbose:
                print('Dropping signature %s from sample %s' % (least_contributing_signature, sample))
                print('Number of signatures left:',len(significant_signatures.columns.tolist()))
            # recalculate the base similarity with a new set of signatures
            _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
            base_similarity = stat_info[similarity_index]
        else:
            if verbose:
                print('All cleaned up!')
            all_weak_signatures_removed = True

        # check if only one significant signature is left
        if len(significant_signatures.columns.tolist())==1:
            if verbose:
                print('Only one signature left:',significant_signatures.columns.tolist())
            all_weak_signatures_removed = True

    # calculate the new similarity with a final set of signatures
    _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    final_similarity = stat_info[similarity_index]

    if verbose:
        print('Finished the removal loop:')
        print('Significant signatures after removal loop:')
        print(significant_signatures.columns.tolist())
        print('Similarity after removal loop:', final_similarity)
        print('*'*30)

    return significant_signatures, final_similarity

def add_strong_signatures(selected_mutations, initial_signatures, all_available_signatures, strong_threshold = 0.05, similarity_index=-3, verbose = False):
    """
    Perform a loop to remove all weak signatures with given penalty (weak_threshold).
    The NNLS is intially applied on a given input set of signatures and mutations.
    The similarity of the reconstructed profile (linear combination of input signatures)
    to the initial sample is calculated, called base similarity (column with index similarity_index in stat_info dataframe).
    The actual loop to add strong signatures is executed, in which all excluded so far signatures
    that increase the similarity by more than a given threshold (strong_threshold) are included again.
    The output is the final list of remaining signatures.

    Parameters:
    ----------
    selected_mutations: list of ints
        List of mutation counts for a single sample. Example: list of 96 integers,
        for a 96-context SBS mutations.

    initial_signatures: pandas dataframe
        Dataframe of input signatures with the index of mutation categories, with the order
        that has to be the same as for the input mutations.

    all_available_signatures: pandas dataframe
        Dataframe of all available signatures with the index of mutation categories, with the order
        that has to be the same as for the input mutations.

    strong_threshold: float
        Similarity increase threshold to include strongest signatures (default: 0.05)

    similarity_index: into
        Index of the similarity metric (e.g. 3rd from the end of stat_info column list: -3)

    verbose: boolean
        Verbosity flag: lots of output for debugging if set to True

    Returns:
    -------
    significant_signatures: list of strings
        The list of remaining signatures passing the optimisation loop.
    -------
    """

    significant_signatures = copy.deepcopy(initial_signatures)
    # calculate the base similarity with an input set of signatures
    _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    # global similarities_L2, similarities_L2_normed_by_1st, difference_in_norm
    base_similarity = stat_info[similarity_index]
    # remaining signatures to loop through in search for most contributing ones
    remaining_signatures = all_available_signatures.drop(initial_signatures.columns, axis=1)

    if verbose:
        print('Starting the addition loop.')
        print('Current significant signatures:')
        print(significant_signatures.columns.tolist())
        print('Current base similarity:', base_similarity)

    # a loop to add strong signatures if there is anything to add
    if remaining_signatures.empty:
        all_strong_signatures_added = True
    else:
        all_strong_signatures_added = False

    while all_strong_signatures_added == False:
        signatures_contribution = {}
        if verbose:
            print('Current significant signatures:')
            print(significant_signatures.columns.tolist())

        # loop through the remaining signatures to add contributing ones
        for signature in remaining_signatures.columns.tolist():
            signatures_to_run = copy.deepcopy(significant_signatures)
            signatures_to_run[signature] = pd.Series(remaining_signatures[signature], index = signatures_to_run.index)
            if verbose:
                print('Running with signature:', signature)
            _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, signatures_to_run)
            new_similarity = stat_info[similarity_index]
            contribution = new_similarity - base_similarity
            signatures_contribution[signature] = contribution

        if verbose:
            print('Signatures contribution for sample', sample)
            print(signatures_contribution)

        most_contributing_signature = max(signatures_contribution, key = signatures_contribution.get)
        if signatures_contribution[most_contributing_signature]>strong_threshold:
            if verbose:
                print('Adding signature %s to sample %s' % (most_contributing_signature, sample))
            significant_signatures[most_contributing_signature] = pd.Series(remaining_signatures[most_contributing_signature], index = significant_signatures.index)
            remaining_signatures = remaining_signatures.drop(most_contributing_signature, axis=1)
            # recalculate the base similarity with a new set of signatures
            _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
            base_similarity = stat_info[similarity_index]
            if remaining_signatures.empty:
                if verbose:
                    print('All strong signatures added!')
                all_strong_signatures_added = True
        else:
            if verbose:
                print('All strong signatures added!')
            all_strong_signatures_added = True

    # sort signatures
    significant_signatures = significant_signatures.reindex(sorted(significant_signatures.columns), axis=1)

    # calculate the new similarity with a final set of signatures
    _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    final_similarity = stat_info[similarity_index]

    if verbose:
        print('Finished the addition loop.')
        print('Significant signatures after addition loop:')
        print(significant_signatures.columns.tolist())
        print('Similarity after addition loop:', final_similarity)
        print('*'*30)

    return significant_signatures, final_similarity

def optimise_signatures(selected_mutations, initial_signatures, all_available_signatures, strategy='removal', weak_threshold=0.01, strong_threshold=0.05, similarity_index=-3, loops_limit=100, verbose=False):
    """
    Perform signature optimisation for NNLS attribution method.
    The method is outlined in the PCAWG paper, with the main idea as follows:
    The NNLS is intially applied on a given input set of signatures and mutations.
    The similarity of the reconstructed profile (linear combination of input signatures)
    to the initial sample is calculated, called base similarity (column with index similarity_index in stat_info dataframe).
    Afterwards, a loop of addition and/or removal of signatures is executed, with set penalties, until convergence is reached.
    The strategy can be defined as 'removal', 'addition' or 'add-remove', determining the nature of loops.
    The output is the final list of remaining signatures.

    Parameters:
    ----------
    selected_mutations: list of ints
        List of mutation counts for a single sample. Example: list of 96 integers,
        for a 96-context SBS mutations.

    initial_signatures: pandas dataframe
        Dataframe of input signatures with the index of mutation categories, with the order
        that has to be the same as for the input mutations.

    all_available_signatures: pandas dataframe
        Dataframe of all available signatures with the index of mutation categories, with the order
        that has to be the same as for the input mutations.

    strategy: str
        String parameter defining the strategy: 'removal' (default), 'addition' or 'add-remove',
        determining the method of executing signature addition and/or removal loops.

    weak_threshold: float
        Similarity decrease threshold to exclude weakest signatures (default: 0.01)

    strong_threshold: float
        Similarity increase threshold to include strongest signatures (default: 0.05)

    similarity_index: int
        Index of the similarity metric (e.g. 3rd from the end of stat_info column list: -3)

    loops_limit: int
        Maximum number of add-remove iterations. If reached, optimisation is finished (default: 100).

    verbose: boolean
        Verbosity flag: lots of output for debugging if set to True

    Returns:
    -------
    final_signatures: list of strings
        The final list of remaining signatures passing the optimisation.
    -------
    """

    significant_signatures = copy.deepcopy(initial_signatures)

    # drop all signatures to start from zero
    # significant_signatures.drop(significant_signatures.columns, axis=1, inplace=True)

    # calculate the base similarity with an input set of signatures
    _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    base_similarity = stat_info[similarity_index]


    if verbose:
        print('Initial set of signatures:')
        print(significant_signatures.columns.tolist())
        print('Initial (base) similarity:', base_similarity)

    # Signature optimisation routine starts here
    if strategy=='removal':
        significant_signatures, converging_similarity = remove_weak_signatures(selected_mutations, significant_signatures,
                                                    weak_threshold = weak_threshold, similarity_index = similarity_index, verbose = verbose)
    elif strategy=='addition':
        significant_signatures, converging_similarity = add_strong_signatures(selected_mutations, significant_signatures, all_available_signatures,
                                                    strong_threshold = strong_threshold, similarity_index = similarity_index, verbose = verbose)
    elif strategy=='add-remove':
        convergence_delta = 1
        converging_similarity = 1
        loop_counter = 0
        while convergence_delta > 0 and loop_counter < loops_limit:
            loop_counter += 1
            if converging_similarity:
                base_similarity = converging_similarity
                significant_signatures, converging_similarity = add_strong_signatures(selected_mutations, significant_signatures, all_available_signatures,
                                                        strong_threshold = strong_threshold, similarity_index = similarity_index, verbose = verbose)
                significant_signatures, converging_similarity = remove_weak_signatures(selected_mutations, significant_signatures,
                                                        weak_threshold = weak_threshold, similarity_index = similarity_index, verbose = verbose)
            convergence_delta = abs(converging_similarity-base_similarity)
            if loop_counter >= loops_limit:
                warnings.warn("Maximum number of iterations (%i) reached, stopping optimisation." % loops_limit)
            if verbose:
                print('Loop', loop_counter, 'done.')
                print('Convergence delta: ', convergence_delta)
    else:
        raise ValueError('Unknown strategy: please choose from removal, addition and add-remove')

    final_signatures = significant_signatures

    if verbose:
        if strategy=='add-remove':
            print('Number of add-remove loops:', loop_counter)
        print('Final set of signatures:')
        print(final_signatures.columns.tolist())
        print('Final similarity:', converging_similarity)

    return final_signatures

def perform_signature_attribution(selected_mutations, signatures, normalise_mutations=False, verbose=False):
    """
    Apply the NNLS attribution method using NNLS function from scipy.optimize:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html

    Parameters:
    ----------
    selected_mutations: list of ints
        List of mutation counts for a single sample. Example: list of 96 integers,
        for a 96-context SBS mutations. The order has to be identical to input signatures

    signatures: pandas dataframe
        Dataframe of input signatures with the index of mutation categories, with the order
        that has to be the same as for the input mutations

    normalise_mutations: boolean
        Boolean flag, if set to true the mutation counts normalised to the input mutational burden
        will be returned. False by default.

    verbose: boolean
        Verbosity flag: lots of output for debugging if set to True

    Returns:
    -------
    normalised_weights, mutation_numbers, stat_info: tuple
    normalised_weights: list of floats
        The list weights normalised by the sum of weights (i.e. all add up to one).
        Corresponds to a probability of all signatures to contribute to the given sample
    mutation_numbers: list of floats
        The list of normalised weights mentioned above multiplied by total mutational burden in case if
        normalise_mutations flag is set to True. Otherwise, mutation counts as extracted by NNLS.
    stat_info: a list of floats
        A list of statistical info variables, as follows:
        input_mutational_burden: Total mutational burden of input sample
        rss: residual sum of squares (RSS),
        chi2: Chi-2 of the fit,
        r2: determination coefficient (R2)
        cosine_similarity: cosine similarity of the reconstructed profile with the input sample
        correlation: correlation of the reconstructed profile with the input sample
        chebyshev_similarity: Chebyshev similarity for mentioned profiles
        L1/L2/L3 similarity: 1 minus L1/L2/L3 Minkowski distances between mentioned profiles
        jensenshannon_similarity: Jensen-Shannon similarity of mentioned profiles
    -------
    """
    if signatures.empty:
        if verbose:
            print('Zero signatures provided to NNLS.')
        return np.nan, np.nan, np.nan, np.nan, [sum(selected_mutations), np.nan, np.nan, np.nan, 0, 0, 0, 0, 0, 0]

    reg = nnls(signatures, selected_mutations)
    weights = reg[0]

    normalised_weights = weights/sum(weights)
    if normalise_mutations:
        mutation_numbers = normalised_weights*sum(selected_mutations)
    else:
        mutation_numbers = weights

    # calculate RSS, chi2, r2, similarity metrics
    observed = selected_mutations
    fitted = np.matmul(signatures.values, weights)
    residuals = observed - fitted
    chi2 = sum(residuals*residuals/fitted)

    mean = np.mean(observed)
    tot_squares = (observed-mean)*(observed-mean)
    tss = sum(tot_squares)
    rss = sum(residuals*residuals)
    r2 = 1-rss/tss

    input_mutational_burden = sum(selected_mutations)
    cosine_similarity = calculate_similarity(observed, fitted)
    correlation = calculate_similarity(observed, fitted, metric='Correlation')
    chebyshev_similarity = calculate_similarity(observed, fitted, metric='Chebyshev', normalise = True)
    L1_similarity = calculate_similarity(observed, fitted, metric='L1', normalise = True)
    # L2_similarity = calculate_similarity(observed, fitted, metric='L2', normalise = True)
    L2_similarity = calculate_similarity(observed, fitted, metric='L2_normalised_by_first')
    L3_similarity = calculate_similarity(observed, fitted, metric='L3', normalise = True)
    jensenshannon_similarity = calculate_similarity(observed, fitted, metric='jensen-shannon')

    stat_info = [input_mutational_burden, rss, chi2, r2, cosine_similarity, correlation, chebyshev_similarity, L1_similarity, L2_similarity, L3_similarity, jensenshannon_similarity]

    if verbose:
        if chi2>10e10:
            print('************* High chi2 sample *************')
        print('Signatures:')
        print(signatures.columns.tolist())
        print('[input_mutational_burden, rss, chi2, r2, cosine_similarity, correlation, chebyshev_similarity, L1_similarity, L2_similarity, L3_similarity, jensenshannon_similarity]:',stat_info)
        print('Observed:', observed)
        print('Sum observed:', np.sum(observed))
        print('Fitted:', fitted)
        print('Sum Fitted:', np.sum(fitted))
        print('Residuals:', residuals)

    return normalised_weights, mutation_numbers, fitted, residuals, stat_info

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='SIM',
                        help="set the dataset name (SIM by default)")
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='',
                        help="set mutation type (SBS, DBS, ID, SV, CNV)")
    parser.add_argument("-c", "--context", dest="context", default=96, type=int,
                        help="set SBS context (96, 192, 288, 1536)")
    parser.add_argument("-i", "--input_path", dest="input_path", default='input_mutation_tables/',
                        help="set path to input mutation tables")
    parser.add_argument("-s", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                        help="set path to signature tables")
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                        help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_argument("-o", "--output_path", dest="output_path", default='output_tables/',
                        help="set path to save output tables")
    parser.add_argument("-x", "--optimise_signatures", dest="optimise_signatures", action="store_true",
                        help="perform signature optimisation (remove weak signatures, add strong ones)")
    parser.add_argument("-W", "--weak_threshold", dest="weak_threshold", default=0.01, type=float,
                        help="Similarity decrease threshold to exclude weakest signatures in optimisation (default: 0.01)")
    parser.add_argument("-S", "--strong_threshold", dest="strong_threshold", default=0.01, type=float,
                        help="Similarity increase threshold to include strongest signatures in optimisation (default: 0.01)")
    parser.add_argument("-N", "--normalise_mutations", dest="normalise_mutations", action="store_true",
                        help="Normalise mutation counts to the input mutational burden in NNLS output")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="verbosity flag for debugging (lots of output)")
    parser.add_argument("-n", "--number", dest="number_of_samples", default=-1, type=int,
                        help="limit the number of samples to analyse (all by default)")
    parser.add_argument("-B", "--bootstrap", dest="bootstrap", action="store_true",
                        help="Reshuffle input mutation tables for bootstrapping")
    parser.add_argument("--bootstrap_method", dest="bootstrap_method", default='binomial',
                        help="choose a method for bootstrapping (perturbing) samples (classic, binomial, multinomial)")
    parser.add_argument("--add_suffix", dest="add_suffix", action="store_true",
                        help="add suffixes to output folders (useful during optimisation scan analysis)")
    parser.add_argument("--optimisation_strategy", dest="optimisation_strategy", default='removal',
                        help="set optimisation strategy (removal by default, addition or add-remove)")

    args = parser.parse_args()

    dataset_name = args.dataset_name
    mutation_type = args.mutation_type
    context = args.context
    input_path = args.input_path
    signature_tables_path = args.signature_tables_path
    signatures_prefix = args.signatures_prefix
    output_path = args.output_path + '/' + dataset_name

    if args.add_suffix:
        output_path += '_%i_NNLS' % context
        if args.optimise_signatures:
            output_path += '_%.4f_%.4f' % (args.weak_threshold, args.strong_threshold)
        else:
            output_path += '_unoptimised'

    make_folder_if_not_exists(output_path)

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (or '-t DBS', '-t ID').")
    elif mutation_type not in ['SBS', 'DBS', 'ID', 'SV', 'CNV']:
        raise ValueError("Error: Unknown mutation type: %s. Known types: SBS, DBS, ID, SV, CNV" % mutation_type)

    if mutation_type=='SBS':
        if context==96:
            signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=[0,1])
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_path, dataset_name, dataset_name, context), sep=None, index_col=[0,1])
        elif context in [192, 288]:
            signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type, context), sep=None, index_col=[0,1,2])
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_path, dataset_name, dataset_name, context), sep=None, index_col=[0,1,2])
        elif context==1536:
            signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type, context), sep=None, index_col=0)
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_path, dataset_name, dataset_name, context), sep=None, index_col=0)
        else:
            raise ValueError("Context %i is not supported." % context)
    elif mutation_type=='DBS':
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=0)
        input_mutations = pd.read_csv('%s/%s/WGS_%s.dinucs.csv' % (input_path, dataset_name, dataset_name), sep=None, index_col=0)
    elif mutation_type=='ID':
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=0)
        input_mutations = pd.read_csv('%s/%s/WGS_%s.indels.csv' % (input_path, dataset_name, dataset_name), sep=None, index_col=0)
    else: # SV and CNV
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=0)
        input_mutations = pd.read_csv('%s/%s/WGS_%s.%s.csv' % (input_path, dataset_name, dataset_name, mutation_type), sep=None, index_col=0)

    print("Performing NNLS for %s dataset, %s mutation type." % (dataset_name, mutation_type))
    if mutation_type=='SBS':
        print("SBS context: %i" % context)
    if args.optimise_signatures:
        print("Optimised NNLS method, weak/strong thresholds: %f/%f" % (args.weak_threshold, args.strong_threshold))
    if args.bootstrap:
        print("Perturbing the input mutation sample, bootstrap (simulation) method: %s" % args.bootstrap_method)

    if args.bootstrap and args.bootstrap_method != "bootstrap_residuals":
        input_mutations = bootstrap_mutation_table(input_mutations, method=args.bootstrap_method)

    # limit the number of samples to analyse (if specified by -n option)
    if args.number_of_samples!=-1:
        input_mutations = input_mutations.iloc[:,0:args.number_of_samples]

    samples = input_mutations.columns

    n_types = signatures.shape[0]
    num_ref_sigs = signatures.shape[1]

    sel_sig_nums=list(range(0,num_ref_sigs))
    # sel_sig_nums=[0,4,26,44]
    # sel_sig_nums=range(0,44)+range(46,65)  # all 65 signatures except SBS40

    signature_columns_list = signatures.columns[sel_sig_nums].tolist()
    print("Analysing signatures:", signature_columns_list)
    # print("Analysing samples:", input_mutations.columns.tolist())

    output_weights = pd.DataFrame(index=samples, columns=signature_columns_list)
    output_mutations = pd.DataFrame(index=samples, columns=signature_columns_list)
    output_stat_info = pd.DataFrame(index=samples, columns=['Mutational burden', 'RSS', 'Chi2', 'R2', 'Cosine similarity', 'Correlation', 'Chebyshev similarity', 'L1 similarity', 'L2 similarity', 'L3 similarity', 'Jensen-Shannon similarity'])

    output_weights.index.name = output_mutations.index.name = output_stat_info.index.name = 'Sample'

    start_time = time.process_time()

    # align indexes to prevent signatures/samples mismatch
    # input_mutations.sort_index(inplace = True)
    # signatures.sort_index(inplace = True)
    input_mutations = input_mutations.reindex(signatures.index)

    # residuals and fitted values dataframes
    residuals_dataframe = pd.DataFrame(index=input_mutations.index, columns=input_mutations.columns)
    fitted_dataframe = pd.DataFrame(index=input_mutations.index, columns=input_mutations.columns)
    if args.bootstrap and args.bootstrap_method=="bootstrap_residuals":
        # load pre-created residuals and fitted values dataframes (needed for bootstrap)
        # so far only for SBS 192 context
        residuals_dataframe = pd.read_csv('%s/%s/output_%s_%s_residuals.csv' % (input_path, dataset_name, dataset_name, mutation_type), index_col=[0,1,2])
        fitted_dataframe = pd.read_csv('%s/%s/output_%s_%s_fitted_values.csv' % (input_path, dataset_name, dataset_name, mutation_type), index_col=[0,1,2])
        input_mutations = bootstrap_mutation_table(input_mutations, method=args.bootstrap_method, fitted=fitted_dataframe, residuals=residuals_dataframe)

    # print(input_mutations.index.tolist())
    # print(signatures.index.tolist())
    for sample in input_mutations.columns:
        print("Processing sample", sample)
        selected_mutations = input_mutations[sample].tolist()
        if not sum(selected_mutations)>0:
            warnings.warn("Zero mutations: skipping")
            continue

        initial_signatures = signatures.iloc[:,sel_sig_nums]

        if args.optimise_signatures:
            final_signatures = optimise_signatures(selected_mutations, initial_signatures, signatures,
                            strategy = args.optimisation_strategy, weak_threshold = args.weak_threshold,
                            strong_threshold = args.strong_threshold, verbose = args.verbose)
        else:
            # skip signature optimisation
            final_signatures = initial_signatures

        if not final_signatures.empty:
            normalised_weights, mutation_numbers, fitted, residuals, stat_info = perform_signature_attribution(selected_mutations, final_signatures, normalise_mutations=args.normalise_mutations, verbose=args.verbose)
            signatures_to_fill = final_signatures.columns
            output_weights.loc[sample, signatures_to_fill] = normalised_weights
            output_mutations.loc[sample, signatures_to_fill] = mutation_numbers
            output_stat_info.loc[sample] = stat_info
            residuals_dataframe[sample] = residuals
            fitted_dataframe[sample] = fitted

    end_time = time.process_time()
    process_elapsed_time = end_time - start_time
    print("Signature attribution took %.2f s (%.2f s per sample)" % (process_elapsed_time, process_elapsed_time/len(samples.tolist())))

    # replace NANs by zeros
    output_weights = output_weights.fillna(0.0)
    output_mutations = output_mutations.fillna(0.0)
    output_stat_info = output_stat_info.fillna(0.0)

    # write to files
    output_weights.to_csv(output_path + '/output_%s_%s_weights_table.csv' % (dataset_name, mutation_type))
    output_mutations.to_csv(output_path + '/output_%s_%s_mutations_table.csv' % (dataset_name, mutation_type))
    output_stat_info.to_csv(output_path + '/output_%s_%s_stat_info.csv' % (dataset_name, mutation_type))
    residuals_dataframe.to_csv(output_path + '/output_%s_%s_residuals.csv' % (dataset_name, mutation_type))
    fitted_dataframe.to_csv(output_path + '/output_%s_%s_fitted_values.csv' % (dataset_name, mutation_type))
    if not args.bootstrap:
        residuals_dataframe.to_csv('%s/%s/output_%s_%s_residuals.csv' % (input_path, dataset_name, dataset_name, mutation_type))
        fitted_dataframe.to_csv('%s/%s/output_%s_%s_fitted_values.csv' % (input_path, dataset_name, dataset_name, mutation_type))
