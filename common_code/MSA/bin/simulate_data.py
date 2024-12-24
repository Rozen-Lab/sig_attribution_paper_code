import argparse
import os
import warnings
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from common_methods import make_folder_if_not_exists
from math import ceil

# signatures to generate if not bootstrapping and input signature activities table (-B option)
signatures_to_generate = {
    # Dictionary with signature burdens to be generated.
    # Each signature is assigned with a list of type [mu, sigma], where mu is the mean
    # and sigma is the standard deviation of the normal distribution used for simulating
    # the mutational burden attributed to the signature.
    # Signature names (strings) have to be present in the input signature tables.
    'SBS1':[150, 100],
    'SBS2':[150, 130],
    'SBS5':[0, 50],
    'SBS13':[150, 130],
    'SBS18':[100, 100],
    'SBS22':[0, 50],
    'SBS40':[200, 300],
    'DBS1':[140,270],
    'DBS2':[2000,1000],
    'ID1':[140,270],
    'ID2':[2000,1000]
}

# additional signatures to inject if -I option is specified (useful with -B option)
signatures_to_inject = {
    # Dictionary with signature burdens to be injected in addition to modelled/bootstrapped signatures.
    # Each signature is assigned with a list of type [mu, sigma, p], where mu is the mean
    # and sigma is the standard deviation of the normal distribution used for simulating
    # the mutational burden attributed to the signature, and p is the probability of injection to a given sample.
    # Signature burden means and standard deviations are specified as ratios of mutation burdens of the samples to be injected in.
    # Signature names (strings) have to be present in the input signature tables.
    'SBS4':[0.05, 0.005, 0.1]
}

def plot_mutational_burden(mutational_burden, mu=None, sigma=None, title='Total', x_label='Mutation count', y_label='Normalised number of samples', savepath = './burden.pdf'):
    """
    Plot a histogram of the input mutational burden, plot the Gaussian on top if
    the mean (mu) and standard deviation (sigma) are provided, otherwise fit a
    Gaussian to the histogram.

    Parameters:
    ----------
    mutational_burden: array_like
        list of numbers (ints or floats)

    mu: float (Optional)
        Mean of the normal distribution to plot

    sigma: float (Optional)
        Standard deviation of the normal distribution to plot

    Returns:
    -------
    Nothing. Creates and saves a plot in specified path.
    -------
    """
    make_folder_if_not_exists(savepath.rsplit('/',1)[0])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, patches = plt.hist(mutational_burden, density=True, bins=20, color='blue', histtype='step')
    if mu and sigma:
        # create a Gaussian with given mu and sigma
        y = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) )
        text = '$\\mu$=%.2f, $\\sigma=$%.2f' %  (mu, sigma)
    else:
        # fit a Gaussian to input distribution
        (mu, sigma) = stats.norm.fit(mutational_burden)
        y = stats.norm.pdf( bins, mu, sigma)
        text = 'Fitted Gaussian:\n $\\mu$=%.2f, $\\sigma=$%.2f' %  (mu, sigma)

    number_of_positive_samples = sum(x > 0 for x in mutational_burden)
    text += '\n %i/%i samples' % (number_of_positive_samples, len(mutational_burden))

    ax.plot(bins, y, 'r--', linewidth=2)
    # remove indentation from latex-based text
    ax.text(0.65, 0.85, text, fontsize=9, transform=ax.axes.transAxes)

    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(title, fontsize=14, pad=10)
    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='',
                        help="set mutation type (SBS, DBS, ID, SV, CNV)")
    parser.add_argument("-c", "--context", dest="context", default=96, type=int,
                        help="set SBS context (96, 192, 288, 1536)")
    parser.add_argument("-s", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                        help="set path to signature tables")
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                        help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_argument("-o", "--output_path", dest="output_path", default='input_mutation_tables/',
                        help="set path to save output simulatied mutation tables")
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='SIM',
                        help="set the dataset name ('SIM' by default)")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="verbosity flag for debugging (lots of output)")
    parser.add_argument("-n", "--number_of_samples", dest="number_of_samples", default=-1, type=int,
                        help="set the number of samples to generate (-1 by default means 10 times number of input samples but at least 1000)")
    parser.add_argument("-z", "--noise", dest="add_noise", action="store_true",
                        help="Add noise of type specified by noise_type parameter")
    parser.add_argument("--noise_type", dest="noise_type", default='poisson',
                        help="Choose the type of noise: Poisson, Gaussian or negative_binomial (Poisson variation by default)")
    parser.add_argument("-Z", "--noise_sigma", dest="noise_sigma", default=5, type=float,
                        help="Set standard deviation (in percentage of sample mutation burden) of Gaussian noise if used (5 percent by default)")
    parser.add_argument("--zero_inflation_threshold", dest="zero_inflation_threshold", default=0.01, type=float,
                        help="set the relative threshold below which all simulated signature activities are set to zero (0.01 by default)")
    parser.add_argument("-r", "--random", dest="random_signatures", action="store_true",
                        help="Use randomly generated signatures instead of PCAWG reference ones")
    parser.add_argument("-N", "--number_of_random_sigs", dest="number_of_random_sigs", default=5, type=int,
                        help="Number of random signatures to consider")
    parser.add_argument("-B", "--bootstrap_input_signature_activities_table", dest="bootstrap_input_signature_activities_table", action="store_true",
                        help="Bootstrap (reshuffle with replacement) the input signature activities table instead of using signatures_to_generate dictionary")
    parser.add_argument("-i", "--input_table", dest="input_table", default='',
                        help="set path to input signature activities table to bootstrap")
    parser.add_argument("-I", "--inject_signatures", dest="inject_signatures", action="store_true",
                        help="inject additional signatures as specified in signatures_to_inject dictionary in the script")
    args = parser.parse_args()

    mutation_type = args.mutation_type
    dataset_name = args.dataset_name
    context = args.context
    signature_tables_path = args.signature_tables_path
    signatures_prefix = args.signatures_prefix
    output_path = args.output_path + '/' + dataset_name
    number_of_samples = args.number_of_samples
    random_signatures = args.random_signatures
    number_of_random_sigs = args.number_of_random_sigs

    # calculate number of simulated samples by default
    if number_of_samples==-1:
        if not args.input_table:
            warnings.warn("The input signature activities table not provided with -i option. Attempting to run 1000 simulations.")
            number_of_samples = 1000
        else:
            # calculating the number of samples in the input signature activities table
            number_of_input_lines = sum(1 for line in open(args.input_table))
            # rounded factor of ten but not less than 1000
            number_of_samples = 1000*ceil(number_of_input_lines/100.0)

    make_folder_if_not_exists(output_path)

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (or '-t DBS', '-t ID').")
    elif mutation_type not in ['SBS','DBS','ID','SV','CNV']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID, SV, CNV" % mutation_type)

    if mutation_type=='SBS':
        output_filename = '%s/WGS_%s.%i.csv' % (output_path, dataset_name, context)
        if context==96:
            reference_signatures = pd.read_csv('%s/%s_%s_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type), index_col=[0,1])
        elif context in [192, 288]:
            reference_signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type, context), index_col=[0,1,2])
        elif context==1536:
            reference_signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type, context), index_col=0)
        else:
            raise ValueError("Context %i is not supported." % context)
    elif mutation_type=='DBS':
        output_filename = '%s/WGS_%s.dinucs.csv' % (output_path, dataset_name)
        reference_signatures = pd.read_csv('%s/%s_%s_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
    elif mutation_type=='ID':
        output_filename = '%s/WGS_%s.indels.csv' % (output_path, dataset_name)
        reference_signatures = pd.read_csv('%s/%s_%s_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
    else: # SV and CNV
        output_filename = '%s/WGS_%s.%s.csv' % (output_path, dataset_name, mutation_type)
        reference_signatures = pd.read_csv('%s/%s_%s_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type), index_col=0)

    for signature in reference_signatures.columns:
        if not np.isclose(reference_signatures.sum()[signature], 1, rtol=1e-2):
            raise ValueError("Probabilities for signature %s do not add up to 1: %.3f" % (signature, reference_signatures.sum()[signature]))

    generated_signature_burdens = {}
    if args.bootstrap_input_signature_activities_table:
        if not args.input_table:
            parser.error("Please provide the input signature activities table for bootstrap with -i option.")
        input_table = pd.read_csv(args.input_table, index_col=0, sep=None)
        bootstrapped_table = input_table.sample(n=number_of_samples, replace=True)
        for signature in input_table.columns:
            print('Generating signature burden for', signature)
            generated_signature_burdens[signature] = bootstrapped_table[signature].to_numpy()
            generated_signature_burdens[signature] = np.nan_to_num(generated_signature_burdens[signature])
            plot_mutational_burden(generated_signature_burdens[signature], title=signature,
            savepath='%s/%s_plots/generated_burden_%s.pdf' % (output_path, dataset_name, signature))
    else:
        if random_signatures:
            random_signature_indexes = random.sample(reference_signatures.columns.to_list(), number_of_random_sigs)
            for random_sig in random_signature_indexes:
                mu = np.random.normal(2000, 1000)
                sigma = np.random.normal(500, 100)
                mu = mu if mu>=0 else 0
                sigma = sigma if sigma>=0 else 0
                # draw burdens from normal distribution
                generated_burdens = np.random.normal(mu, sigma, number_of_samples)
                # replace negative values of the Gaussian with zeros
                generated_burdens[generated_burdens<0] = 0
                # plot the generated burden
                plot_mutational_burden(generated_burdens, mu, sigma, random_sig,
                                    savepath='%s/%s_plots/generated_burden_%s.pdf' % (output_path, dataset_name, random_sig))
                generated_signature_burdens[random_sig] = generated_burdens
        else:
            # generate burdens for each signature from the 'signatures_to_generate' dictionary
            for signature in signatures_to_generate.keys():
                if mutation_type not in signature:
                    continue
                mu = signatures_to_generate[signature][0]
                sigma = signatures_to_generate[signature][1]
                # draw burdens from normal distribution
                generated_burdens = np.random.normal(mu, sigma, number_of_samples)
                # replace negative values of the Gaussian with zeros
                generated_burdens[generated_burdens<0] = 0
                # plot the generated burden
                plot_mutational_burden(generated_burdens, mu, sigma, signature,
                                    savepath='%s/%s_plots/generated_burden_%s.pdf' % (output_path, dataset_name, signature))
                generated_signature_burdens[signature] = generated_burdens

    samples_range = range(0,number_of_samples)

    generated_weights = pd.DataFrame(0, index=samples_range, columns=reference_signatures.columns)
    generated_mutations = pd.DataFrame(0, index=reference_signatures.index, columns=samples_range)

    # obtain generated mutational burdens
    generated_mutational_burdens = []
    for i in samples_range:
        mutational_burden = sum(signature_burden[i] for signature_burden in list(generated_signature_burdens.values()))
        generated_mutational_burdens.append(mutational_burden)
    generated_mutational_burdens = np.asarray(generated_mutational_burdens)
    generated_mutational_burdens = np.nan_to_num(generated_mutational_burdens)

    # optionally, inject additional signatures:
    if args.inject_signatures:
        for signature in signatures_to_inject.keys():
            if mutation_type not in signature:
                continue
            if not args.bootstrap_input_signature_activities_table and not random_signatures and signature in signatures_to_generate:
                raise ValueError("Signature %s can not be generated and injected simultaneously" % (signature))
            mu = signatures_to_inject[signature][0]*mutational_burden
            sigma = signatures_to_inject[signature][1]*mutational_burden
            p = signatures_to_inject[signature][2]
            # draw burdens from normal distribution, multiplied by binomial distribution for p probability of injection
            generated_burdens = np.random.normal(mu, sigma, number_of_samples) * np.random.binomial(1, p, number_of_samples)
            # draw burdens from exponential instead:
            # generated_burdens = np.random.exponential(500, number_of_samples) * np.random.binomial(1, p, number_of_samples)
            # replace negative values of the Gaussian with zeros
            generated_burdens[generated_burdens<0] = 0
            # plot the generated burden
            plot_mutational_burden(generated_burdens, mu, sigma, signature,
                                savepath='%s/%s_plots/generated_burden_%s.pdf' % (output_path, dataset_name, signature))
            if signature in generated_signature_burdens.keys():
                generated_signature_burdens[signature] += generated_burdens
            else:
                generated_signature_burdens[signature] = generated_burdens

    # loop to fill mutation tables:
    for i in samples_range:
        weights = {}
        for signature in generated_signature_burdens.keys():
            if mutational_burden!=0:
                weights[signature] = generated_signature_burdens[signature][i]/generated_mutational_burdens[i]
                # replace signature activities of <1% with 0
                if weights[signature] < args.zero_inflation_threshold:
                    weights[signature] = 0
            else:
                weights[signature] = 0
            generated_weights.loc[i,signature] = weights[signature]
            generated_mutations[i] += reference_signatures[signature]*weights[signature]

        # multiply by mutational burden
        generated_mutations[i] = generated_mutational_burdens[i]*generated_mutations[i]

        # optionally, add noise:
        if args.add_noise:
            if args.noise_type=="gaussian" or args.noise_type=="Gaussian" or args.noise_type=="normal" or args.noise_type=="Normal":
                # # absolute stdev implementation
                # noise_term = np.random.normal(0, args.noise_sigma, len(generated_mutations[i]))
                # generated_mutations[i] += noise_term
                # # relative stdev implementation
                noise_stdev = generated_mutations[i]*args.noise_sigma/100
                noisy_data = np.random.normal(generated_mutations[i], noise_stdev)
                generated_mutations[i] = noisy_data
            elif args.noise_type=="poisson" or args.noise_type=="Poisson":
                for category_index in range(len(generated_mutations[i])):
                    generated_mutations.iloc[category_index, i] = np.random.poisson(generated_mutations.iloc[category_index, i])
                    #generated_mutations[i][category_index] = np.random.poisson(generated_mutations[i][category_index])
            elif args.noise_type=="negative_binomial" or args.noise_type=="Negative_binomial":
                noise_term = np.random.negative_binomial(2, 0.5, len(reference_signatures[signature].values))
                generated_mutations[i] += noise_term
            # make sure there are no negative mutation counts
            generated_mutations.loc[generated_mutations[i]<0, i] = 0

        # treating nans with zero (happens when mutation burden is zero in older pandas versions)
        generated_mutations[i].fillna(0, inplace=True)

        # rounding and converting to integer counts
        generated_mutations[i] = round(generated_mutations[i],0)
        generated_mutations[i] = generated_mutations[i].astype(int)

    # treating nans with zero for generated weights
    generated_weights.fillna(0, inplace=True)
    # plot total mutational burden distribution:
    plot_mutational_burden(generated_mutational_burdens, savepath='%s/%s_plots/generated_burden_%s_total.pdf' % (output_path, dataset_name, mutation_type) )

    # plot relative signature weights:
    for signature in generated_weights.columns:
        print('Plotting weights for signature', signature)
        plot_mutational_burden(np.nan_to_num(generated_weights[signature].to_numpy()), title=signature, x_label = 'Relative weight',
            savepath='%s/%s_plots/generated_weight_%s.pdf' % (output_path, dataset_name, signature))

    # save dataframes
    generated_mutations.to_csv(output_filename)
    generated_weights.to_csv(output_filename.replace('.csv','.weights.csv'))
