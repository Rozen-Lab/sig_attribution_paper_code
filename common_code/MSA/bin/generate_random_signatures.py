from optparse import OptionParser
import os, copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from common_methods import calculate_similarity
from itertools import combinations

def make_folder_if_not_exists(folder):
    if not os.path.exists(folder):
        try:
            os.makedirs(folder)
        except:
            warnings.warn("Could not create a folder ", folder)
#
def normalise_mutations(input_data):
    # normalise the number of mutations by the total number of mutations
    data = copy.deepcopy(input_data)
    total_mutational_burden = data.sum()
    data = data/total_mutational_burden
    return data

def plot_similarities(input_signatures, method='Cosine', savepath = './similarities.pdf'):
    make_folder_if_not_exists(savepath.rsplit('/',1)[0])

    data = copy.deepcopy(input_signatures)
    signatures_list = data.columns

    similarities = []

    for combination in combinations(signatures_list, 2):
        signature_A, signature_B = combination
        if metric in ['L1', 'L2', 'L3', 'Manhattan', 'manhattan', 'Chebyshev', 'chebyshev']:
            similarity = calculate_similarity(data[signature_A], data[signature_B], method, normalise=True)
        else:
            similarity = calculate_similarity(data[signature_A], data[signature_B], method)
        similarities.append(similarity)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    n, bins, patches = plt.hist(similarities, density=True, bins=20, color='blue', histtype='step')
    # fit a Gaussian to input distribution
    (mu, sigma) = stats.norm.fit(similarities)
    y = stats.norm.pdf( bins, mu, sigma)
    text = 'Fitted Gaussian:\n $\\mu$=%.2f, $\\sigma=$%.2f' %  (mu, sigma)

    ax.plot(bins, y, 'r--', linewidth=2)
    # remove indentation from latex-based text
    ax.text(0.65, 0.9, text, fontsize=9, transform=ax.axes.transAxes)

    ax.set_ylabel("Percentage", fontsize=12)
    ax.set_xlabel("Similarity", fontsize=12)
    ax.set_title(method + ' similarities (N = %i)' % len(signatures_list), fontsize=14, pad=10)
    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-t", "--mutation_type", dest="mutation_type", default='',
                      help="set mutation type (SBS, DBS, ID, SV, CNV)")
    parser.add_option("-c", "--context", dest="context", default=96, type='int',
                      help="set SBS context (96, 192, 1536)")
    parser.add_option("-s", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                      help="set path to PCAWG signature tables (for reference)")
    parser.add_option("-o", "--output_path", dest="output_path", default='signature_tables/',
                      help="set path to save output simulated signatures")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      help="verbosity flag for debugging (lots of output)")
    parser.add_option("-n", "--number", dest="number_of_signatures", default=100, type='int',
                      help="set the number of signatures to generate (100 by default)")
    parser.add_option("-z", "--noise", dest="add_noise", action="store_true",
                      help="Add Gaussian white noise")
    parser.add_option("-Z", "--noise_sigma", dest="noise_sigma", default=2, type='float',
                      help="Set the Gaussian standard deviation of the white noise (2 by default)")
    parser.add_option("-S", "--sparse", dest="sparse", action="store_true",
                      help="Generate sparse signatures (pick a level-1 category at random)")
    (options, args) = parser.parse_args()

    mutation_type = options.mutation_type
    context = options.context
    signature_tables_path = options.signature_tables_path
    output_path = options.output_path
    number_of_signatures = options.number_of_signatures

    make_folder_if_not_exists(output_path)

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (or '-t DBS', '-t ID').")
    elif mutation_type not in ['SBS','DBS','ID','SV','CNV']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID, SV, CNV" % mutation_type)

    if mutation_type=='SBS':
        if context==96:
            output_filename = '%s/sigRandom_%s_signatures.csv' % (output_path, mutation_type)
            reference_signatures = pd.read_csv('%s/sigProfiler_%s_signatures.csv' %
                                            (signature_tables_path, mutation_type), index_col=[0,1])
        elif context in [192, 288]:
            output_filename = '%s/sigRandom_SBS_%i_signatures.csv' % (output_path, context)
            reference_signatures = pd.read_csv('%s/sigProfiler_SBS_%i_signatures.csv' %
                                            (signature_tables_path, context), index_col=[0,1,2])
        elif context==1536:
            output_filename = '%s/sigRandom_SBS_%i_signatures.csv' % (output_path, context)
            reference_signatures = pd.read_csv('%s/sigProfiler_SBS_%i_signatures.csv' %
                                            (signature_tables_path, context), index_col=0)
        else:
            raise ValueError("Context %i is not supported." % context)
    else:
        output_filename = '%s/sigRandom_%s_signatures.csv' % (output_path, mutation_type)
        reference_signatures = pd.read_csv('%s/sigProfiler_%s_signatures.csv' %
                                            (signature_tables_path, mutation_type), index_col=0)

    signatures_range = range(0,number_of_signatures)
    generated_signatures = pd.DataFrame(0, index=reference_signatures.index, columns=signatures_range)

    for i in signatures_range:
        # generated_signatures[i] = np.random.normal(0, 10, len(generated_signatures[i].values))
        generated_signatures[i] = np.random.poisson(lam=3, size=len(generated_signatures[i].values))
        generated_signatures.loc[generated_signatures[i]<0, i] = 0
        # print(generated_signatures[i].values)

        if options.add_noise:
            noise_term = np.random.normal(0, options.noise_sigma, len(generated_signatures[i].values))
            generated_signatures[i] += noise_term
            # make sure there are no negative mutation counts
            generated_signatures.loc[generated_signatures[i]<0, i] = 0

    if options.sparse:
        if not mutation_type=='SBS':
            warnings.warn("Sparse signatures are not supported for non-SBS signatures yet. Generating homogeneous signatures.")
        else:
            if context==96:
                main_category_index_level = 0
            elif context in [192, 288]:
                main_category_index_level = 1
            else:
                raise ValueError("Context %i is not supported." % context)
        main_categories = generated_signatures.index.unique(main_category_index_level)
        # print(main_categories)
        for i in signatures_range:
            # pick a random main category (e.g. C>A, T>C, etc)
            category = main_categories[np.random.randint(0, len(main_categories))]
            # print('Picked category:', category)
            # set all other categories to zero (to generate sparseness)
            generated_signatures.loc[generated_signatures.index.get_level_values(main_category_index_level) != category, i] = 0

    # normalise the signatures to one
    for i in signatures_range:
        generated_signatures[i] = normalise_mutations(generated_signatures[i])
        # print(generated_signatures[i].values)

    # # save dataframes
    generated_signatures.to_csv(output_filename)
    metrics = ['Jensen-Shannon', 'L1', 'L2', 'L3', 'Chebyshev', 'Correlation', 'Cosine']
    for metric in metrics:
        plot_similarities(generated_signatures, metric, savepath = output_path + '/similarities/' + metric + '.pdf')
