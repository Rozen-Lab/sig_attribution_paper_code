from argparse import ArgumentParser
import copy
import pandas as pd
import numpy as np
from scipy.spatial import distance
from common_methods import calculate_similarity
import matplotlib.pyplot as plt

def plot_array_as_histogram(arrays, labels, title, savepath='./hist.pdf'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # bins = [0.025*k for k in range(42)]
    bins = [0.01*k for k in range(101)]
    for array, label in zip(arrays, labels):
        plt.hist(array, label = label, histtype='step', bins = bins)

    legend = plt.legend(loc="upper left")
    plt.gca().add_artist(legend)

    ax.set_title(title, fontsize=14, pad=10)
    # ax.set_xlabel('Attribution ($\%$)', fontsize=12)
    # ax.set_ylabel('Number of attributions', fontsize=12)
    # ax.set_xticks([0.1*k for k in range(11)])

    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_matrix", dest="input_matrix", default='input_mutation_tables/SIM_ESCC/WGS_SIM_ESCC.192.csv',
                        help="set path to the input matrix to calculate similarities from")
    parser.add_argument("-o", "--output_path", dest="output_path", default='./hist.pdf',
                        help="set path to save the plot")
    parser.add_argument("-n", "--number", dest="number_of_samples", default=-1, type=int,
                        help="limit the number of samples to analyse (all by default)")
    parser.add_argument("-N", "--normalise", dest="normalise", action="store_true",
                        help="Normalise vectors by their L2 norm")

    args = parser.parse_args()
    input_matrix = pd.read_csv(args.input_matrix, index_col=[0,1,2])

    # limit the number of samples to analyse (if specified by -n option)
    if args.number_of_samples!=-1:
        input_matrix = input_matrix.iloc[:,0:args.number_of_samples]

    number_of_samples = len(input_matrix.columns)
    L1_sims = []
    L2_sims = []
    L3_sims = []
    Chebyshev_sims = []
    corr_sims = []
    cosine_sims = []
    JS_sims = []
    L2_norm_sims = []

    for i in input_matrix.columns:
        for j in input_matrix.columns:
            if i==j:
                continue
            # L2_norm_sims.append(calculate_similarity(input_matrix[i],input_matrix[j], metric='eucl_norm_squared', normalise = args.normalise))
            L1_sims.append(calculate_similarity(input_matrix[i],input_matrix[j], metric='L1', normalise = args.normalise))
            L2_sims.append(calculate_similarity(input_matrix[i],input_matrix[j], metric='L2', normalise = args.normalise))
            L3_sims.append(calculate_similarity(input_matrix[i],input_matrix[j], metric='L3', normalise = args.normalise))
            Chebyshev_sims.append(calculate_similarity(input_matrix[i],input_matrix[j], metric='Chebyshev', normalise = args.normalise))
            corr_sims.append(calculate_similarity(input_matrix[i],input_matrix[j], metric='correlation', normalise = args.normalise))
            cosine_sims.append(calculate_similarity(input_matrix[i],input_matrix[j], normalise = args.normalise))
            JS_sims.append(calculate_similarity(input_matrix[i],input_matrix[j], metric='jensen-shannon', normalise = args.normalise))

    if args.normalise:
        plot_array_as_histogram([corr_sims, cosine_sims, JS_sims, L1_sims, L2_sims, L3_sims, Chebyshev_sims], ['Correlation', 'Cosine', 'JS', 'L1', 'L2', 'L3', 'Chebyshev'], 'Similarities', savepath = args.output_path)
    else:
        plot_array_as_histogram([corr_sims, cosine_sims, JS_sims], ['Correlation', 'Cosine', 'JS'], 'Similarities', savepath = args.output_path)
