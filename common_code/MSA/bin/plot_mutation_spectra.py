from argparse import ArgumentParser
import pandas as pd
from common_methods import make_folder_if_not_exists
from common_methods import categorise_mutations, normalise_mutations
from common_methods import mutational_burden, calculate_errors, merge_plots
from common_methods import plot_array_as_histogram
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
from scipy import stats

mutation_category_labels = {
    #SBS
    'C>A':'C$>$A',
    'C>G':'C$>$G',
    'C>T':'C$>$T',
    'T>A':'T$>$A',
    'T>C':'T$>$C',
    'T>G':'T$>$G',
    #DBS
    'AC>':'AC$>$NN',
    'AT>':'AT$>$NN',
    'CC>':'CC$>$NN',
    'CG>':'CG$>$NN',
    'CT>':'CT$>$NN',
    'GC>':'GC$>$NN',
    'TA>':'TA$>$NN',
    'TC>':'TC$>$NN',
    'TG>':'TG$>$NN',
    'TT>':'TT$>$NN',
    #ID
    'DEL_C_1':'C',
    'DEL_T_1':'T',
    'INS_C_1':'C',
    'INS_T_1':'T',
    'DEL_repeats_2':'2',
    'DEL_repeats_3':'3',
    'DEL_repeats_4':'4',
    'DEL_repeats_5+':'5+',
    'INS_repeats_2':'2',
    'INS_repeats_3':'3',
    'INS_repeats_4':'4',
    'INS_repeats_5+':'5+',
    'DEL_MH_2':'2',
    'DEL_MH_3':'3',
    'DEL_MH_4':'4',
    'DEL_MH_5+':'5+',
    #SV
    'clustered_del':'',
    'clustered_tds':'',
    'clustered_inv':'',
    'clustered_trans':'',
    'non-clustered_del':'',
    'non-clustered_tds':'',
    'non-clustered_inv':'',
    'non-clustered_trans':'',
    #CNV
    '0:homdel':'',
    '1:LOH':'',
    '2:LOH':'',
    '3-4:LOH':'',
    '5-8:LOH':'',
    '9+:LOH':'',
    '2:het':'',
    '3-4:het':'',
    '5-8:het':'',
    '9+:het':'',
}

mutation_categories = {
    'SBS':['C>A','C>G','C>T','T>A','T>C','T>G'],
    'DBS':['AC>','AT>','CC>','CG>','CT>','GC>','TA>','TC>','TG>','TT>'],
    'ID':['DEL_C_1', 'DEL_T_1', 'INS_C_1', 'INS_T_1', 'DEL_repeats_2', 'DEL_repeats_3', 'DEL_repeats_4', 'DEL_repeats_5+',
    'INS_repeats_2', 'INS_repeats_3', 'INS_repeats_4', 'INS_repeats_5+', 'DEL_MH_2', 'DEL_MH_3', 'DEL_MH_4', 'DEL_MH_5+'],
    'SV':['clustered_del', 'clustered_tds', 'clustered_inv', 'clustered_trans', 'non-clustered_del', 'non-clustered_tds', 'non-clustered_inv', 'non-clustered_trans'],
    'CNV':['0:homdel', '1:LOH', '2:LOH', '3-4:LOH', '5-8:LOH', '9+:LOH', '2:het', '3-4:het', '5-8:het', '9+:het']
}

mutation_colours = {
    'SBS':['skyblue','black','firebrick','gray','lightgreen','salmon'],
    'DBS':['skyblue','royalblue','lightgreen','green','salmon','firebrick','navajowhite','darkorange','plum','blueviolet'],
    'ID':['navajowhite','darkorange','lightgreen','green','lightcoral','coral','red','firebrick','lightskyblue','skyblue',
    'deepskyblue','royalblue','lavender','plum','mediumpurple','blueviolet'],
    'SV':['navajowhite','darkorange','lightgreen','green','lightcoral','coral','red','firebrick'],
    'CNV':['navajowhite','darkorange','lightgreen','green','lightcoral','coral','red','firebrick','lightskyblue','skyblue'],
}

def lighten_colour(colour, amount=0.5):
    """
    Lightens the given colour by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib colour string, hex string, or RGB tuple.

    Examples:
    >> lighten_colour('g', 0.3)
    >> lighten_colour('#F034A3', 0.6)
    >> lighten_colour((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[colour]
    except:
        c = colour
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def make_condensed_spectrum_plot(data, mutation_type='SBS', title='', y_label='', show_errors=False, relative=False, strand_bias=False, add_non_transcribed=False, savepath="./test.pdf", text='', additional_text=''):
    make_folder_if_not_exists(savepath.rsplit('/',1)[0])

    categories = mutation_categories[mutation_type]
    colours = mutation_colours[mutation_type]

    if len(categories)!=len(colours):
        raise ValueError("Length of categories/colours lists do not match:", categories, colours)

    mutations = categorise_mutations(data, categories, mutation_type=mutation_type, condensed=True, strand_bias=strand_bias, add_non_transcribed=add_non_transcribed)
    # latex handling
    if mutation_type=='ID':
        mutations.index = [category.replace('_','-') for category in list(mutations.index)]
    else:
        mutations.index = [mutation_category_labels[category] for category in list(mutations.index)]

    # annotate strand bias before normalisation, otherwise no significance can be found using fractions
    # strand_bias_present = False
    # if strand_bias:
    #     mutations, strand_bias_present, _ = annotate_strand_bias(mutations)

    if show_errors:
        errors = calculate_errors(mutations, normalise = relative)

    # normalise if needed
    if relative:
        mutations = normalise_mutations(mutations)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=14, pad=10)

    # only show integer y axis labels, unless relative or strand bias (horizontal plot)
    if not relative and not strand_bias:
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    if strand_bias:
        if add_non_transcribed:
            strand_colours = ['green', 'red', 'blue']
        else:
            strand_colours = ['red', 'blue']
        mutations.iloc[::-1, ::-1].plot.barh(ax = ax, xerr = errors if show_errors else None, edgecolor = "none", color = strand_colours)
        # reverse the order in legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], frameon=False)
        plt.xlabel(y_label, fontsize=12)
    else:
        mutations.plot.bar(ax = ax, yerr = errors if show_errors else None, edgecolor = "none", color = colours, legend=False)
        plt.ylabel(y_label, fontsize=12)

    # if strand_bias_present:
    #     plt.text(0.657, 1.01, '$\\star$ indicates significant strand bias', fontsize=10, transform=ax.transAxes)
    if text or additional_text:
        # add some text
        plt.text(1.1, 0.55, text, fontsize=12, transform=ax.transAxes)
        plt.text(1.1, 0.24, additional_text, fontsize=12, transform=ax.transAxes)

        fig.set_tight_layout({'rect':[0.01, 0.03, 0.8, 0.95]})
    else:
        fig.set_tight_layout(True)

    plt.savefig(savepath, transparent=True)
    plt.close()

def make_spectrum_plot(data, mutation_type='SBS', title='', y_label='', show_errors=False, relative=False, strand_bias = False, add_non_transcribed=False, condensed=False, savepath="./test.pdf", text='', additional_text=''):
    if condensed:
        make_condensed_spectrum_plot(data, mutation_type=mutation_type, title=title, y_label=y_label, show_errors = show_errors, relative=relative, strand_bias = strand_bias, add_non_transcribed=add_non_transcribed, savepath=savepath, text=text, additional_text=additional_text)
        return
    make_folder_if_not_exists(savepath.rsplit('/',1)[0])

    categories = mutation_categories[mutation_type]
    colours = mutation_colours[mutation_type]

    if len(categories)!=len(colours):
        raise ValueError("Length of categories/colours lists do not match:", categories, colours)

    mutation_dict = categorise_mutations(data, categories, mutation_type=mutation_type, strand_bias=strand_bias, add_non_transcribed=add_non_transcribed)

    # # annotate strand bias before normalisation, otherwise no significance can be found using fractions
    # strand_bias_present = False
    # if strand_bias:
    #     for category in categories:
    #         mutation_dict[category], strand_bias_present_in_category, _ = annotate_strand_bias(mutation_dict[category])
    #         if strand_bias_present_in_category: strand_bias_present = True

    if show_errors:
        errors = calculate_errors(mutation_dict, normalise = relative)
    # normalise mutation dictionary if needed
    if relative:
        mutation_dict = normalise_mutations(mutation_dict)

    f, axes = plt.subplots(1, len(categories), sharey=True, figsize=(20, 3))
    f.suptitle(title, fontsize=14, y=0.88)

    for category, axis, colour in zip(categories, axes, colours):
        if strand_bias:
            if add_non_transcribed:
                strand_colours = [colour, lighten_colour(colour), 'green']
            else:
                strand_colours = [colour, lighten_colour(colour)]
            mutation_dict[category].plot.bar(ax=axis, yerr = errors[category] if show_errors else None,
                    color=strand_colours, ecolor="darkgrey", edgecolor = "none").legend(frameon=False)
        else:
            mutation_dict[category].plot.bar(ax=axis, color=colour, edgecolor="none",
                                     ecolor="darkgrey", yerr=errors[category] if show_errors else None)

        # latex handling
        category_to_plot = mutation_category_labels[category]

        if colour in ['blue','black','firebrick','green','royalblue','blueviolet']:
            label_text_colour = 'white'
        else:
            label_text_colour = 'black'

        axis.set_title(category_to_plot, fontsize=10, color = label_text_colour)
        axis.xaxis.label.set_visible(False)

        # Create a Rectangle patch
        rect = patches.Rectangle((0,1),1,0.2,linewidth=1,clip_on=False,edgecolor=colour,facecolor=colour,transform=axis.transAxes)
        # Add the patch to the Axes
        axis.add_patch(rect)

        # horizontal x-axis labels for indels
        if mutation_type=='ID':
            axis.xaxis.set_tick_params(rotation=0)


    axes[0].set_ylabel(y_label, fontsize=12)
    if not relative:
        axes[0].yaxis.set_major_locator(MaxNLocator(integer=True))

    # some additional text
    # if strand_bias_present:
        # plt.text(0.36, 1.155, '$\\star$ indicates significant strand bias', fontsize=10, transform=axes[-1].transAxes)
    if text or additional_text:
        # add some text
        plt.text(1.1, 0.55, text, fontsize=12, transform=axes[-1].transAxes)
        plt.text(1.1, -0.24, additional_text, fontsize=12, transform=axes[-1].transAxes)

        f.set_tight_layout({'rect':[0.01, 0.03, 0.9, 0.95]})
    else:
        f.set_tight_layout(True)

    # handling the ticks
    plt.minorticks_off()
    # for ax in axes:
    #     ax.tick_params(axis=u'both', which=u'both',length=0)
    #     ax.spines['top'].set_visible(False)
    #     ax.spines['right'].set_visible(False)
    #     ax.spines['bottom'].set_visible(False)
    #     ax.spines['left'].set_visible(False)

    # additional subtitles for indels
    if mutation_type=='ID':
        # leave more space at the top to accomodate the additional titles
        f.subplots_adjust(top=0.72)

        # save the axes bounding boxes for later use
        ext = []
        for i in range(len(categories)):
            ext.append([axes[i].get_window_extent().x0, axes[i].get_window_extent().width ])

        # from the axes bounding boxes calculate the optimal position of the column spanning title
        inv = f.transFigure.inverted()
        first_del_width = ext[0][0]+(ext[1][0]+ext[1][1]-ext[0][0])/2.
        first_del_center = inv.transform( (first_del_width, 1) )
        first_ins_width = ext[2][0]+(ext[3][0]+ext[3][1]-ext[2][0])/2.
        first_ins_center = inv.transform( (first_ins_width, 1) )
        second_del_width = ext[4][0]+(ext[7][0]+ext[7][1]-ext[4][0])/2.
        second_del_center = inv.transform( (second_del_width, 1) )
        second_ins_width = ext[8][0]+(ext[11][0]+ext[11][1]-ext[8][0])/2.
        second_ins_center = inv.transform( (second_ins_width, 1) )
        MH_del_width = ext[12][0]+(ext[15][0]+ext[15][1]-ext[12][0])/2.
        MH_del_center = inv.transform( (MH_del_width, 1) )

        # set column spanning title
        # the first two arguments to figtext are x and y coordinates in the figure system (0 to 1)
        plt.figtext(first_del_center[0],0.88,"1 bp deletion", va="center", ha="center", size=13)
        plt.figtext(first_del_center[0],0.05,"Homopolymer length", va="center", ha="center", size=13)
        plt.figtext(first_ins_center[0],0.88,"1 bp insertion", va="center", ha="center", size=13)
        plt.figtext(first_ins_center[0],0.05,"Homopolymer length", va="center", ha="center", size=13)
        plt.figtext(second_del_center[0],0.88,"$>$1bp deletions at repeats \n (Deletion length)", va="center", ha="center", size=13)
        plt.figtext(second_del_center[0],0.05,"Number of repeat units", va="center", ha="center", size=13)
        plt.figtext(second_ins_center[0],0.88,"$>$1bp insertions at repeats \n (Insertion length)", va="center", ha="center", size=13)
        plt.figtext(second_ins_center[0],0.05,"Number of repeat units", va="center", ha="center", size=13)
        plt.figtext(MH_del_center[0],0.88,"Deletions with microhomology \n (Deletion length)", va="center", ha="center", size=13)
        plt.figtext(MH_del_center[0],0.05,"Microhomology length", va="center", ha="center", size=13)

    plt.savefig(savepath, transparent=True)
    plt.close()

def plot_residuals(input_array, title = '', savepath = './test_gaus.pdf'):
    make_folder_if_not_exists(savepath.rsplit('/',1)[0])

    fig = plt.figure()
    ax = fig.add_subplot(111)

    n, bins, patches = plt.hist(input_array, normed=True, bins=50, color='blue', histtype='step')
    # fit a Gaussian to input distribution
    (mu, sigma) = stats.norm.fit(input_array)
    y = stats.norm.pdf( bins, mu, sigma)
    text = 'Fitted Gaussian:\n $\\mu$=%.2f, $\\sigma=$%.2f' %  (mu, sigma)

    ax.plot(bins, y, 'r--', linewidth=2)
    # remove indentation from latex-based text
    ax.text(0.65, 0.9, text, fontsize=9, transform=ax.axes.transAxes)

    ax.set_ylabel("Normalised N of SNVs", fontsize=12)
    ax.set_xlabel("Residuals", fontsize=12)
    ax.set_title(title, fontsize=14, pad=10)
    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_folder", dest="input_folder", default='input_mutation_tables/',
                        help="set path to datasets with input mutation tables")
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='',
                        help="set mutation type (SBS, DBS, ID, SV, CNV)")
    parser.add_argument("-c", "--context", dest="context", default=96, type=int,
                        help="set SBS context (96, 192, 288, 1536)")
    parser.add_argument("-s", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                        help="set path to signature tables")
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                        help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_argument("-o", "--output_folder", dest="output_folder", default='plots/',
                        help="set path to save plots")
    parser.add_argument("-r", "--relative", dest="relative", action="store_true",
                        help="make plots with relative (not absolute) mutation counts")
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='ESCC',
                        help="set the dataset name (PCAWG, TCE)")
    parser.add_argument("-e", "--show_errors", dest="show_errors", action="store_true",
                        help="show Poisson errors")
    parser.add_argument("-b", "--strand_bias", dest="strand_bias", action="store_true",
                        help="plot spectra with transcribed and untranscribed strands using the 192 mutation types (SBS)")
    parser.add_argument("-n", "--non_transcribed", dest="non_transcribed", action="store_true",
                        help="Add non transcribed region in plots, in addition to strands")
    parser.add_argument("-C", "--condensed", dest="condensed", action="store_true",
                        help="Make condensed plots with fewer categories")
    parser.add_argument("--number", dest="number_of_spectra", default=-1, type=int,
                        help="limit the number of spectra to plot (all by default)")
    parser.add_argument("--samples", nargs='+', dest="samples", default = [],
                        help="Make plots for the specified list of samples, e.g. -S PD37458a PD37459a")
    parser.add_argument("-f", "--plot_fitted_spectra", dest="plot_fitted_spectra", action="store_true",
                        help="Make plots for the fitted spectra (after attribution)")
    parser.add_argument("-R", "--plot_residuals", dest="plot_residuals", action="store_true",
                        help="Make residuals plots (difference between input and fitted spectra after attribution)")
    parser.add_argument("-H", "--plot_residuals_histograms", dest="plot_residuals_histograms", action="store_true",
                        help="Make histograms of summed residuals (to use together with -R)")
    parser.add_argument("-S", "--plot_signatures", dest="plot_signatures", action="store_true",
                        help="Make signature plots (signature path and prefix arguments must be verified)")

    args = parser.parse_args()

    dataset_name = args.dataset_name
    mutation_type = args.mutation_type
    context = args.context
    signature_tables_path = args.signature_tables_path
    signatures_prefix = args.signatures_prefix
    input_folder = args.input_folder
    output_folder = args.output_folder + '/' + dataset_name + '/' + mutation_type + '/'
    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (DBS, ID).")
    elif mutation_type not in ['SBS', 'DBS', 'ID', 'SV', 'CNV']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID, SV, CNV" % mutation_type)

    if args.strand_bias and mutation_type!='SBS':
        raise ValueError("Can not plot strand bias information for %s mutation type." % mutation_type)
    elif args.strand_bias and mutation_type=='SBS' and context==96:
        raise ValueError("Strand information not available in 96 context." % mutation_type)

    if args.non_transcribed and not args.strand_bias:
        raise ValueError("Non-transcribed option must be used together with strand bias option (-b).")
    elif args.non_transcribed and not context in [288, 384]:
        raise ValueError("Non-transcribed region information is not available in %i context." % context)

    if args.plot_residuals_histograms and not args.plot_residuals:
        raise ValueError("Residuals histograms (-H) option must be used together with residuals option (-R).")

    if (args.plot_signatures and args.plot_residuals) or (args.plot_signatures and args.plot_fitted_spectra) or (args.plot_residuals and args.plot_fitted_spectra):
        raise ValueError("Please only use one of the following options: -R, -f, -S (residuals, fitted or signature plotting).")

    if args.plot_signatures:
        print("Making signature plots for %s mutation type..." % (mutation_type))
    elif args.plot_fitted_spectra:
        print("Making fitted spectra plots for dataset %s, %s mutation type..." % (dataset_name, mutation_type))
    elif args.plot_residuals:
        print("Making residuals spectra plots for dataset %s, %s mutation type..." % (dataset_name, mutation_type))
    else:
        print("Making mutation spectra plots for dataset %s, %s mutation type..." % (dataset_name, mutation_type))

    # index column treatment based on mutation type
    index_col = 0
    if mutation_type=='SBS':
        if context==96:
            index_col = [0,1]
        elif context in [192, 288]:
            index_col = [0,1,2]
        elif context==1536:
            index_col = 0
        else:
            raise ValueError("Context %i is not supported." % context)

    if args.plot_fitted_spectra:
        input_spectra = pd.read_csv('%s/%s/output_%s_%s_fitted_values.csv' % (input_folder, dataset_name, dataset_name, mutation_type), sep=None, index_col=index_col)
    elif args.plot_residuals:
        input_spectra = pd.read_csv('%s/%s/output_%s_%s_residuals.csv' % (input_folder, dataset_name, dataset_name, mutation_type), sep=None, index_col=index_col)
    elif args.plot_signatures:
        if mutation_type=='SBS' and context!=96:
            input_spectra = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type, context), sep=None, index_col=index_col)
        else:
            input_spectra = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=index_col)
    else:
        if mutation_type=='SBS':
            input_spectra = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_folder, dataset_name, dataset_name, context), sep=None, index_col=index_col)
        elif mutation_type=='DBS':
            input_spectra = pd.read_csv('%s/%s/WGS_%s.dinucs.csv' % (input_folder, dataset_name, dataset_name), sep=None, index_col=index_col)
        elif mutation_type=='ID':
            input_spectra = pd.read_csv('%s/%s/WGS_%s.indels.csv' % (input_folder, dataset_name, dataset_name), sep=None, index_col=index_col)
        else: # SV and CNV
            input_spectra = pd.read_csv('%s/%s/WGS_%s.%s.csv' % (input_folder, dataset_name, dataset_name, mutation_type), sep=None, index_col=index_col)

    if args.strand_bias:
        strand_bias_subfolder = 'TSB'
        strand_bias_suffix = '_TSB'
        if args.non_transcribed:
            strand_bias_subfolder += 'N'
            strand_bias_suffix += 'N'
    else:
        strand_bias_subfolder = ''
        strand_bias_suffix = ''

    if args.condensed:
        condensed_subfolder = 'condensed'
        condensed_suffix = '_condensed'
    else:
        condensed_subfolder = ''
        condensed_suffix = ''

    # limit the number of samples to analyse (if specified by --number option)
    if args.number_of_spectra!=-1:
        input_spectra = input_spectra.iloc[:,0:args.number_of_spectra]

    if args.samples:
        list_of_spectra = args.samples
    else:
        list_of_spectra = input_spectra.columns

    title_suffix = dataset_name.replace('_','/') + ', ' + mutation_type

    spectra_to_plot = input_spectra
    if args.plot_fitted_spectra:
        title_suffix += ', fitted'
        output_folder += 'fitted_spectra/'
    elif args.plot_residuals:
        title_suffix += ', residuals'
        output_folder += 'residuals/'
    elif args.plot_signatures:
        title_suffix = mutation_type
        output_folder += 'signature_plots/'

    # plotting residual histograms
    if args.plot_residuals and args.plot_residuals_histograms:
        # all residuals
        plot_residuals(spectra_to_plot.sum(), 'All residuals', savepath = output_folder + 'all_residuals.pdf')
        # context-dependent residuals
        if context==288 and mutation_type=='SBS':
            # 96-context residuals
            simplified_spectra = spectra_to_plot
            simplified_spectra = simplified_spectra.loc['T']+simplified_spectra.loc['U']+simplified_spectra.loc['N']
            simplified_spectra.index = simplified_spectra.index.map('_'.join)
            categories = simplified_spectra.index
            for category in categories:
                plot_residuals(simplified_spectra.loc[category], category, savepath = output_folder + '96_context/residuals_' + category + '.pdf')

            merged_filename = output_folder + '96_context_residuals.pdf'
            merge_plots('%s/96_context/*pdf' % (output_folder), merged_filename)

            # 6-context residuals
            simplified_spectra = spectra_to_plot.droplevel(2)
            simplified_spectra = simplified_spectra.groupby(level=1).sum()
            categories = mutation_categories[mutation_type]
            for category in categories:
                plot_residuals(simplified_spectra.loc[category], category, savepath = output_folder + '6_context/residuals_' + category + '.pdf')

            merged_filename = output_folder + '6_context_residuals.pdf'
            merge_plots('%s/6_context/*pdf' % (output_folder), merged_filename)

            # 6-context per strand residuals
            for strand in ['T', 'U', 'N']:
                simplified_spectra = spectra_to_plot.loc[strand]
                simplified_spectra = simplified_spectra.droplevel(1)
                simplified_spectra = simplified_spectra.groupby(level=0).sum()
                categories = mutation_categories[mutation_type]
                for category in categories:
                    plot_residuals(simplified_spectra.loc[category], strand +'-strand, ' + category,
                                   savepath = output_folder+ '6_context_per_strand/residuals_' + strand + '_' + category + '.pdf')

            merged_filename = output_folder + '6_context_per_strand_residuals.pdf'
            merge_plots('%s/6_context_per_strand/*pdf' % (output_folder), merged_filename)

    if args.plot_signatures:
        type = 'signatures'
        y_label = 'Relative contribution [%]'
    else:
        type = 'mutation_spectra'
        y_label = 'Mutation number'
        if args.relative:
            type += '_relative'
            y_label = 'Relative contribution [%]'

    for spectrum in list_of_spectra:
        burden = mutational_burden(spectra_to_plot[spectrum], absolute=False)
        if args.plot_signatures:
            title = '%s signature (%s)' % (spectrum, title_suffix)
        else:
            title = '%s sample: %i mutations (%s)' % (spectrum, burden, title_suffix)

        make_spectrum_plot(spectra_to_plot[spectrum], mutation_type=mutation_type, title=title, y_label=y_label, relative=args.relative,
                            show_errors = args.show_errors, strand_bias=args.strand_bias, add_non_transcribed=args.non_transcribed, condensed=args.condensed,
                            savepath=output_folder + '/%s/%s/%s/%s.pdf' % (type, condensed_subfolder, strand_bias_subfolder, spectrum))

    # plot composite spectra
    if args.relative:
        # make a normalised composite spectrum (normalisation by mutational burden for each spectrum, add up, then normalise by 100%)
        spectra_to_plot = spectra_to_plot.apply(lambda x: x/abs(sum(x)), axis=0)
    # otherwise, just add up spectra
    make_spectrum_plot(spectra_to_plot.sum(axis=1), mutation_type=mutation_type, title='Composite spectrum (%s)' % title_suffix, y_label=y_label, relative=args.relative,
                       show_errors = args.show_errors, strand_bias=args.strand_bias, add_non_transcribed=args.non_transcribed, condensed=args.condensed,
                       savepath=output_folder + '/composite_%s/composite_spectrum%s%s.pdf' % (type, condensed_suffix, strand_bias_suffix))

    # merge the sorted pdfs into a single one
    merged_filename = '%s/all_%s%s%s.pdf' % (output_folder, type, condensed_suffix, strand_bias_suffix)
    merge_plots('%s/%s/%s/%s/*pdf' % (output_folder, type, condensed_subfolder, strand_bias_subfolder), merged_filename)
