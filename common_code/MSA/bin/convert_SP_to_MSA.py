""" convert_SP_to_MSA.py
Module to convert mutation tables from SigProfilerMatrixGenerator or signatures
tables from SigProfilerExtractor output paths specified by -i parameter, to
comma-separated, multi-indexed tables. Output files are saved in directory
specifiable by -o flag.
"""

import os
import glob
import copy
import pandas as pd
from argparse import ArgumentParser
from common_methods import make_folder_if_not_exists

def compare_index(first, second):
    if not first.index.equals(second.index):
        print('Converted index:', first.index.to_list())
        print('Target index:', second.index.to_list())
        raise ValueError("Index mismatch, check your input data.")
    return

def convert_index(input_dataframe, context=96):
    input_table = copy.deepcopy(input_dataframe)
    if context==96:
        input_table = input_table.sort_index(level=0)
        for element in input_table.index:
            sub = element.split('[', 1)[1].split(']')[0]
            replaced_element = element.replace('['+sub+']', sub[0])
            replaced_element = sub + ':' + replaced_element
            input_table.rename(index={element:replaced_element}, inplace=True)
        input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
    elif context==192 or context==384:
        if context==192:
            input_table = input_table[~input_table.index.str.contains("B:")]
            input_table = input_table[~input_table.index.str.contains("N:")]
        for element in input_table.index:
            sub = element.split('[', 1)[1].split(']')[0]
            replaced_element = element.replace('['+sub+']', sub[0])
            replaced_element = replaced_element.replace(':',':' + sub + ':')
            input_table.rename(index={element:replaced_element}, inplace=True)
        input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
    elif context==288:
        for element in input_table.index:
            sub = element.split('[', 1)[1].split(']')[0]
            replaced_element = element.replace('['+sub+']', sub[0])
            replaced_element = replaced_element.replace(':',':' + sub + ':')
            input_table.rename(index={element:replaced_element}, inplace=True)
        input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
    return input_table

if __name__ == '__main__':
    parser = ArgumentParser(description='Convert SigProfilerMatrixGenerator output to MSA format')
    parser.add_argument("-i", "--input_folder", dest="input_path",
                      help="set path to SigProfiler output (matrix generator or extractor)")
    parser.add_argument("-S", "--reindex_signatures", dest="reindex_signatures", action="store_true",
                      help="reindex signature tables (assume SP extractor output)")
    parser.add_argument("-C", "--COSMIC", dest="COSMIC", action="store_true",
                      help="assume COSMIC signatures (forced 96 context in SBS)")
    parser.add_argument("-d", "--dataset_name", dest="dataset_name", default='',
                        help="set dataset name to use in converted filenames")
    parser.add_argument("-t", "--mutation_types", nargs='+', dest="mutation_types", default=['SBS','DBS','ID'],
                      help="set mutation types, e.g. -t SBS DBS ID (default)")
    parser.add_argument("-c", "--contexts", nargs='+', dest="contexts", type=int, default=[96, 288, 1536],
                      help="set SBS contexts e.g. -c 96 288 1536 (default). Supported contexts: 96, 192, 288, 1536")
    parser.add_argument("-o", "--output_path", dest="output_path", default='./',
                        help="set output path for converted tables (default: ./)")
    # parser.add_argument("-e", "--exome", dest="exome", action="store_true",
    #                   help="Treat input matrices of the exome regions of the genome, otherwise assumme WGS")
    parser.add_argument("-s", "--signature_path", dest="signature_tables_path", default='signature_tables',
                      help="set path to signature tables to extract indexes (default: signature_tables)")
    parser.add_argument("-p", "--input_signatures_prefix", dest="input_signatures_prefix", default='sigProfiler',
                      help="set prefix in signature filenames to extract indexes (sigProfiler by default)")
    parser.add_argument("-n", "--reindexed_signatures_prefix", dest="reindexed_signatures_prefix", default='sigProfilerNew',
                      help="set prefix for reindexed signatures filenames (sigProfilerNew by default)")
    parser.add_argument("-E", "--use_extractor_for_mutation_tables", dest="use_extractor_for_mutation_tables", action="store_true",
                      help="use Samples.txt file from SigProfiler extractor output, rather than matrix generator output")

    options = parser.parse_args()
    input_path = options.input_path
    dataset_name = options.dataset_name
    mutation_types = options.mutation_types
    contexts = options.contexts
    output_path = options.output_path
    signature_tables_path = options.signature_tables_path
    input_signatures_prefix = options.input_signatures_prefix
    reindexed_signatures_prefix = options.reindexed_signatures_prefix

    if not input_path:
        raise ValueError("Please specify the input path to SigProfiler tables using -i flag.")

    if not options.reindex_signatures:
        if not dataset_name:
            raise ValueError("Please specify an arbitrary dataset name for use in MSA pipeline execution (see -h for help)")
        print('Converting mutation tables for dataset', dataset_name, ' mutation types', mutation_types, ', considering SBS contexts', contexts)
        print('SigProfilerMatrixGenerator output to parse:', input_path)
    else:
        if options.COSMIC:
            signatures_type = 'COSMIC'
        else:
            signatures_type = 'De-Novo'
        print('Converting %s signature tables from SigProfilerExtractor output %s' % (signatures_type, input_path))

    for mutation_type in mutation_types:
        if mutation_type not in ['SBS', 'DBS', 'ID', 'SV', 'CNV']:
            raise ValueError("Unsupported mutation type: %s. Supported types: SBS, DBS, ID, SV, CNV" % mutation_type)

    signatures = {}
    for mutation_type in mutation_types:
        if mutation_type=='SBS':
            for context in contexts:
                if context == 96:
                    index_col = [0,1]
                    signatures[mutation_type + str(context)] = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, input_signatures_prefix, mutation_type), sep=',', index_col=index_col)
                elif context in [192, 288]:
                    index_col = [0,1,2]
                    signatures[mutation_type + str(context)] = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, input_signatures_prefix, mutation_type, context), sep=',', index_col=index_col)
                elif context == 1536:
                    index_col = 0
                    signatures[mutation_type + str(context)] = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, input_signatures_prefix, mutation_type, context), sep=',', index_col=index_col)
        else:
            signatures[mutation_type] = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, input_signatures_prefix, mutation_type), sep=',', index_col=0)

    if options.reindex_signatures:
        # reindex signatures
        for mutation_type in mutation_types:
            input_signatures = glob.glob(input_path + '/%s*/Suggested_Solution/*%s*/Signatures/*Signatures.txt' % (mutation_type, signatures_type))
            if not input_signatures:
                raise ValueError("Can't find any signature tables of type %s, mutation type %s in input path %s" % (mutation_type, signatures_type, input_path) )
            for signature_table in input_signatures:
                if mutation_type=='SBS':
                    for context in contexts:
                        if not str(context) in signature_table:
                            continue
                        signature_table_to_reindex = pd.read_csv(signature_table, sep='\t', index_col=0)
                        print('Converting signature table %s (%s mutation type, %i context)' % (signature_table, mutation_type, context))
                        template_context = context if not options.COSMIC else 96
                        reindexed_signatures = convert_index(signature_table_to_reindex, context=template_context)
                        template_signatures = signatures[mutation_type + str(template_context)]
                        reindexed_signatures = reindexed_signatures.reindex(template_signatures.index)
                        compare_index(reindexed_signatures, template_signatures)
                        reindexed_signatures_filename = '%s/%s_%s_%i_signatures.csv' % (output_path, reindexed_signatures_prefix, mutation_type, context)
                        if context==96:
                            reindexed_signatures_filename = reindexed_signatures_filename.replace('_%i' % context,'')
                        reindexed_signatures.to_csv(reindexed_signatures_filename, sep = ',')
                        print('Done. Check the output signature table:', reindexed_signatures_filename)
                else:
                    if mutation_type=='DBS':
                        context = 78
                    elif mutation_type=='ID':
                        context = 83
                    elif mutation_type=='SV':
                        context = 32
                    elif mutation_type=='CNV':
                        context = 48
                    if not str(context) in signature_table:
                        continue
                    signature_table_to_reindex = pd.read_csv(signature_table, sep='\t', index_col=0)
                    print('Converting signature table %s (%s mutation type, %i context)' % (signature_table, mutation_type, context))
                    # simply overwrite index for non-SBS mutation types
                    reindexed_signatures = signature_table_to_reindex
                    template_signatures = signatures[mutation_type]
                    reindexed_signatures.index = template_signatures.index
                    compare_index(reindexed_signatures, template_signatures)
                    reindexed_signatures_filename = '%s/%s_%s_signatures.csv' % (output_path, reindexed_signatures_prefix, mutation_type)
                    reindexed_signatures.to_csv(reindexed_signatures_filename, sep = ',')
                    print('Done. Check the output signature table:', reindexed_signatures_filename)
    else:
        # reindex mutation tables
        make_folder_if_not_exists(output_path + '/' + dataset_name)
        for mutation_type in mutation_types:
            if options.use_extractor_for_mutation_tables:
                input_files = glob.glob(input_path + '/%s*/Samples.txt' % mutation_type)
            else: # assume matrix generator output
                input_files = glob.glob(input_path + '/%s/*%s*' % (mutation_type, mutation_type))
            if not input_files:
                raise ValueError("Can't find any files of mutation type %s in input path %s" % (mutation_type, input_path) )
            for file in input_files:
                if mutation_type=='SBS':
                    for context in contexts:
                        if context==192 and not '384' in file:
                            continue
                        if context!=192 and not str(context) in file:
                            continue
                        input_table = pd.read_csv(file, sep='\t', index_col=0)
                        print('Converting:', mutation_type, context, file)
                        input_table = convert_index(input_table, context=context)
                        signature_table = signatures[mutation_type + str(context)]
                        if signature_table is not None:
                            input_table = input_table.reindex(signature_table.index)
                            compare_index(input_table, signature_table)
                        new_filename = output_path + '/%s/WGS_%s.%i.csv' % (dataset_name, dataset_name, context)
                        if context==192:
                            new_filename = new_filename.replace('384','192')
                        input_table.to_csv(new_filename, sep = ',')
                else:
                    if mutation_type=='DBS' and not '78' in file:
                        continue
                    if mutation_type=='ID' and not '83' in file:
                        continue
                    if mutation_type=='SV' and not '32' in file:
                        continue
                    if mutation_type=='CNV' and not '48' in file:
                        continue
                    # simply overwrite index for other mutation types (equality assumption)
                    print('Converting:', mutation_type, file)
                    input_table = pd.read_csv(file, sep='\t', index_col=0)
                    input_table.index = signatures[mutation_type].index
                    if mutation_type == 'DBS':
                        new_filename = output_path + '/%s/WGS_%s.dinucs.csv' % (dataset_name, dataset_name)
                    elif mutation_type == 'ID':
                        new_filename = output_path + '/%s/WGS_%s.indels.csv' % (dataset_name, dataset_name)
                    else: # SV and CNV
                        new_filename = output_path + '/%s/WGS_%s.%s.csv' % (dataset_name, dataset_name, mutation_type)
                    input_table.to_csv(new_filename, sep = ',')
        print('Done, please check the outputs:', output_path + '/' + dataset_name)
