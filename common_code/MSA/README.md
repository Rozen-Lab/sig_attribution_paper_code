# Mutational Signature Attribution pipeline

![logo](MSA.png)

Mutational signature attribution analysis, featuring automised optimisation study with simulated data.

## Introduction
The purpose of this Readme is to provide a guide for the quick start of using MSA. An extensive Wiki page detailing the usage of this tool can be found [here](https://gitlab.com/s.senkin/MSA/-/wikis/home).

## Running with Nextflow
The best way to run the code is by using [Nextflow](https://www.nextflow.io/).
Once you have installed Nextflow, run the test job locally or on your favourite cluster:

```
nextflow run https://gitlab.com/s.senkin/MSA -profile docker
```

If you don't have [docker](https://www.docker.com/) installed, you can also use [conda](https://conda.io) or [singularity](https://sylabs.io/singularity/) profiles (docker and singularity are recommended).

The pipeline should run and produce all the results automatically. You can also retrieve the code ([see below](https://gitlab.com/s.senkin/MSA#getting-started)) in order to adjust all the inputs and parameters. In the [run_auto_optimised_analysis.nf](run_auto_optimised_analysis.nf) file various parameters can be specified.

## Running on SigProfiler output

MSA natively supports [SigProfilerExtractor](https://github.com/AlexandrovLab/SigProfilerExtractor) and [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator) outputs.

The simplest way to run is as follows:

```
nextflow run https://gitlab.com/s.senkin/MSA -profile docker --dataset SP_test \
                              --SP_extractor_output_path /full/path/to/SP_extractor_output/
```

If [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator) output is provided, it will take priority over the [SigProfilerExtractor](https://github.com/AlexandrovLab/SigProfilerExtractor) one for input mutation matrices:

```
nextflow run https://gitlab.com/s.senkin/MSA -profile docker --dataset SP_test \
                              --SP_matrix_generator_output_path /full/path/to/SP_matrix_generator_output/ \
                              --SP_extractor_output_path /full/path/to/SP_extractor_output/
```

## Options

All parameters are described in the dedicated [wiki page](https://gitlab.com/s.senkin/MSA/-/wikis/Parameters-description-table). Most general parameters are listed below.

### General parameters

| Parameters  | Default value | Description |
|-----------|-------------|-------------|
| --help | null | print usage and optional parameters |
| --SP_matrix_generator_output_path | null | optionally use SigProfilerMatrixGenerator output from specified **full path** |
| --SP_extractor_output_path | null | optionally use SigProfilerExtractor output from specified **full path** to attribute signatures extracted by SigProfiler |
| --dataset | SIM_test | set the name of the dataset. If no SigProfiler output is provided, the matrices must exist in params.input_tables folder (see example) |
| --input_tables | $baseDir/input_mutation_tables | if not using SigProfiler outputs, **full path** to input mutation tables, the repository one is used by default |
| --signature_tables | $baseDir/signature_tables | if not using SigProfiler outputs, **full path** to input signature tables, the repository one is used by default |
| --signature_prefix | sigProfiler | if not using SigProfiler outputs, prefix of signature files to use, must be located in signature_tables folder (e.g. sigProfiler, sigRandom) |
| --output_path | . | output path for plots and tables |
| --mutation_types | \['SBS', 'DBS', 'ID'\] | mutation types to analyse. Only one can be specified from command line, or a list in the run_auto_optimised_analysis.nf file |
| --number_of_samples | -1 | number of samples to analyse (-1 means all available) |
| --SBS_context | 96 | SBS context to use (96, 192, 288 or 1536) |
| --COSMIC_signatures | false | if set to true, COSMIC signatures are used form SigProfiler output, otherwise de-novo ones are used |

## Running manually

### Getting started

Retrieve the code:
```
git clone https://gitlab.com/s.senkin/MSA.git
cd MSA
```

Run the fully automised pipeline with optimisation on the test dataset (or adjust parameters/input matrices accordingly, specifying the **full path** to them):
```
nextflow run run_auto_optimised_analysis.nf -profile docker --dataset SIM_test --input_tables $PWD/input_mutation_tables --output_path test
```

Alternatively, to run the pipeline without optimisation, using fixed penalties on the test dataset:
```
nextflow run run_analysis.nf -profile docker --dataset SIM_test --weak_threshold 0.02 --input_tables $PWD/input_mutation_tables --output_path test
```

### Setting up dependencies

If you can not run *nextflow*, you can still run some basic analysis manually (scripts in the *./bin* folder).
Dependencies so far are: *pandas*, *numpy*, *scipy*, *matplotlib* and *seaborn*. If you don't have them, the easiest way is to set up the virtual environment using [conda](https://conda.io) package manager:

```
conda env create -f environment.yml
```

This only needs to be done once. Afterwards, just activate the environment whenever needed:

```
source activate msa
```

Alternatively, you can use *docker* yourself with the *Dockerfile* provided, or use ready-made images ([docker](https://hub.docker.com/r/ssenkin/msa/tags) or [singularity](https://cloud.sylabs.io/library/ssenkin/default/msa)).


### Simulating data (if not using automised pipeline)

* [input_mutation_tables/SIM](input_mutation_tables/SIM) folder contains a set produced with existing (PCAWG or COSMIC) signatures. [np.random.normal](https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.normal.html) function was used to generate normal distributions of mutational burdens corresponding to each PCAWG signature mentioned in the [signatures_to_generate](bin/simulate_data.py#L9) dictionary in the script, containing Gaussian means and standard deviations for each signature.
* Note that the distributions are not strictly Gaussian since negative numbers of burdens are replaced by zeros
* To reproduce the simulated set of samples with reshuffled *SBS1/5/22/40* PCAWG signatures, one can run the following script (without the *-r* option):
```
python bin/simulate_data.py -t SBS -c 96 -n 100 -s signature_tables
```

* [input_mutation_tables/SIMrand](input_mutation_tables/SIMrand) folder contains a set of 100 simulated samples for 96/192 contexts SBS, as well as dinucs and indels, where each sample contains contributions from **5** randomly selected signatures out of **100** Poisson-generated signatures. To reproduce (e.g. for 96-context SBS, 100 signatures and samples), run:
```
python bin/generate_random_signatures.py -t SBS -c 96 -n 100
python bin/simulate_data.py -r -t SBS -c 96 -n 100 -d SIMrand
```
In both scripts, a normal distribution can be used to generate white noise using *-z* option, with a Gaussian centred around **0** for each category of mutations, with standard deviation set by *-Z* option (**2** by default). Additional flags can be viewed in the code or using *-h* option in each script.

The file format produced is the same as that of the PCAWG dataset.
