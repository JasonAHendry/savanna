<p align="center"><img src="misc/savanna_logo-v2.png" width="500"></p>

## Overview
*Savanna* is a tool for analysing targeted nanopore sequencing data. It supports a workflow of basecalling, demultiplexing, and downstream analysis $-$ locally or using High-performance computing (HPC). A variety of analysis pipelines are available.

Please note *Savanna* is still in early stage development and your feedback is welcome.

## Features
- [x] Basecalling with [*Dorado*](https://github.com/nanoporetech/dorado)
- [x] Sample demultiplexing with [*Dorado*](https://github.com/nanoporetech/dorado) or *Guppy*
- [x] Read mapping with [*Minimap2*](https://github.com/lh3/minimap2)
- [x] Sample quality control and amplicon coverage evaluation
- [x] Variant calling with [*bcftools*](https://github.com/samtools/bcftools)
- [x] Support for multiple species 

## Installation
### Docker
<details>
  
#### Requires
* [*Docker*](https://www.docker.com/)
  
#### Steps
```
docker pull jasonahendry/savanna:0.0
```
This will download an image that already has `dorado`, `savanna`, and all dependencies pre-installed. Unfortunately it is a bit more cumbersome to run from the command line:

```
docker run -w `pwd` -v `pwd`:`pwd` jasonahendry/dorado:0.0 savanna
```
</details>

### From source
<details>
  
#### Requires
- The version control software [*Git*](https://github.com/git-guides/install-git)
- The package manager [*Conda*](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [*Mamba*](https://mamba.readthedocs.io/en/latest/installation.html) 
  - Mamba is faster and is recommended
- [*Dorado*](https://github.com/nanoporetech/dorado) in must be installed and inside of `$PATH` for `savanna basecall`
- [*Dorado*](https://github.com/nanoporetech/dorado) or [*Guppy*](https://community.nanoporetech.com) must be installed for `savanna demultiplex`

#### Steps
**1.  Clone the repository:**
```
git clone https://github.com/JasonAHendry/savanna.git
cd savanna
```

**2.  Install other depedendencies with conda:**
```
conda env create -f environments/run.yml
```
or equivalently, with mamba:
```
mamba env create -f environments/run.yml
```
\
**3. Install `savanna` and remaining dependencies:**
```
pip install -e .
```
\
**4. Test your installation.**
In the terminal, you should see available commands by typing:
```
savanna --help
```
</details>

## Basic usage
*Savanna* has four main subcommands which can be viewed by typing `savanna --help`:
```
Usage: savanna [OPTIONS] COMMAND [ARGS]...

  Analyse targeted nanopore sequencing data for genomic surveillance

Options:
  --help  Show this message and exit.

Commands:
  download     Download reference genomes.
  basecall     POD5 to FASTQ.
  demultiplex  FASTQ to per-sample FASTQ.
  analyse      Per-sample FASTQ to results.
```

**A. Download your reference genome** 

In most cases you will want to start by downloading your reference genome(s) of interest. For example, to download the *P. falciparum* reference genome, you would run:
```
savanna download -r Pf3D7
```
Both the FASTA files and GFF files for the reference will be downloaded.
> **Note:** You will need a stable internet connection for this step.

**B. Analyse results** 

If you have already basecalled and demultiplexed your data (e.g. using MinKNOW), then the next step will be to analyse the data using `savanna analyse`. As an example, the following command will analyse the example data for NOMADS8 provided in the github repository:
```
savanna analyse \
-e 0000-00-00_expt1 \
-f example_data/expt1/fastq_pass \
-m example_data/expt1/metadata/sample_info.csv \
-r example_data/expt1/metadata/nomads8.amplicons.bed \
--pipeline plasmo 
```
Here is a breakdown of key *flags*:
| Flag | Description | Required / Optional |
| ---    | --- | --- |
| ` -e ` | Name of the experiment, used as output directory name. E.g. '2023-05-12_exptA'. | Required |
| ` -f `   | Path to directory containing demultiplexed FASTQ files (e.g. '<path>/<to>/fastq_pass'). Typically produced by MinKNOW, dorado, guppy or with savanna demultiplex. | Required |
| ` -m `   | Path to metadata CSV file containing barcode and sample information. Required to contain `barcode` and `sample_id` columns; can optionally contain other columns of relevance. See [here](https://github.com/JasonAHendry/savanna/blob/master/example_data/expt1/metadata/sample_info.csv) for an example. | Required |
| ` -r `   | Path to BED file specifying genomic regions of interest. See [here](https://github.com/JasonAHendry/savanna/blob/master/example_data/expt1/metadata/nomads8.amplicons.bed) for an example.  | Required |
| ` -p `  | Name of the pipeline to be run. Default is `plasmo` for *P. falciparum*. | Optional |
| ` -b `  | Analyse only a single barcode from the experiment, indicated by an integer. E.g. to analyse `barcode03` you would include `-b 3` | Optional |
| ` -s `  | Only run experiment-wide summary. Mainly useful for running *Savanna* in HPC environments. | Optional |


## Testing
Example data is present in [example_data](https://github.com/JasonAHendry/savanna/tree/master/example_data) and example scripts are present in [scripts](https://github.com/JasonAHendry/savanna/tree/master/scripts).

## Development
#### Creating a new analysis module
* Create a new directory inside of src/savanna/analyse
* Implement a `BarcodeAnalysis` subclass
  * This specifies what the analysis will do for a single barcode
  * Pass parameters of interest to the initialisation method
* Implement a `ExperimentAnalysis` subclass
  * This will automatically run a `BarcodeAnalysis` across an entire experiment
  * Optionally allows for outputs to be summarised across barcodes, and plots created
* Use your `ExperimentAnalysis` in an existing or new `Pipeline` subclass
* Invoke the pipeline using the `--pipeline` flag of `savanna analyse`

## Acknowledgements
This work was funded by the Bill and Melinda Gates Foundation (INV-003660, INV-048316).

