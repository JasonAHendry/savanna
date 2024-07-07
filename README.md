# savanna
*Savanna* packages several bioinformatic pipelines for amplicon sequencing data generated with nanopore.

## Features
* Basecalling with [*Dorado*](https://github.com/nanoporetech/dorado)
* Sample demultiplexing with [*Dorado*](https://github.com/nanoporetech/dorado) or *Guppy*
* Read mapping with [*Minimap2*](https://github.com/lh3/minimap2) and amplicon coverage evaluation
* Variant calling with [bcftools](https://github.com/samtools/bcftools)

## Installation
### Docker
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

### From source
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

## Basic usage
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

## Testing
Example data is present in [example_data](https://github.com/JasonAHendry/savanna/tree/master/example_data).

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

