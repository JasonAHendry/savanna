# savanna
*Savanna* packages several bioinformatic pipelines for amplicon sequencing data generated with nanopore.

## Features
* Basecalling with [*Dorado*](https://github.com/nanoporetech/dorado)
* Sample demultiplexing with [*Dorado*](https://github.com/nanoporetech/dorado) or *Guppy*
* Read mapping with [*Minimap2*](https://github.com/lh3/minimap2) and amplicon coverage evaluation
* Variant calling

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

#### Steps
**1.  Clone the repository:**
```
git clone https://github.com/JasonAHendry/savanna.git
cd savanna
```
\
**2. Install `dorado`.**
Follow the installation instructions [here](https://github.com/nanoporetech/dorado).
- The `Dorado` wrapper will automatically download the appropriate basecalling models for you
\

**3. Install `guppy`.**
\
`dorado` does not currently support higher stringency demultiplexing, so we are still using `guppy`. In order to install `guppy`, you need to log onto the [nanopore community](https://community.nanoporetech.com) and naviigate to the Software  Downlaods section. At the bottom of the paged in Archived Software, you wll see Guppy v6.5.7 available for download.

\
**4.  Install other depedendencies with conda:**
```
conda env create -f environments/run.yml
```
or equivalently, with mamba:
```
mamba env create -f environments/run.yml
```
\
**5. Install `savanna` and remaining dependencies:**
```
pip install -e .
```
\
**6. Test your installation.**
In the terminal, you should see available commands by typing:
```
savanna --help
```

## Basic usage
First, download your desired reference genome. For example...
```
savanna download -r Pf3D7
```
downloads 3D7 reference genome of *P. falciparum*.

Basecall data using dorado with the command:
```
savanna basecall
```

Analyse demultiplexed FASTQ data with the command:
```
savanna analyse \
 -e <your_experiment_name> \
 -m <path/to/metadata.csv> \
 -b <path/to/your_regions.bed>
 --pipeline plasmo
```

\
Flag information:
- `-e`: Your experiment name. For example, '2023-06-12_seq-nomads8'. This is used to create the output directory.
- `-m`: Path to a metadata CSV file. 
  - This file *must* have a `barcode` and `sample_id` column.
  - Both `barcode` and `sample_id` columns must have only unique entries.
  - See `example_data/metadata/sample_info.csv` for an example.
- `-b`: Path to a BED file defining your regions of interest.
  - Note the BED file specification [here](https://en.wikipedia.org/wiki/BED_(file_format)).
  - This BED file *must* contain a forth column for region name information, and the names must be unique. They are used to generate plots and compute summary statistics.
  - See `example_data/beds/nomads8.amplicons.bed` for an example.

For a full running example look at [scripts/run_example.sh](https://github.com/JasonAHendry/savanna/blob/master/scripts/run_example.sh)

## Testing
Example data is present in [example_data](https://github.com/JasonAHendry/savanna/tree/master/example_data).

