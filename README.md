# savanna
## Overview
`savanna` supports in-depth analysis of nanopore amplicon sequencing data.

## Install

#### Requirements
To install `savanna`, you will need:
- The version control software [git](https://github.com/git-guides/install-git)
- The package manager [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation.html) 
  - Mamba is faster and is recommended

#### Steps
**1.  Clone the repository:**
```
git clone https://github.com/JasonAHendry/nomadic3.git
cd nomadic3
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

**A. Download your reference genome** 

`savanna` performs real-time mapping to a reference genome. Start by downloading the reference genome for your target organism, e.g.:
```
savanna gather -r Pf3D7
```
For the 3D7 reference genome of *P. falciparum*.

\
**B. Run `nomadic realtime`**

Once sequencing has started, you can perform real-time analysis using the `nomadic realtime` command as follows:

```
savanna run \
 -e <your_experiment_name> \
 -m <path/to/metadata.csv> \
 -b <path/to/your_regions.bed>
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

For a full running example look at `scripts/run_realtime.sh`

## Testing
I have created a script to simulate a small nanopore sequencing run, that allows you to test `nomadic3 realtime` without having an actual sequencing experiment running. To try this, first run: 

```
./scripts/run_realtime.sh
```

Then, in a second terminal window, and run:
```
python scripts/simulate_sequencing.py
```




