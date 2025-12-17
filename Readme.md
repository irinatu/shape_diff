# SHAPE Reactivity Comparison Tool

This script compares two RNA SHAPE(-MaP) reactivity profiles provided in `.map` format
(e.g. *in vitro* vs *in cellulo*) and identifies regions with significantly different
reactivities. The method is described by Smola et al. (In-cell RNA structure probing with SHAPE-MaP. Nature Protocols 2018 13:6, 13, 1181–1195.)

## Dependencies

### Python packages

The following Python libraries are required:

- **Python ≥ 3.10**
- **numpy**
- **pandas**
- **matplotlib**
- **scipy**
- **varnaapi** *(optional, required only if `-f` option is used)*

They can be installed using:

```bash
pip install numpy pandas matplotlib scipy varnaapi
```
### External software

VARNA

For RNA secondary structure visualization, the script uses VARNA via the Python interface varnaapi.
VARNA must be installed separately as a .jar file
The path to the VARNA .jar file must be provided via an environment variable called VARNA_JAR. This needs to be done once, and then the script can be used without specifying the path again.
macOS 
```bash
echo 'export VARNA_JAR=/path/to/VARNAv3-94a.jar' >> ~/.zshrc
source ~/.zshrc
```

HPC / bashrc
```bash
echo 'export VARNA_JAR=/path/to/VARNAv3-94a.jar' >> ~/.bashrc
source ~/.bashrc
```
You can verify that the variable is set correctly by running:
```bash
echo $VARNA_JAR
```

The -f / --figures option will fail gracefully if VARNA_JAR is not set.
VARNA figure generation (-f option) works best for RNAs up to ~200 nt
More information about VARNA:
https://varna.lri.fr

## Notes on environment
The script was developed and tested on macOS with the following Python and library versions:
Python  3.10.4
numpy 2.2.6
pandas 2.2.3
matplotlib 3.9.2
scipy 1.15.3 
varnaapi 1.0.0

Other operating systems are supported, provided that Java and VARNA are correctly installed
Missing SHAPE reactivities must be encoded as -999

## Features

- Smoothing of SHAPE reactivities and standard errors (rolling window = 3 nt)
- Calculation of:
  - ΔSHAPE (difference in smoothed reactivities)
  - Z-factor–like metric
  - Standardized score (S)
- Identification of significantly different regions using a sliding window criterion
- Spearman correlation between two SHAPE profiles
- Generation of:
  - tab-delimited output files
  - `.shape` files for VARNA visualization
  - comparison plots
  - optional VARNA structure figures

## Input files

- SHAPE reactivity files in `.map` format (4 columns):
  1. nucleotide index
  2. reactivity value (`-999` for missing data)
  3. standard error
  4. nucleotide
- RNA secondary structure in dot-bracket format (for VARNA visualization)

## Usage

```bash
python shape_diff_weeks.py \
  -v examples/input/cond1.map \
  -c examples/input/cond2.map \
  -d test.dot \
  -t first \
  -e second \
  -o -s -f

