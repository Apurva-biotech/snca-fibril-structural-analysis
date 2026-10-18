# Structural Analysis of Alpha-Synuclein Fibrils

## Overview
This project uses Python and Biopython to examine the structure and
inter-chain contacts of human alpha-synuclein (SNCA), a protein associated
with Parkinson's disease.

The analysis uses PDB entry 2N0A, an alpha-synuclein fibril structure
determined by solid-state NMR.

## Objectives
- Validate the alpha-synuclein sequence and parsed structure.
- Quantify inter-chain C-alpha contacts at defined distance cutoffs.
- Compare contact counts across neighboring-chain interfaces.
- Explore contact participation by sequence position and region.
- Produce reproducible CSV tables and figures.

## Methods
- Structure parsing: Biopython PDBParser.
- Contact atom: alpha carbon (C-alpha).
- Inter-chain contact: a C-alpha pair from different chains with a distance
  strictly below the selected cutoff.
- Cutoff sensitivity: 6, 8, and 10 angstroms.
- Regional analysis: N-terminal tail (1-30), fibril core (31-100),
  and C-terminal tail (101-140). These are exploratory analysis ranges.
- Interface analysis: all unique chain pairs; neighboring interfaces
  are A-B through I-J.

## Preliminary Results
At an 8 angstrom C-alpha cutoff:
- 2,633 inter-chain residue-pair contacts were recorded.
- The nine neighboring interfaces had a mean of 290.11 contacts.
- The sample standard deviation was 5.44 contacts.
- The coefficient of variation was 1.88%.

Direct calculation from the PDB gave:
- Below 6 angstroms: 1,219 contacts.
- Below 8 angstroms: 2,633 contacts.
- Below 10 angstroms: 5,430 contacts.

These are descriptive measurements of this structural model, not independent
biological observations.

## Project Structure
- `SNCA.fasta`: sequence input.
- `SNCA_verified.fasta`: sequence verification output.
- `SNCA_fibril_2N0A.pdb`: structure input.
- `scripts/`: analysis scripts.
- `results/`: generated CSV files and figures.

## Requirements
- Python 3
- Biopython
- Matplotlib

Activate the project's virtual environment before running scripts.

## Limitations
- C-alpha proximity is a geometric proxy and does not establish a hydrogen
  bond, side-chain interaction, binding affinity, or biological function.
- Contact counts depend on the distance cutoff and structural representation.
- Repeated chains within one model are not independent biological replicates.
- The selected residue-region boundaries are approximate and should be
  justified against structural literature before biological interpretation.
- These analyses alone cannot establish a mechanism of Parkinson's disease
  or identify validated drug-binding hotspots.

## Reproducibility
Run scripts from the project root. For example:

```bash
python scripts/analyze_fibril.py
python scripts/contact_cutoff_sensitivity.py
python scripts/interface_analysis.py
python scripts/region_analysis.py
python scripts/sequence_hotspots.py
python scripts/plot_sequence_hotspots.py
```

## Data Source
RCSB Protein Data Bank: https://www.rcsb.org/structure/2N0A

## Status
Exploratory structural bioinformatics project. Analyses and interpretations
are subject to further validation.
