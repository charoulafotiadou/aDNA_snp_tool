# aDNA_snp_tool
A set of tools for identifying and analyzing phylogenetically informative SNPs across nodes using multiple sequence alignments and read-level data. Designed for ancient DNA and mitochondrial datasets, it enables detection of node-specific SNPs and read-based evaluation of support, deamination, and coverage patterns.

## Features

- Identify **node-specific SNPs** from a multiple sequence alignment and node mapping.
- Analyze BAM files to assess:
  - Read-level support for SNPs
  - Deamination damage patterns (C→T / G→A)
  - Coverage and directionality (forward/reverse reads)
- Export SNP summaries and detailed read-level reports.

---

## Installation

### Requirements

- Python 3.7+
- `pandas`
- `pysam`
- `biopython`

### Install dependencies

```bash
pip install pandas pysam biopython
