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
```
---

## Usage

**1. Identify Node-Specific SNPs**
Script: **identify_snps_per_node.py**
```bash
python3 identify_snps_per_node.py <alignment.fasta> <node_mapping.json> <outgroup> <output_snps.tsv>
```
- alignment.fasta: Multiple sequence alignment including outgroup reference.
- node_mapping.json: JSON dictionary mapping nodes to sample IDs.
- outgroup: Name of the outgroup.
- output_snps.tsv: Output file with SNPs per node (Position, Reference, Derived).

**2. Analyse Reads Covering SNPs**
Script: **analyse_snp_reads.py**
Given a BAM file and a node SNP definition file, this script analyzes read support, damage, and coverage.
```bash
python3 analyze_snp_reads.py <file.bam> <node_snps.tsv> <output_prefix>
```
- file.bam: Aligned sequencing reads (e.g., for ancient mtDNA).
- node_snps.tsv: Output from identify_snps_per_node.py.
- output_prefix: Prefix for output files.
Outputs:
- <prefix>_snp_read_details.csv: Per-read information (base, match, direction, deamination, etc.)
- <prefix>_snp_summary.csv: Summary per node (coverage, % derived, etc.)
- <prefix>_snp_read_analysis_log_*.log: Log file with processing details

## File Format Details
**node_mapping.json** example
```json
{
  "NodeA": ["Sample1", "Sample2"],
  "NodeB": ["Sample3", "Sample4"]
}
```
**node_snps.tsv** output format
|Node|Position|Reference|DerivedAllele|
|NodeA|152|T|C|
|NodeB|263|A|G|

## Example Workflow
```bash
# Step 1: Identify SNPs per node
python identify_snps_per_node.py alignment.fasta node_mapping.json snps.tsv

# Step 2: Analyze BAM for read-level SNP support
python analyse_snp_reads.py sample.bam snps.tsv results/sample
```

## Applications
- Ancient DNA analysis (e.g., identifying haplogroups from degraded mtDNA)
- Phylogenetic assignment
- Validation of SNP authenticity (damage vs. read variation vs. contamination)
