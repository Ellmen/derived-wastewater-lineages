# Learning novel viral lineages from wastewater sequencing data

This is the repository for running NMF on wastewater sequencing data to determine lineage definitions from mixed samples.
A manuscript describing the details of the method is coming soon.

## Instructions

This describes the workflow for a typical SARS-CoV-2 run (other viruses are similar)

1. Download or create a run containing Gromstole "coverage" and "mapped" csvs into the `data/sars-cov-2/runs` directory. A thin wrapper for running the Gromstole alignment on fastqs is provided in the `preprocess` directory but may require some changes.
2. Set the virus, number of lineages, fasta name, and runs (including the new run) variables in `find_lineages.py`. On the first run, all subsequent steps should be uncommented, however the mutation frequency matrix and learned nmf vectors are saved so subsequent runs can comment out the generation steps if the data is unchanged.
3. Run `python find_lineages.py`. This will create the mutations frequency matrix, learn the NMF vectors, and create a fasta: `data/sars-cov-2/[fasta_name].fasta` where `fasta_name` is specified in `find_lineages.py`.
4. (Optional) Analyze the data using a tool like pangolin or nextclade to determine which lineages have been discovered.
