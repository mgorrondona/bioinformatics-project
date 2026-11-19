# bioinformatics-project
Linking gene expression at peak bleaching to survival in Acropora cervicornis

Which genes expressed at peak bleaching (October 2024) in Acropora cervicornis are associated with later survival (March 2025), independent of depth?

Status: DRAFT for [course name]. STAR alignment and DESeq2 are complete; annotation of candidate genes and figure polishing are in progress.

Repository contents
Path	What it is
README.md	This file
ai_usage.md	Log of AI use (prompts, what was used, how it was checked)
scripts/	All analysis code (see "Run order")
metadata/sample_table.csv	The 32 samples analyzed: genotype, depth, fate, notes
results/	Count matrix, mapping QC table, DESeq2 results tables
figures/	PCA, volcano plot, heat map
env/	environment.yml (conda) and sessionInfo.txt (R packages)
