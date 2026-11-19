# Linking gene expression at peak bleaching to survival in *Acropora cervicornis*

Which genes expressed at peak bleaching (October 2024) in *Acropora cervicornis* are associated with later survival (March 2025), independent of depth?

**Status:** DRAFT for [course name]. STAR alignment and DESeq2 are complete; annotation of candidate genes and figure polishing are in progress.

## Repository contents

| Path | What it is |
|---|---|
| `README.md` | This file |
| `ai_usage.md` | Log of AI use (prompts, what was used, how it was checked) |
| `scripts/` | All analysis code (see "Run order") |
| `metadata/sample_table.csv` | The 32 samples analyzed: genotype, depth, fate, notes |
| `results/` | Count matrix, mapping QC table, DESeq2 results tables |
| `figures/` | PCA, volcano plot, heat map |
| `env/` | `environment.yml` (conda) and `sessionInfo.txt` (R packages) |

## Data

- RNA-seq: 48 paired-end (2 x 151 bp) October 2024 samples, trimmed with fastp (Mydlarz lab; unpublished).
- Reference: *A. cervicornis* genome and annotation (Selwyn & Vollmer, 2023) [add the accession you used].
- Raw reads, BAM files and the STAR index are **not** in this repo (size and data ownership).

## Run order

1. `scripts/config.sh`: set paths and activate the conda environment
2. `scripts/00_make_samplesheet.sh`: build `samples.tsv` from the trimmed reads (BL_* samples only)
3. `scripts/01_star_index.sbatch`: build the STAR index from the FASTA and GTF
4. `scripts/02_star_align.sbatch`: align each sample, with `--quantMode GeneCounts` (SLURM array)
5. `scripts/03_merge_counts.py`: merge per-sample counts into `acer_counts.tsv` and write `acer_mapping_qc.tsv`
6. `scripts/04_deseq2.R`: DESeq2 (Models A and B), PCA, volcano plot, heat map

[Also list your STAR R Markdown file here and say what it covers.]

```bash
# example: on the HPC, from the scripts folder
bash 00_make_samplesheet.sh
sbatch 01_star_index.sbatch
sbatch --array=1-48%10 02_star_align.sbatch
python3 03_merge_counts.py <star_out_dir> <output_prefix>
# then, locally in R, with acer_counts.tsv and sample_table.csv in the working folder:
Rscript 04_deseq2.R
```

## Analysis summary

- **Samples:** 48 sequenced; 32 analyzed. Excluded: 5 with under 1M gene-assigned reads, source-colony samples with no depth, samples with ambiguous genotype (294/300?), and samples absent from the metadata.
- **Fate:** a sample was labeled "Alive" if a March 2025 sample from the same genotype and depth existed, otherwise "Dead" (a group-level proxy, not per-fragment survival).
- **Models (DESeq2):** A: `~ depth + fate`; B: `~ genotype + depth + fate`. Genes at BH-adjusted p < 0.05 were called significant.
- **Key result so far:** [one sentence on P5673_016782].

## Software versions

STAR [2.7.11b], Python [version], pandas [version], R [version], DESeq2 [version], apeglm [version], ggplot2 [version], pheatmap [version]. Full lists: `env/environment.yml`, `env/sessionInfo.txt`.

## Data and ethics notes

Sequencing data and metadata belong to the Mydlarz lab and are unpublished. [State whether count/metadata files are included and that the lab approved sharing.]

## AI use

See `ai_usage.md`.

## Author

Mason Gorrondona, [course/institution]. [Contact.]
