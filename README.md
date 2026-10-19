# SynDrugMD

SynDrugMD is an annotation-based framework relating drug combination synergy to synthetic-lethality (SL) and pathway information.

**Index Terms:** Drug combination synergy, synthetic-lethality annotation, pathway relationships

## Overview

- The reported collection contains **555,002 records** across **3,861 drugs** and **229 cell lines**.
- Raw measurements were reprocessed using a common protocol; collected multi-omics resources were not computationally integrated.
- The reported selected SL subset contains **14,041 annotation-positive records** among **22,881 synergistic records (61.37%)**.
- Pathway analyses use **216,030 target-annotated records**.
- **OPC (Overlap Pathway Count)** measures distinct common pathways, while **RPC (Regulatory Pathway Count)** counts direct pathway relations.
- The analyses report group counts, mean differences, cluster-aware intervals and multiplicity-adjusted comparisons, retaining 17 original tissue/lineage labels.
- Raw OPC/RPC differences attenuated after source, target/pathway-size and degree standardization; intervals for the fully standardized normalized metrics included zero.
- Threshold profiles report both conditional class proportions and retained counts.

These results describe size-sensitive annotation associations and do not establish mechanisms, independently validated prediction, or clinical efficacy.

## Repository Structure

- `data/` — integrated datasets: drug combinations, cell lines, drug targets, gene expression / copy number, pathway and protein annotations, and the human SL reference (`Human_SL.csv`).
- `analysis/` — analysis code and result tables, including the synthetic-lethality (SL) subset analyses and the OPC/RPC pathway-metric comparisons (`OPC_RPC/`).
- `figures/` — figures and plotting notebooks for the reported statistics.
- `SynDrugMD_statistics_revised.pdf` — the paper.
- `requirements.txt` — Python dependencies.

Large CSV files are tracked with [Git LFS](https://git-lfs.com/); clone with `git lfs pull` after cloning.

## Reference

SynDrugMD: An Annotation-Based Framework for Synthetic Lethality and Pathway Analysis of Drug Combination Synergy
