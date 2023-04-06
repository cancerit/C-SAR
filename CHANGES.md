# Changes

## 1.3.7

- Update RCRISPR to 1.2.2.0 (allows numeric sample names in count and fold change matrices)

## 1.3.6

- Update BAGEL config scaled_bf_infile_gene_column_index to 2

## 1.3.5

- Requires RCRISPR >= 1.2.1.0
- Allow QC PCA plots to use shape when groups not enabled
- Use BAGEL sgRNA BF file to calculate gene-level BFs (average of sgRNA BFs per gene)

## 1.3.4

Add contrast (e.g. treatment_vs_plasmid) and input type (e.g. counts or fc) to filename suffix for intermediate QC BAGEL statistics files.

- With NNMD (gene level fold changes): `bagel_classification_summary.treatment_vs_plasmid.fc.tsv`
- Without NNMD (sgrna counts): `bagel_classification_summary.treatment_vs_plasmid.counts.tsv`

## 1.3.3

RRA processing enriched/depleted files optional for MAGeCK

## 1.3.2

Update RCRISPR to 1.1.1.0

## 1.3.1

Teething issue with first tagged CI build

## 1.3.0

First public release.  Linked to [RCRISPR](https://github.com/cancerit/RCRISPR) 1.1.0.0.

## 1.2.0

- Update functions for nnmd and glass delta
- Round glass and nnmd to 3dp
- Update to ubuntu:20.04 due to R issues
  - Additionally improve R package installs by adding versions
