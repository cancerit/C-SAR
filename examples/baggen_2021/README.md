# Example dataset: Baggen et al. 2021

In 2021, Baggen et al. performed a genome-wide CRISPR screen on SARS-CoV-2 infected human Huh7 cells, transfected with the Brunello lentiCRISPRv2 library to identify host factors for infection.

> Baggen, J., Persoons, L., Vanstreels, E. et al.
> Nat Genet 53, 435–444 (2021).
> https://doi.org/10.1038/s41588-021-00805-2

## Analysing the example dataset with C-SAR

### Clone the C-SAR GitHub repository

To run [C-SAR](https://github.com/cancerit/C-SAR) on the example dataset, you will first need to clone the C-SAR repository:

```
git clone --recurse-submodules --branch main https://github.com/cancerit/C-SAR.git
```

### Set path for example directory

To make the analysis location agnostic, we suggest setting an environment variable pointing at the `examples/baggen_2021` directory.

```
ANALYSIS_DIR=<path to example directory>
```

### Running C-SAR

[C-SAR](https://github.com/cancerit/C-SAR) can be driven from a single configuration file `baggen_low_stringency.config` which contains the path to the individual counts (`reformatted_counts`), sample metadata (`baggen_low_stringency_manifest.tsv`), library (`Brunello_library_Baggen_SupplTable3.reformatted.csar.tsv`) and essential (`essential_genes.txt`)/non-essential (`non_essential_genes.txt`) input files as well as parameters for each of the analysis stages.

Once installed, to run C-SAR:

```
c-sar -c ${ANALYSIS_DIR}/baggen_low_stringency.config -w ${ANALYSIS_DIR}/work
```

*Note: an internal C-SAR version (1.2.0) was used for CRISPR and Beyond: Perturbations at Scale to Underatand Genomes 2021 (available on request).*

## How was the data prepared?

### Download FASTQ from SRA

FASTQ files for the relevant samples were downloaded from the SRA using [SRA-toolkit version 2.11.0](https://github.com/ncbi/sra-tools/archive/refs/tags/2.11.0.tar.gz).

Once installed, to download the relevant FASTQ files (listed in `sra-accessions.csv`):

```
mkdir ${ANALYSIS_DIR}/fastq
while read acc; do
  echo "Downloading $acc"
  fastq-dump --gzip -O ${ANALYSIS_DIR}/fastq $acc
done <${ANALYSIS_DIR}/sra_accessions.csv
```

This will write gzipped FASTQ files into a directory called `fastq`.

### Brunello lentiCRISPRv2 library

Brunello lentiCRISPRv2 library containing 76,441 guides was downloaded from [AddGene](https://www.addgene.org/static/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt).

For the library to be used with [pyCROQUET](https://github.com/cancerit/pycroquet) and [C-SAR](https://github.com/cancerit/C-SAR), the library needs to be reformatted.

#### pyCROQUET

The library used with pyCROQUET is `Brunello_library_Baggen_SupplTable3.reformatted.tsv` which has the following format:

```
##library-type: single
##library-name: Brunello_Baggen_SupplTable3
##species: human
##assembly: ???
##gene-build-source: ???
##gene-build-version: ???
#id	sgrna_ids	sgrna_seqs	gene_pair_id
1	A1BG_1	CATCTTCTTTCACCTGAACG	A1BG
2	A1BG_2	CTCCGGGGAGAACTCCGGCG	A1BG
3	A1BG_3	TCTCCATGGTGCATCAGCAC	A1BG
4	A1BG_4	TGGAAGTCCACTCCACTCAG	A1BG
5	A2M_1	ACTGCATCTGTGCAAACGGG	A2M
```

#### C-SAR

The library used with C-SAR is `Brunello_library_Baggen_SupplTable3.reformatted.csar.tsv` which has the following format:

```
id	sgrna_ids	sgrna_seqs	gene_pair_id
1	A1BG_1	CATCTTCTTTCACCTGAACG	A1BG
2	A1BG_2	CTCCGGGGAGAACTCCGGCG	A1BG
3	A1BG_3	TCTCCATGGTGCATCAGCAC	A1BG
4	A1BG_4	TGGAAGTCCACTCCACTCAG	A1BG
5	A2M_1	ACTGCATCTGTGCAAACGGG	A2M
```

It's essentially the pyCROQUET-formatted library with the headers stripped to leave a single header row. In the C-SAR configuration file  (`baggen_low_stringency.config`) we use `library_id_column_index = 2` and `library_gene_column_index = 4` to tell C-SAR that the guide identifiers and genes are in column 2 and 4 respectively.

### Quantification with pyCROQUET

#### Formatting FASTQ for use with pyCROQUET

A pre-release version of [pyCROQUET](https://github.com/cancerit/pycroquet) was used to quantify each sample (i.e. the number of reads assigend to each guide). However, [version 1.1.1](https://github.com/cancerit/pycroquet/releases/tag/1.1.1) can be used to replicate with functionality remaining unchanged.

Pre-quantification, a small modification to the FASTQ read headers is required (if they are not standard Illumina headers):

```
mkdir ${ANALYSIS_DIR}/reformatted_fastqs
for i in ${ANALYSIS_DIR}/fastq/*; do
    f=$(basename "$i")
    zcat $i | sed "s/ /:/g" | sed "s/length=151/0000:00000:00000/g" > "${ANALYSIS_DIR}reformatted_fastq/$f";
done
```

This will write the reformatted FASTQ files to `reformatted_fastq`.

#### Quantifying samples with pyCROQUET

Below is an example pyCROQUET command for SRR13255539:

```
pycroquet single-guide -g ${ANALYSIS_DIR}/Brunello_library_Baggen_SupplTable3.reformatted.tsv -q ${ANALYSIS_DIR}/reformatted_fastq/SRR13255539.fastq.gz -s SRR13255539 -o ${ANALYSIS_DIR}/counts/SRR13255539 -n"
```

To run with multiple CPU, add `-c <cpu>`> to the command above e.g. `-c 4` to use 4 CPUs.

### Essential and non-essential gene lists

BAGEL essential [CEG](https://github.com/hart-lab/bagel/blob/master/CEGv2.txt) and non-essential [NEG](https://github.com/hart-lab/bagel/blob/master/NEGv1.txt) gene lists were downloaded for mouse from the [BAGEL2 repository](https://github.com/hart-lab/bagel) (commit ref: 53388ad) as `essential_genes.txt` and `non_essential_genes.txt` respectively.
