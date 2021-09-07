# Example data

## Count matrix

Count matrix (`reads_hap1.txt`) is from the hart-lab BAGEL repository (commit f9eedca):

[https://github.com/hart-lab/bagel/blob/master/reads_hap1.txt](https://github.com/hart-lab/bagel/blob/master/reads_hap1.txt)

To download:

```
curl -OJL https://raw.githubusercontent.com/hart-lab/bagel/master/reads_hap1.txt
```

## Library

TKOv3 can be downloaded from Addgene:

https://media.addgene.org/cms/filer_public/71/a8/71a81179-7a62-4d75-9b53-236e6f6b7d4d/tkov3_guide_sequence.xlsx

Library will need to be opened and saved as a TSV.

## Commands to generate formatted inputs

Library (`TKOv3.tsv`) and counts (`reads_hap1.txt`) are joined by guide sequence (unique) and then formatted for use with the pipeline.

```
library=$(join -1 2 -2 1 -o 1.1,1.2,1.3,2.3,2.4,2.5,2.6 <(tail -n +2 TKOv3.tsv | sort -k2) <(tail -n +2 reads_hap1.txt | sort -k1))
awk -F" " 'BEGIN{OFS="\t"; print "sgRNA_ID","Gene","chr","start","end" } { split( $3,c,"[:|_|\-]"); gsub("chr","",c[1]); print $3,$1,c[1],c[2],c[3] }' <(echo "$library") > library.tsv
awk -F" " 'BEGIN{OFS="\t"; print "sgRNA_ID","Gene","HAP1_T0" } { split( $3,c,"[:|_|\-]"); gsub("chr","",c[1]); print $3,$1,$4 }' <(echo "$library") > 'HAP1_T0.tsv'
awk -F" " 'BEGIN{OFS="\t"; print "sgRNA_ID","Gene","HAP1_T18A" } { split( $3,c,"[:|_|\-]"); gsub("chr","",c[1]); print $3,$1,$5 }' <(echo "$library") > 'HAP1_T18A.tsv'
awk -F" " 'BEGIN{OFS="\t"; print "sgRNA_ID","Gene","HAP1_T18B" } { split( $3,c,"[:|_|\-]"); gsub("chr","",c[1]); print $3,$1,$6 }' <(echo "$library") > 'HAP1_T18B.tsv'
awk -F" " 'BEGIN{OFS="\t"; print "sgRNA_ID","Gene","HAP1_T18C" } { split( $3,c,"[:|_|\-]"); gsub("chr","",c[1]); print $3,$1,$7 }' <(echo "$library") > 'HAP1_T18C.tsv'
```

## Sample information

Manually generated based on the count matrix and BAGEL manuscript information. Raw read totals not provided.
