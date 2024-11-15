#!/usr/bin/env bash

# Gather singleton count data for major SV resources for NRG SV review

WRKDIR=/Users/ryan/Desktop/Collins/Talkowski/misc/NRG_SV_review_Fall2020/NRG_SV_review_figures/NRG_SV_review_figure_data/public_vcfs/
cd $WRKDIR
for subdir in singleton_tables deletion_beds ref; do
  if ! [ -e $WRKDIR/$subdir ]; then
    mkdir $WRKDIR/$subdir
  fi
done



# Curate hg19 Gencode exons, genes, promoters, and UTRs
wget \
  -P $WRKDIR/ref/ \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
zcat $WRKDIR/ref/gencode.v19.annotation.gtf.gz \
| fgrep -w protein_coding \
| awk -v OFS="\t" '{ if ($3=="UTR") print $1, $4, $5 }' \
| sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | bgzip -c \
> $WRKDIR/ref/hg19.UTRs.bed.gz
tabix -p bed -f $WRKDIR/ref/hg19.UTRs.bed.gz
zcat $WRKDIR/ref/gencode.v19.annotation.gtf.gz \
| fgrep -w protein_coding \
| awk -v OFS="\t" '{ if ($3=="exon") print $1, $4, $5 }' \
| sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - \
| bedtools subtract -a - -b $WRKDIR/ref/hg19.UTRs.bed.gz | bgzip -c \
> $WRKDIR/ref/hg19.exons.bed.gz
tabix -p bed -f $WRKDIR/ref/hg19.exons.bed.gz
zcat $WRKDIR/ref/gencode.v19.annotation.gtf.gz \
| fgrep -w protein_coding \
| awk -v OFS="\t" '{ if ($3=="transcript") print $1, $4, $5 }' \
| sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | bgzip -c \
> $WRKDIR/ref/hg19.genes.bed.gz
tabix -p bed -f $WRKDIR/ref/hg19.genes.bed.gz
zcat $WRKDIR/ref/gencode.v19.annotation.gtf.gz \
| fgrep -w protein_coding \
| awk -v OFS="\t" '{ if ($3=="transcript") print $1, $4, $5, $7 }' \
| awk -v OFS="\t" -v pdist=1000 '{ if ($4=="+") print $1, $2-pdist, $2; else print $1, $3, $3+pdist }' \
| sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | bgzip -c \
> $WRKDIR/ref/hg19.promoters.bed.gz
tabix -p bed -f $WRKDIR/ref/hg19.promoters.bed.gz
# Reformat all hg19 files to GRCh37
for element in exons UTRs genes promoters; do
  zcat $WRKDIR/ref/hg19.$element.bed.gz \
  | sed 's/^chr//g' | bgzip -c \
  > $WRKDIR/ref/b37.$element.bed.gz
done



# Curate hg38 Gencode exons, genes, promoters, and UTRs
wget \
  -P $WRKDIR/ref/ \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
zcat $WRKDIR/ref/gencode.v44.annotation.gtf.gz \
| fgrep -w protein_coding \
| awk -v OFS="\t" '{ if ($3=="UTR") print $1, $4, $5 }' \
| sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | bgzip -c \
> $WRKDIR/ref/hg38.UTRs.bed.gz
tabix -p bed -f $WRKDIR/ref/hg38.UTRs.bed.gz
zcat $WRKDIR/ref/gencode.v44.annotation.gtf.gz \
| fgrep -w protein_coding \
| awk -v OFS="\t" '{ if ($3=="exon") print $1, $4, $5 }' \
| sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - \
| bedtools subtract -a - -b $WRKDIR/ref/hg38.UTRs.bed.gz | bgzip -c \
> $WRKDIR/ref/hg38.exons.bed.gz
tabix -p bed -f $WRKDIR/ref/hg38.exons.bed.gz
zcat $WRKDIR/ref/gencode.v44.annotation.gtf.gz \
| fgrep -w protein_coding \
| awk -v OFS="\t" '{ if ($3=="transcript") print $1, $4, $5 }' \
| sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | bgzip -c \
> $WRKDIR/ref/hg38.genes.bed.gz
tabix -p bed -f $WRKDIR/ref/hg38.genes.bed.gz
zcat $WRKDIR/ref/gencode.v44.annotation.gtf.gz \
| fgrep -w protein_coding \
| awk -v OFS="\t" '{ if ($3=="transcript") print $1, $4, $5, $7 }' \
| awk -v OFS="\t" -v pdist=1000 '{ if ($4=="+") print $1, $2-pdist, $2; else print $1, $3, $3+pdist }' \
| sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | bgzip -c \
> $WRKDIR/ref/hg38.promoters.bed.gz
tabix -p bed -f $WRKDIR/ref/hg38.promoters.bed.gz





# Prep header for output file
echo -e "#source\tsvtype\tcriteria\tN_all\tN_singleton\tN_polymorphic" \
> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.tsv





# gnomAD v2
vcf=$WRKDIR/gnomad_v2.1_sv.sites.vcf.gz
head -n1 $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.tsv \
> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.gnomad.tsv
# All SVs (baseline)
bcftools query \
  --format '%INFO/AC\n' \
  --regions $( seq 1 22 | paste -s -d, ) \
  --include "AC>0 & FILTER = \"PASS\"" \
  $vcf \
| awk -v OFS="\t" -v dset="gnomAD_v2.1" -v criteria="all" -v svtype="ALL" \
  '{ if($1==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.gnomad.tsv
# Per SV type
for SVTYPE in DEL DUP INS INV CPX; do
  bcftools query \
    --format '%INFO/AC\n' \
    --regions $( seq 1 22 | paste -s -d, ) \
    --include "AC>0 & INFO/SVTYPE=\"$SVTYPE\" & FILTER = \"PASS\"" \
    $vcf \
  | awk -v OFS="\t" -v dset="gnomAD_v2.1" -v criteria="all" -v svtype=$SVTYPE \
    '{ if($1==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }'
done \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.gnomad.tsv
# Deletions split by context
bcftools query \
  --format '%CHROM\t%POS\t%END\t%AC\n' \
  --regions $( seq 1 22 | paste -s -d, ) \
  --include "AC>0 & INFO/SVTYPE=\"DEL\" & FILTER = \"PASS\"" \
    $vcf \
| sort -Vk1,1 -k2,2n -k3,3n | bgzip -c \
> $WRKDIR/deletion_beds/gnomad.DEL_wAC.bed.gz
tabix -p bed -f $WRKDIR/deletion_beds/gnomad.DEL_wAC.bed.gz
# Exonic
bedtools intersect -wa -u \
  -a $WRKDIR/deletion_beds/gnomad.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/b37.exons.bed.gz \
| awk -v OFS="\t" -v dset="gnomAD_v2.1" -v criteria="exonic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.gnomad.tsv
# UTRs (excluding exons)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/gnomad.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/b37.exons.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/b37.UTRs.bed.gz \
| awk -v OFS="\t" -v dset="gnomAD_v2.1" -v criteria="utr" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.gnomad.tsv
# Introns (excluding exons & UTRs)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/gnomad.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/b37.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.UTRs.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/b37.genes.bed.gz \
| awk -v OFS="\t" -v dset="gnomAD_v2.1" -v criteria="intronic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.gnomad.tsv
# Promoter (excluding exon, UTR, and intron)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/gnomad.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/b37.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.UTRs.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.genes.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/b37.promoters.bed.gz \
| awk -v OFS="\t" -v dset="gnomAD_v2.1" -v criteria="promoter" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.gnomad.tsv
# Intergenic (excluding all of the above)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/gnomad.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/b37.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.UTRs.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.genes.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.promoters.bed.gz \
| awk -v OFS="\t" -v dset="gnomAD_v2.1" -v criteria="intergenic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.gnomad.tsv





# 1000G phase 3
vcf=$WRKDIR/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
head -n1 $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.tsv \
> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000GPh3.tsv
# All SVs (baseline)
bcftools view \
  --include "FILTER = \"PASS\" | FILTER = \".\"" \
  -m2 -M2 \
  --regions $( seq 1 22 | paste -s -d, ) \
  $vcf \
| bcftools query \
  --format '%INFO/AC\n' \
  --include "AC>0" \
  - \
| awk -v OFS="\t" -v dset="1000GPh3" -v criteria="all" -v svtype="ALL" \
  '{ if($1==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000GPh3.tsv
# Per SV type
for SVTYPE in DEL DUP INV; do
  bcftools view \
    --include "FILTER = \"PASS\" | FILTER = \".\"" \
    -m2 -M2 \
    --regions $( seq 1 22 | paste -s -d, ) \
    $vcf \
  | bcftools query \
    --format '%INFO/AC\n' \
    --include "AC>0 & INFO/SVTYPE ~ \"$SVTYPE\"" \
    - \
  | awk -v OFS="\t" -v dset="1000GPh3" -v criteria="all" -v svtype=$SVTYPE \
    '{ if($1==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }'
done \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000GPh3.tsv
# Need to do insertions separately
bcftools view \
  --include "FILTER = \"PASS\" | FILTER = \".\"" \
  -m2 -M2 \
  --regions $( seq 1 22 | paste -s -d, ) \
  $vcf \
| bcftools view \
  --include "INFO/SVTYPE = \"ALU\" | INFO/SVTYPE = \"INS\" | INFO/SVTYPE = \"LINE1\" | INFO/SVTYPE = \"SVA\"" \
  - \
| bcftools query \
  --format '%INFO/AC\n' \
  --include "AC>0" \
  - \
| awk -v OFS="\t" -v dset="1000GPh3" -v criteria="all" -v svtype="INS" \
  '{ if($1==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000GPh3.tsv
# Deletions split by context
bcftools view \
  --include "AC>0 & INFO/SVTYPE ~ \"DEL\"" \
  --regions $( seq 1 22 | paste -s -d, ) \
  -m2 -M2 \
  $vcf \
| bcftools query \
  --format '%CHROM\t%POS\t%END\t%AC\n' \
  --include "FILTER = \"PASS\" | FILTER = \".\"" \
  - \
| sort -Vk1,1 -k2,2n -k3,3n | bgzip -c \
> $WRKDIR/deletion_beds/1000GPh3.DEL_wAC.bed.gz
tabix -p bed -f $WRKDIR/deletion_beds/1000GPh3.DEL_wAC.bed.gz
# Exonic
bedtools intersect -wa -u \
  -a $WRKDIR/deletion_beds/1000GPh3.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/b37.exons.bed.gz \
| awk -v OFS="\t" -v dset="1000GPh3" -v criteria="exonic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000GPh3.tsv
# UTRs (excluding exons)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/1000GPh3.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/b37.exons.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/b37.UTRs.bed.gz \
| awk -v OFS="\t" -v dset="1000GPh3" -v criteria="utr" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000GPh3.tsv
# Introns (excluding exons & UTRs)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/1000GPh3.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/b37.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.UTRs.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/b37.genes.bed.gz \
| awk -v OFS="\t" -v dset="1000GPh3" -v criteria="intronic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000GPh3.tsv
# Promoter (excluding exon, UTR, and intron)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/1000GPh3.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/b37.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.UTRs.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.genes.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/b37.promoters.bed.gz \
| awk -v OFS="\t" -v dset="1000GPh3" -v criteria="promoter" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000GPh3.tsv
# Intergenic (excluding all of the above)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/1000GPh3.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/b37.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.UTRs.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.genes.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/b37.promoters.bed.gz \
| awk -v OFS="\t" -v dset="1000GPh3" -v criteria="intergenic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000GPh3.tsv





# 1000G 30X
vcf=$WRKDIR/HGSV.WGS.wAFs.sites.vcf.gz
head -n1 $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.tsv \
> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000G30X.tsv
# All SVs (baseline)
bcftools query \
  --format '%INFO/AC\n' \
  --regions $( seq 1 22 | awk '{ print "chr"$1 }' | paste -s -d, ) \
  --include "AC>0 & FILTER = \"PASS\"" \
  $vcf \
| awk -v OFS="\t" -v dset="1000G30X" -v criteria="all" -v svtype="ALL" \
  '{ if($1==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000G30X.tsv
# Per SV type
for SVTYPE in DEL DUP INV CPX; do
  bcftools query \
    --format '%INFO/AC\n' \
    --regions $( seq 1 22 | awk '{ print "chr"$1 }' | paste -s -d, ) \
    --include "AC>0 & INFO/SVTYPE=\"$SVTYPE\" & FILTER = \"PASS\"" \
    $vcf \
  | awk -v OFS="\t" -v dset="1000G30X" -v criteria="all" -v svtype=$SVTYPE \
    '{ if($1==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }'
done \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000G30X.tsv
# Must do insertions separately
bcftools view \
  --regions $( seq 1 22 | awk '{ print "chr"$1 }' | paste -s -d, ) \
  --include "INFO/SVTYPE = \"INS\" | INFO/SVTYPE = \"MEI\"" \
  $vcf \
| bcftools query \
  --format '%INFO/AC\n' \
  --include "AC>0 & FILTER = \"PASS\"" \
  - \
| awk -v OFS="\t" -v dset="1000G30X" -v criteria="all" -v svtype="INS" \
  '{ if($1==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000G30X.tsv
# Deletions split by context
bcftools query \
  --format '%CHROM\t%POS\t%END\t%AC\n' \
  --regions $( seq 1 22 | awk '{ print "chr"$1 }' | paste -s -d, ) \
  --include "AC>0 & INFO/SVTYPE=\"DEL\" & FILTER = \"PASS\"" \
    $vcf \
| sort -Vk1,1 -k2,2n -k3,3n | bgzip -c \
> $WRKDIR/deletion_beds/1000G30X.DEL_wAC.bed.gz
tabix -p bed -f $WRKDIR/deletion_beds/1000G30X.DEL_wAC.bed.gz
# Exonic
bedtools intersect -wa -u \
  -a $WRKDIR/deletion_beds/1000G30X.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/hg38.exons.bed.gz \
| awk -v OFS="\t" -v dset="1000G30X" -v criteria="exonic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000G30X.tsv
# UTRs (excluding exons)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/1000G30X.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/hg38.exons.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/hg38.UTRs.bed.gz \
| awk -v OFS="\t" -v dset="1000G30X" -v criteria="utr" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000G30X.tsv
# Introns (excluding exons & UTRs)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/1000G30X.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/hg38.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.UTRs.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/hg38.genes.bed.gz \
| awk -v OFS="\t" -v dset="1000G30X" -v criteria="intronic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000G30X.tsv
# Promoter (excluding exon, UTR, and intron)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/1000G30X.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/hg38.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.UTRs.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.genes.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/hg38.promoters.bed.gz \
| awk -v OFS="\t" -v dset="1000G30X" -v criteria="promoter" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000G30X.tsv
# Intergenic (excluding all of the above)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/1000G30X.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/hg38.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.UTRs.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.genes.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.promoters.bed.gz \
| awk -v OFS="\t" -v dset="1000G30X" -v criteria="intergenic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.1000G30X.tsv





# CCDG
vcf=$WRKDIR/Abel_2020.Build38.public.v2.vcf.gz
head -n1 $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.tsv \
> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.CCDG.tsv
# All SVs (baseline)
bcftools query \
  --format '%INFO/AC\n' \
  --regions $( seq 1 22 | awk '{ print "chr"$1 }' | paste -s -d, ) \
  --include "AC>0 & FILTER = \"PASS\" & INFO/SVTYPE !=\"BND\"" \
  $vcf \
| awk -v OFS="\t" -v dset="CCDG" -v criteria="all" -v svtype="ALL" \
  '{ if($1==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.CCDG.tsv
# Per SV type
for SVTYPE in DEL DUP MEI INV; do
  if [ $SVTYPE == "MEI" ]; then
    svt_label="INS"
  else
    svt_label=$SVTYPE
  fi
  bcftools query \
    --format '%INFO/AC\n' \
    --regions $( seq 1 22 | awk '{ print "chr"$1 }' | paste -s -d, ) \
    --include "AC>0 & INFO/SVTYPE=\"$SVTYPE\" & FILTER = \"PASS\"" \
    $vcf \
  | awk -v OFS="\t" -v dset="CCDG" -v criteria="all" -v svtype=$svt_label \
    '{ if($1==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }'
done \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.CCDG.tsv
# Deletions split by context
bcftools query \
  --format '%CHROM\t%POS\t%END\t%AC\n' \
  --regions $( seq 1 22 | awk '{ print "chr"$1 }' | paste -s -d, ) \
  --include "AC>0 & INFO/SVTYPE=\"DEL\" & FILTER = \"PASS\"" \
    $vcf \
| sort -Vk1,1 -k2,2n -k3,3n | bgzip -c \
> $WRKDIR/deletion_beds/CCDG.DEL_wAC.bed.gz
tabix -p bed -f $WRKDIR/deletion_beds/CCDG.DEL_wAC.bed.gz
# Exonic
bedtools intersect -wa -u \
  -a $WRKDIR/deletion_beds/CCDG.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/hg38.exons.bed.gz \
| awk -v OFS="\t" -v dset="CCDG" -v criteria="exonic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.CCDG.tsv
# UTRs (excluding exons)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/CCDG.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/hg38.exons.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/hg38.UTRs.bed.gz \
| awk -v OFS="\t" -v dset="CCDG" -v criteria="utr" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.CCDG.tsv
# Introns (excluding exons & UTRs)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/CCDG.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/hg38.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.UTRs.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/hg38.genes.bed.gz \
| awk -v OFS="\t" -v dset="CCDG" -v criteria="intronic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.CCDG.tsv
# Promoter (excluding exon, UTR, and intron)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/CCDG.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/hg38.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.UTRs.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.genes.bed.gz \
| bedtools intersect -wa -u -a - -b $WRKDIR/ref/hg38.promoters.bed.gz \
| awk -v OFS="\t" -v dset="CCDG" -v criteria="promoter" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.CCDG.tsv
# Intergenic (excluding all of the above)
bedtools intersect -wa -v \
  -a $WRKDIR/deletion_beds/CCDG.DEL_wAC.bed.gz \
  -b $WRKDIR/ref/hg38.exons.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.UTRs.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.genes.bed.gz \
| bedtools intersect -wa -v -a - -b $WRKDIR/ref/hg38.promoters.bed.gz \
| awk -v OFS="\t" -v dset="CCDG" -v criteria="intergenic" -v svtype="DEL" \
  '{ if($4==1){SINGLE+=1}else{POLY+=1} }END{ print dset, svtype, criteria, SINGLE + POLY, SINGLE, POLY }' \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.CCDG.tsv



# Combine all four studies into single table
cat $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.*.tsv \
| sort -Vk2,2 -k1,1V -k3,3V | uniq | fgrep -v "#" \
>> $WRKDIR/singleton_tables/NRG_SV_review.singleton_counts.tsv