# selecting only biallelic SNPs:
bcftools view --types snps -m 2 -M 2 --threads 4 input.vcf.gz -Oz -o output_biallelic.vcf.gz

# filter the file for genotype quality
bcftools filter -i 'FMT/GQ >= 20' output_biallelic.vcf.gz -o filtered_GQ20.vcf

#1.- Calculate missing data with plink v2.0:
plink2 --vcf filtered_GQ20.vcf --missing vcols=chrom,pos,nmiss,fmiss --allow-extra-chr --out missing.out

#2.- Calculate read depth with vcftools v0.1.16:
vcftools --gzvcf vcfile.vcf.gz --site-mean-depth --out RD

#3.- R script to visualize the results:
library(ggplot2)
       args <- commandArgs(trailingOnly = TRUE)
Sites <- args[1]
Ind <- args[2]
depth <- args[3]
Sites_T <- read.table(Sites, head=T, comment.char="") 
Ind_T <- read.table(Ind, head=T, comment.char="")
depth_T <- read.table(depth, head=T, comment.char="")
jpeg("Missing_site.jpg")
ggplot(Sites_T, aes(x=F_MISS, color="red")) + geom_histogram(fill="white", bins=100)
dev.off()
jpeg("Missing_ind.jpg")
ggplot(Ind_T, aes(x=F_MISS, color="blue")) + geom_histogram(fill="white", bins=100)
dev.off()
 
jpeg("Readdepth.jpg")
ggplot(depth_T, aes(x=MEAN_DEPTH, color="green")) + geom_histogram(fill="white", bins=100) # + xlim(0,20) maybe you want to play with xlim to visualize better the RD
dev.off()

#I saved this script in the file Plot_missing_RD.R and I’ll ran it like this:
Rscript Plot_missing_RD.R missing.out.vmiss missing.out.smiss RD.ldepth.mean

#4.-Once you decide which sites to remove you could do it like this with bcftools v1.16 :

bcftools view -i 'AVG(FMT/DP)>0.2 & AVG(FMT/DP)<1.5' vcfile.vcf.gz -Oz > vcffile_RD0_4_1_5.vcf.gz

#5.- For the missing data once you have decided about your cut offs. You could create the tables to filter like this:

#!/umesr/bin/env Rscript
library(data.table)
library(ggplot2)
 
args <- commandArgs(trailingOnly = TRUE)
 
Sites <- args[1]
Ind <- args[2]
 
sitescut <- args[3]
indcut <- args[4]
 
Sites_T <- read.table(Sites, head=T, comment.char="")
 
Ind_T <- read.table(Ind, head=T, comment.char="")
 
 Sites_sub <- subset(Sites_T, Sites_T$F_MISS >= sitescut)
 Ind_sub <- subset(Ind_T, Ind_T$F_MISS >= indcut)
 
Sites_sub <- Sites_sub[,c(1,2)]
Ind_sub <-Ind_sub[,c(1,2)]
 
write.table(Sites_sub, "Sites_to_remove.txt", col.names = F, sep = "\t", row.names = F, quote = F)
write.table(Ind_sub, "Ind_to_remove.txt", col.names = F, sep = "\t", row.names = F, quote = F)

#You could run it like this:

Rscript Create_tables.R missing.out.vmiss missing.out.smiss 0.9 0.95
#0.9 is the value to select the sites with more than 90% of missing data that we want to filter and 0.95 to select the individuals with more than 95% of missing data to filter.

#6-To remove the sites, we could do it like this:

bcftools view -T ^Sites_to_remove.txt vcffile_RD0_4_1_5.vcf.gz | bgzip -c > vcffile_RD0_2__1_5_misssites09.vcf.gz
bcftools view -S ^Ind_to_remove.txt vcffile_RD0_2__1_5_misssites09.vcf.gz | bgzip -c > vcffile_RD0_2__1_5_misssites09_missind095.vcf.gz

#And this last line to remove the  individuals To filter sites by allele frequency

bcftools view -i 'MAF > 0.1' vcffile_RD0_2__1_5_misssites09_missind095.vcf.gz | bgzip -c > vcffile_RD0_2__1_5_misssites09_missind095_MAF10.vcf.gz

# you could check MAF values like this and decide if to filter them or not
vcftools --gzvcf vcfile.vcf.gz --freq2 --out MAF --max-alleles 2

#Filter to keep only biallelic sites, ya está arriba
#bcftools view --types snps -m 2 -M 2 --threads 4 vcffile_RD0_2__1_5_misssites09_missind095_MAF10.vcf.gz -Oz -o vcffile_RD0_2__1_5_misssites09_missind095_MAF10_biallelic.vcf.gz

#7.- filter repetitive regions
#I used this software https://github.com/Ensembl/plant-scripts/tree/master/repeats

#I download it like this:
git clone https://github.com/Ensembl/plant-scripts.git

#I ran the program like this:
./Red2Ensembl.py --exe /pathtoscripts/plant-scripts/plant-scripts/lib/Red/bin/Red --cor 9 --msk_file Redmasked.fasta --bed_file Redmasked.bed REF.fasta OUTDIR

#Then you could filtrate the SNPs in repetitive regions like this:
bedtools intersect -a vcffile_RD0_2__1_5_misssites09_missind095_MAF10_biallelic.vcf.gz -b Redmasked.bed -sorted -g REF.bed -v


#REF.bed is a file that has the the name of the chomosomes and its length
#This file looks like this:
#Chr01  159217232
#Chr02  184765313
#Chr03  182121094
#Chr04  221926752
#Chr05  187413619
#The file could be generated like this:
samtools faidx REF.fasta (editado) 
cut -f 1,2 REF.fasta.fai > REF.bed




