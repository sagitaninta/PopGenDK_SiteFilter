# Site Filtering Pipeline

This is a (experimental) documentation on the pipeline used by @popgenDK to analyse whole genome sequences of the non-model species. This pipeline starts at the fastqc files and the subsequent depth will be different depending on the mean depth of the genome. The goal of this pipeline is to produce a BED file that contains all sites to be used for genotype calling and further downstream analyses such as population structure, ROH, D-statistics, and others.

The pipeline consists of three major steps:
- mapping 
- reference quality filtering
- sample quality filtering

## 1. Mapping 

### a. Pre-mapping QC

Run MultiQC on the freshly arrived `.fastq` files.

### b. Mapping

Mapping is usually done to two different reference genomes: distantly related reference genome and closely related reference genome. This allow us to measure the error rate. 

For mapping, we run the fastqc with [`PALEOMIX`](https://github.com/MikkelSchubert/paleomix).

### c. Post-mapping QC

The resulting BAM alignments were filtered to remove unmapped reads, secondary alignments, PCR duplicates, and supplementary alignments, and reads flagged as having failed QC. We furthermore removed alignments with an inferred insert size <50 bp or >1000 bp and reads where less than 50 bp or 50% were mapped to the reference genome. Finally, we removed pairs of reads mapping to different contigs or in an unexpected orientation and reads for which the mate had been removed by any of the above criteria.

## 2. Reference quality filtering

### a. Mappability and RepeatMasker
We estimated the mappability of the reference genome using GENMAP v1.3.072. Here, we used 100 bp k-mers allowing for two mismatches (-K 100 -E 2) and the remaining parameters set to default settings. All sites with a mappability score <1 were excluded from downstream analyses. RepeatMasker v4.1.173 was used to identify repeat elements in the warthog genome assembly, utilising ‘rmblast’ as the search engine and ‘mammal’ as the query species with default settings. Repeat regions identified with RepeatMasker were masked to limit mismapping in these regions. Annotated sex chromosomes and scaffolds that were not assembled into chromosomes were also excluded.

### b. Global sequencing depth

Finally, we removed sites with extreme depth. We estimated the global depth (read count across all samples) for each site using Angsd74 (-minMapQ 30 -minQ 30 -doCounts 1 -doDepth 1 -dumpCounts 1 -maxdepth 4000). This was done separately for each species for all (n =67), unrelated (n = 54) and medium-high depth samples (n =18). Only autosomal chromosomes were included. From the global depth we calculated the upper 1% and lower 3% percentiles and visually inspected the plots before deciding on a threshold for excluding sites with extreme sequencing depth. Only sites that were within the thresholds for both low- and medium-high depth samples were used in the downstream analyses.

### c. Excessive heterozygosity filter

Several factors related to genomic repeats, like the presence of copy number variation between the reference genome and the analyzed samples, can lead to mis-mapping of reads originating from multiple genome locations to a single location. Any differences in these locations will result in inferred sites with an excess of heterozygosity that can be used to identify and mask these regions. We computed a preliminary genotype likelihoods file for common polymorphic sites (MAF R 0.05 and SNP p value < 106) restricted to data with a base quality of at least 30.21 Using these genotype likelihoods as input to PCAngsd,19 we calculated per-site inbreeding coefficients (F), ranging from 1 where all samples are heterozygous to 1 where all samples are homozygous, and performed a Hardy-Weinberg equilibrium (HWE) likelihood ratio test accounting for population structure.71 To obtain individual allele frequencies in PCAngsd, three principal components that account for population structure were considered. Based on the per-site inbreeding coefficients, windows of 100kb around sites with significant excessive heterozygosity estimates (F < 0.95 and p value < 106) were excluded (Figure S5E). In addition, entire scaffolds where either 20% or more of their total sequence was removed or scaffolds with average F value across all sites F < 0.1 were excluded. (Pecnerova et al 2021)

We also removed genomic regions with unusually high heterozygosity to avoid mismapping artefacts driven by multimapping on paralogous and other repetitive regions. We first estimated genotype likelihoods for SNPs using ANGSD with the GATK model (-GL 2), minimum mapping quality of 30 (-minMapQ 30), a minimum base quality of 30 (-minQ 30), a p-value of 1e−6 to call SNPs (-snp_pval 1e−6) and kept only SNPs with minor allele frequency (MAF) > 0.05 (-minmaf 0.05). Genotype likelihoods were then used as input for PCAngsd’sper site Hardy-Weinberg equilibrium (HWE) test75, which estimates inbreeding coefficients (F), and a likelihood ratio test statistic (LRT) for evidence of deviation from HWE, while controlling for population structure. The PCAngsd MAP test75 was also used to select the optimal number of principal components in each case. Sites with F < −0.9 and LRT > 24 were subsequently removed as they may have been driven by mapping artefacts, and therefore all regions within 10 kb from such sites were also discarded. (Balboa et al 2024)

### d. Autosomal and sex-linked scaffold identification

This is required if your reference genome is not assigned to proper chromosomes. To identify sex-linked scaffolds, we first calculated the average depth for each scaffold for each sample normalized by the average depth of the five largest scaffolds. Based on the normalized depths, we then performed a principal component analysis (PCA) and identified two main clusters of samples, likely representing males and females (Figure S5B). Finally, considering that in leopards females have two X chromosomes and males have a single X chromosome, we designated and excluded putative sex-linked scaffolds where the mean difference was greater than 0.4 between the two clusters. In addition, we excluded all scaffolds with an average normalized depth greater than 1.1 or lower than 0.9 across all samples. This resulted in 23 scaffolds being identified as X chromosome-linked (Figure S5D) and excluded from further analyses.

The union of all masks of the reference genome (scaffold size, mappability, RepeatMasker, aligned depth, autosome identification, excessive heterozygosity), henceforth referred to as strictref. The strictref filter was applied to all analyses unless otherwise noted.

## 3. Sample quality filtering

We identified and excluded samples with high sequencing error rates based on the “perfect individual” approach76. The rationale behind this approach is that any sample in the dataset should have equal genetic distance to the outgroup and therefore samples with excess/deficit of derived alleles would be interpreted as errors. As the “perfect individual” we used a high depth individual ... . This sample was processed with ANGSD to create a consensus sequence (-doFasta 2), taking the most commonly observed base as the consensus (-doAncError 2) while setting the base quality to at least 30 (-minQ 30). We chose ... as an outgroup and mapped all samples to the consensus using BWA excluding sex chromosomes, the mitogenome, repeats, and sites with mappability <1. Individuals with high error rates (> 0.001) were removed from downstream analyses.

We then considered relatedness between samples, where we identified and removed potential relatives and duplicated samples using the methodology described in IBSRELATE77. First, we calculated the Site Allele Frequency (SAF) likelihood in ANGSD for each individual. We used the genotype likelihood-based approach assuming HWE (-doSaf 1). This pipeline need an ancestral reference (-anc), a minimum mapping quality of 30 (-minMapQ 30), a minimum base quality of 30 (-minQ 30), and the GATK method (-GL 2). Then, we inferred the two-dimensional site frequency spectra (2D-SFS) pairwise among all possible combinations of individuals. To limit computational time, we limited the number of sites surveyed to the first 50,000 sites. Based on the 2D-SFS, we calculated R0, R1 and KINGrobust kinship77,78, which can be used to identify close familiar relatives. For the analysis, we include certain set of samples in this analysis to account for potentially interspecies duplicates or mislabelled samples. We identified and removed an individual from each pair of first- and second-degree relatives.

## 4. Genotype calling

We used ANGSD21 to estimate genotype likelihoods using the genotype likelihood (GL) model from GATK (-gl 2),65 inferring major and minor from GL data (-doMajorminor 1), estimating allele frequencies from the GL data (-doMaf 1), and applying the strictref filter (-sites). We restricted the analysis to bases with base quality of at least 30. For SNP calling we used a likelihood ratio test (p value = 106) as implemented in ANGSD.21 We only included common alleles (MAF > 0.05).
For the five high-coverage samples, we called genotypes using bcftools (v.1.9).67 We called the genotypes per sample using the multiallelic caller methodology (-m). We disabled base alignment qualities (BAQ; -B) and restricted the analysis to base qualities of at least 30. In addition to the strictref filter (-T), we also excluded sites with a depth of coverage below 10 and heterozygous calls with less than 3 reads support for each allele using plugin setGT to reduce genotype calling errors. (Pecnerova)

Based on the mapped reads, we used bcftools v1.1371 to call genotypes for the medium-high depth individuals requiring minimum base quality 25, mapping quality 30, and ignoring contigs smaller than 100 kb. We used the ‘--per-sample-mF’ flag for the pileup, and the ‘--multiallelic-caller’ for the calling. From the initial calls, we filtered out sites based on mappability, heterozygosity, and depth outliers as described in the following section and removed all indels and multiallelic sites. Finally, any genotype call with less than 10 read depth was set as missing, along with heterozygous calls with less than two reads supporting the minor allele. (Balboa)




