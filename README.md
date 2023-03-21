# jVCFparser
### A command-line parser for VCF files designed for population genetics analyses.

This is a beta version of the jVCFparser command line tool for processing variant call format (VCF) files. The tool uses memory-efficient descriptive statistics to load VCF file data into memory and then perform population genetics calculations on it. Since the tool only stores allele and genotype frequencies, it is able to process large files. Although reading the files may take some time, all calculations are extremely fast. The tool has been tested on VCF versions 4.0 and 4.2.

VCF 4.2 description (23.08.2022): The manual can be accessed on [SAMtools site](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

Sample data used for testing:

~25K SNP (loci) and 180 sample matrix: [Sessile oak SNP dataset; de novo assembly](https://zenodo.org/record/3908963); File size: ~15MB<br>
~50K SNP (loci) and ~20K sample matrix: [SoySNP50K iSelect BeadChip, Wm82.a1](https://soybase.org/snps/); File size: ~3.5GB<br>
~1.3M SNP (loci) and ~2K sample matrix: [1000 genomes project, Phase 3, Chromosome 21](http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/); File size: ~10GB<br>

Allele and genotype counts for reading the files are the following:<br>
- Number of reference allele '0'
- Number of alternative allele '1'
- Number of unique alleles
- Number of homozygote genotypes (e.g. 0/0 or 1/1)
- Number of heterozygote genotypes (e.g. 0/1 or 1/0)
- Number of missing genotypes (e.g. ./.)
- Number of unique genotypes

Currently implemented diversity-ralated descriptive statistics, and calculations (and counts) as follows:<br>
- Number of missing genotypes
- Number of REF alleles
- Number of ALT alleles
- Number of genotypes
- Number of heterozygotes
- Number of homozygotes
- Average number of different genotypes
- Average number of different alleles (Na)
- Average number of effective alleles (Ne)
- Average Observed Heterozygosity (Ho)
- Average Expected Heterozygosity (He)
- Average Unbiased Expected Heterozygosity (uHe)
- Average Shannon's Information Index (SI)
- Average Simpson's Diversity Index (D)
- Average Fixation Index (F)
- Average Allelic Richness (Ar)

