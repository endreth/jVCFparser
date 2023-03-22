# jVCFparser
### A command-line parser for VCF files designed for population genetics analyses.

This is a beta version of the jVCFparser command line tool for processing variant call format (VCF) files. The tool uses memory-efficient descriptive statistics to load VCF file data into memory and then perform population genetics calculations on it. Since the tool only stores allele and genotype frequencies, it is able to process large files. Although reading the files may take some time, all calculations are extremely fast. The tool has been tested on VCF versions 4.0 and 4.2.

VCF 4.2 description (23.08.2022): The manual can be accessed on [SAMtools site](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

### Usage (example)

Get the JAR artifact [HERE](https://github.com/endreth/jVCFparser/blob/7dfacf604013e732b53484171909c0b713c56537/out/artifacts/jVCFparser_jar/jVCFparser.jar)!
```
$ java -jar jVCFparser.jar -f ".\populations.snps.vcf" -agc -oh -eh
```
![image](https://user-images.githubusercontent.com/104054427/226903182-23d5e9aa-1b05-4381-9647-75c658f05bb8.png)

<b>Requirements:</b><br>
GNU/Linux, Microsoft Windows, or macOS<br>
JRE (JDK 11 or later)<br>

<b>Sample data used for testing:</b><br>
~25K SNP (loci) and 180 sample matrix: [Sessile oak SNP dataset; de novo assembly](https://zenodo.org/record/3908963); File size: ~15MB<br>
~50K SNP (loci) and ~20K sample matrix: [SoySNP50K iSelect BeadChip, Wm82.a1](https://soybase.org/snps/); File size: ~3.5GB<br>
~1.3M SNP (loci) and ~2K sample matrix: [1000 genomes project, Phase 3, Chromosome 21](http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/); File size: ~10GB<br>

<b>Allele and genotype counts for reading the files are the following:</b><br>
- Number of reference allele '0'
- Number of alternative allele '1'
- Number of unique alleles
- Number of homozygote genotypes (e.g. 0/0 or 1/1)
- Number of heterozygote genotypes (e.g. 0/1 or 1/0)
- Number of missing genotypes (e.g. ./.)
- Number of unique genotypes

<b>Currently implemented diversity-ralated descriptive statistics, and calculations (and counts) as follows:</b><br>
- Number of missing genotypes
- Number of REF alleles
- Number of ALT alleles
- Number of genotypes
- Number of heterozygotes
- Number of homozygotes
- Average number of different genotypes (Ng)
- Average number of different alleles (Na)
- Average number of effective alleles (Ne)
- Average Observed Heterozygosity (Ho)
- Average Expected Heterozygosity (He)
- Average Unbiased Expected Heterozygosity (uHe)
- Average Shannon's Information Index (SI)
- Average Simpson's Diversity Index (D)
- Average Fixation Index (F)
- Average Allelic Richness (Ar)

<details open>
<summary><b>Calculation details</b> [Formulas used in calculations and their references.]</summary><br>

<b>Average number of different genotypes (Ng):</b>

![avg_ng](https://user-images.githubusercontent.com/104054427/226822434-6e95eb46-cad5-439d-93c3-812ac513f748.png)<br>
N<sub>g</sub> represents the mean number of distinct genotypes across n loci, denoted as g<sub>i</sub> for i = 1,2,...,n. It is computed as the sum of the distinct genotypes at each locus, divided by the total number of loci.<br>

<b>Average number of different alleles (Na):</b>

![avg_na](https://user-images.githubusercontent.com/104054427/226824848-e233419e-a937-470f-a109-bcb61ebb7e9d.png)<br>
N<sub>a</sub> represents the mean number of different alleles across n loci, denoted as a<sub>i</sub> for i = 1,2,...,n. It is computed as the sum of the different alleles at each locus, divided by the total number of loci.<br>

<b>Average number of effective alleles (Ne):</b>

![avg_ne](https://user-images.githubusercontent.com/104054427/226825215-6d3458bb-a2ee-42af-a5c7-d59f426bb27c.png)<br>
N<sub>e</sub> represents the mean number of effective alleles across n genetic loci, denoted as p<sub>i</sub> for i = 1,2,...,n. It is calculated as the inverse of the sum of allele frequencies, divided by the total number of loci. Based on Brown & Weir (1983).<br>

<b>Average Observed Heterozygosity (Ho):</b>

![ho](https://user-images.githubusercontent.com/104054427/226826154-0d24f49c-cd02-4cca-9b80-7beb7b653f2a.png)<br>
H<sub>o</sub> represents the average proportion of heterozygous individuals across n genetic loci. For each locus i = 1,2,...,n, the proportion of heterozygous individuals is computed as the ratio of the number of heterozygotes to the total number of individuals N. Observed heterozygosity is then calculated as the mean of these proportions across all n loci. Based on Hartl & Clark (1997).<br>

<b>Average Expected Heterozygosity (He):</b>

![he](https://user-images.githubusercontent.com/104054427/226826977-b14e2c88-7a60-4c4b-9758-f5d2c284d443.png)<br>
H<sub>e</sub> represents the mean probability that two randomly chosen alleles at a given locus are different, across n genetic loci, denoted as p<sub>i</sub> and q<sub>i</sub> for i = 1,2,...,n. It is calculated as the average of 1 minus the sum of squared allele frequencies. Based on the intra locus gene diversity (H = 1-p<sup>2</sup>-q<sup>2</sup>) derived from the Hardly-Weinberg equilibrium.<br>

<b>Average Unbiased Expected Heterozygosity (uHe):</b>

![uhe](https://user-images.githubusercontent.com/104054427/226839861-b24fbf41-ff12-4b50-aad7-bb5569a9e186.png)<br>
uH<sub>e</sub> represents the mean probability that two randomly chosen alleles at a given locus are different, across n genetic loci, adjusted for sample size and population size bias. It is calculated as the average of 1 minus the sum of squared allele frequencies, adjusted for sample size bias. Based on Peakall & Smouse (2006).<br>

<b>Average Shannon's Information Index (H):</b>

![shannon](https://user-images.githubusercontent.com/104054427/226840623-6fb903e9-9fc4-4596-b692-ee5324240b1c.png)<br>
H represents the Average Shannon's Information Index, defined as the average amount of uncertainty associated with predicting the identity of a randomly chosen allele at a given locus, across n loci. It is calculated as the negative average of the product of the frequency of the i-th allele, p<sub>i</sub>, and the natural logarithm of p<sub>i</sub>. Based on Brown & Weir (1983).<br>

<b>Average Simpson's Diversity Index (D):</b>

![simpson](https://user-images.githubusercontent.com/104054427/226841327-a3ba0679-e77e-40f0-b60a-fadd2132e0e3.png)<br>
D represents the Average Simpson's Diversity Index, defined as the probability that two randomly chosen alleles at a given locus are identical, across n loci. It is calculated as 1 minus the average of the sum of squared allele frequencies. Based on Simpson (1949) and Morris et al. (2014).<br>

<b>Average Fixation Index (F):</b>

![f](https://user-images.githubusercontent.com/104054427/226841959-2d27395c-94f9-4fd4-a964-bbd552b5d38c.png)<br>
F represents the Average Fixation Index, averaged across n loci. It is calculated as the difference between observed heterozygosity (Ho) and expected heterozygosity (He), normalized by expected heterozygosity and averaged across n loci. Based on Hartl & Clark (1997).<br>

<b>Average Allelic Richness (Ar):</b>

![ar](https://user-images.githubusercontent.com/104054427/226842779-bcde0246-d58d-4226-9c0e-eea1b6014b5d.png)<br>
Ar represents the Average Allelic Richness, defined as the expected number of species in a sample of n genotypes selected at random from a collection containing N alleles ("genes") from S loci. It is calculated as the number of alleles observed in a sample of size N<sub>i</sub>, normalized by the sample size N<sub>i</sub> and averaged across S loci. Based on Hurlbert (1971) and El Mousadik & Petit (1996).<br>

</details>

<details>
<summary><b>References</b></summary><br>
Brown, A. H., & Weir, B. S. (1983). Measuring genetic variability in plant populations. Isozymes in plant genetics and breeding, part A, 219-239.<br><br>
El Mousadik, A., & Petit, R. J. (1996). High level of genetic differentiation for allelic richness among populations of the argan tree [Argania spinosa (L.) Skeels] endemic to Morocco. Theoretical and applied genetics, 92, 832-839.<br><br>
Hartl, D. L., & Clark, A. G. (1997). Principles of population genetics (Vol. 116). Sunderland: Sinauer associates.<br><br>
Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique and alternative parameters. Ecology, 52(4), 577-586.<br><br>
Morris, E. K., Caruso, T., Buscot, F., Fischer, M., Hancock, C., Maier, T. S., ... & Rillig, M. C. (2014). Choosing and using diversity indices: insights for ecological applications from the German Biodiversity Exploratories. Ecology and evolution, 4(18), 3514-3524.<br><br>
Peakall, R. O. D., and Peter E. Smouse. "GENALEX 6: genetic analysis in Excel. Population genetic software for teaching and research." Molecular ecology notes 6.1 (2006): 288-295.<br><br>
Simpson, E. H. (1949). Measurement of diversity. nature, 163(4148), 688-688.<br><br>

</details>
