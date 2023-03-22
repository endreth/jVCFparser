package hu.jvcfparser;

import java.io.*;
import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Stream;
import org.apache.commons.math3.util.CombinatoricsUtils;

public class Main {
    public static void main(String[] args) {

        //readVCF("C:\\Users\\thend\\Desktop\\jVCFparser\\soysnp50k_wm82.a1_41317.vcf");
        //readVCF("C:\\Users\\thend\\Desktop\\jVCFparserBeta\\ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf");
        //readVCF("C:\\Users\\thend\\Desktop\\jVCFparserBeta\\populations.snps.vcf");
        Map vcf = readVCF("C:\\Users\\thend\\Desktop\\jVCFparser\\FilesToTestCalculations\\populations.snps.vcf");

        numOfMissingGT(vcf);
        numOfREFAllele(vcf);
        numOfALTAllele(vcf);
        numOfGenotypes(vcf);
        numOfDiffGenotypes(vcf);
        numOfDiffAlleles(vcf);
        numOfEffAlleles(vcf);
        numOfHeterozygotes(vcf);
        numOfHomozygotes(vcf);
        calcObsHeterozygosity(vcf);
        calcExpHeterozygosity(vcf);
        calcUnbiasedExpHeterozygosity(vcf);
        calcShannonsI(vcf);
        calcSimpsonsI(vcf);
        calcFixationI(vcf);
        calcAR(vcf);

        //TODO:
        // Percentage of missing alleles
        // Percentage of missing genotypes
        // PIC (Polymorphism Information Content)
        // Sample (Individual) subsetting
        // Loci (SNP) subsetting
        // Command line parser with parameters
        // Lot other things..
    }

    public static Map readVCF(String filePath) {

        String[] columnNames = new String[0];
        Map<String, List<String>> locusData = new HashMap<>();

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath), StandardCharsets.UTF_8),1024*1024)) {
            String line = null;
            // Parse header
            while ((line = br.readLine()) != null) {
                if (line.charAt(0) == '#' && line.charAt(1) != '#') {
                    columnNames = line.split("\t");
                    break;
                }
            }
            // Create a map to store the column data
            for (int i = 0; i < columnNames.length; i++) {
                if (i < 9){
                    locusData.put(columnNames[i], new ArrayList<>());
                }
            }
            locusData.put("LOCUSSTAT", new ArrayList<String>());

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath), StandardCharsets.UTF_8),1024*1024)) {
            String line = null;
            long start = System.currentTimeMillis();
            // Read the data row by row and store the column data in the map
            while ((line = br.readLine()) != null) {
                // Genotype counts
                AtomicInteger homoCounter = new AtomicInteger(0);   // Counts homozygote allele combinations
                AtomicInteger hetCounter = new AtomicInteger(0);    // Counts heterozygote allele combinations
                AtomicInteger missCounter = new AtomicInteger(0);   // Counts missing allele combinations (.|. or ./.), missing genotypes
                Set<String> uniqueGenotypes = new HashSet<>();              // Counts unique (different) genotypes
                Set<String> uniqueAlleles = new HashSet<>();                // Counts unique (different) alleles

                // Allele counts
                AtomicInteger refCounter = new AtomicInteger(0);
                AtomicInteger altCounter = new AtomicInteger(0);

                if (line.charAt(0) != '#') {
                    String[] values = line.split("\t");
                    for (int i = 0; i < values.length; i++) {
                        if(i < 9){
                            String columnName = columnNames[i];
                            List<String> valuesList = locusData.get(columnName);
                            valuesList.add(values[i]);
                        } else if (i > 8) {
                            String genotype = values[i].split(":")[0].replace("|", "-").replace("/", "-");
                            // Genotype counts
                            if (genotype.equals("0-0") || genotype.equals("1-1")){
                                homoCounter.incrementAndGet();
                            } else if (genotype.equals("1-0") || genotype.equals("0-1")) {
                                hetCounter.incrementAndGet();
                            } else if (genotype.equals(".-.")) {
                                missCounter.incrementAndGet();
                            }
                            // Unique genotypes
                            if (genotype.equals("0-0") || genotype.equals("1-1") || genotype.equals("1-0") || genotype.equals("0-1")){
                                uniqueGenotypes.add(genotype);
                            }
                            // Allele counts
                            String[] alleles = genotype.split("-");
                            Arrays.stream(alleles).forEach(allele -> {
                                if(!allele.equals(".")){
                                    if (allele.equals("0")){
                                        refCounter.incrementAndGet();
                                        uniqueAlleles.add(allele);
                                    }
                                    if (allele.equals("1")){
                                        altCounter.incrementAndGet();
                                        uniqueAlleles.add(allele);
                                    }
                                }
                            });
                        }
                    }
                    String concatenatedValues = String.join(":",
                            String.valueOf(homoCounter.get()),  //1st column [0]
                            String.valueOf(hetCounter.get()),   //2nd column [1]
                            String.valueOf(missCounter.get()),  //3rd column [2]
                            String.valueOf(refCounter.get()),   //4th column [3]
                            String.valueOf(altCounter.get()),   //5th column [4]
                            String.valueOf(uniqueGenotypes.size()), //6th column [5]
                            String.valueOf(uniqueAlleles.size()));  //7th column [6]
                    locusData.get("LOCUSSTAT").add(concatenatedValues);
                }
            }
            // Print locus stats
            var locis = (List<String>) locusData.get("LOCUSSTAT");
            locis.forEach(System.out::println);

            //Print elapsed time
            long finish = System.currentTimeMillis();
            long timeElapsed = finish - start;
            long minutes = (timeElapsed / 1000) / 60;
            long seconds = (timeElapsed / 1000) % 60;
            System.out.println(timeElapsed + " Milliseconds = "+ minutes + " minutes and "+ seconds + " seconds.");

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return locusData;
    }

    public static void numOfMissingGT(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumMissing = new AtomicInteger(0);
        lstat.forEach(locus->{
            int missing = Integer.parseInt(locus.split(":")[2]);
            sumMissing.addAndGet(missing);
        });
        System.out.println("Number of missing genotypes: "+sumMissing.get());
    }

    public static void numOfREFAllele(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumRefAllele = new AtomicInteger(0);
        lstat.forEach(locus->{
            int refAllele = Integer.parseInt(locus.split(":")[3]);
            sumRefAllele.addAndGet(refAllele);
        });
        System.out.println("Number of REF alleles: "+sumRefAllele.get());
    }

    public static void numOfALTAllele(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumAltAllele = new AtomicInteger(0);
        lstat.forEach(locus->{
            int altAllele = Integer.parseInt(locus.split(":")[4]);
            sumAltAllele.addAndGet(altAllele);
        });
        System.out.println("Number of ALT alleles: "+sumAltAllele.get());
    }

    public static void numOfGenotypes(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumHomo = new AtomicInteger(0);
        lstat.forEach(locus->{
            int homozygote = Integer.parseInt(locus.split(":")[0]);
            sumHomo.addAndGet(homozygote);
        });
        AtomicInteger sumHet = new AtomicInteger(0);
        lstat.forEach(locus->{
            int heterozygote = Integer.parseInt(locus.split(":")[1]);
            sumHet.addAndGet(heterozygote);
        });
        int numOfGenotypes = sumHomo.get()+sumHet.get();
        System.out.println("Number of genotypes: "+numOfGenotypes);
    }

    public static void numOfDiffGenotypes(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumUniqueGenotype = new AtomicInteger(0);
        lstat.forEach(locus->{
            int uniqueGenotype = Integer.parseInt(locus.split(":")[5]);
            sumUniqueGenotype.addAndGet(uniqueGenotype);
        });
        double avgNumOfDiffGenotypes = (double) sumUniqueGenotype.get() / lstat.size();
        double rounded = Math.round(avgNumOfDiffGenotypes * 1000.0) / 1000.0;
        System.out.println("Average number of different genotypes (Ng): "+rounded);
    }

    public static void numOfDiffAlleles(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumDiffAlleles = new AtomicInteger(0);
        lstat.forEach(locus->{
            int diffAllele = Integer.parseInt(locus.split(":")[6]);
            sumDiffAlleles.addAndGet(diffAllele);
        });
        double avgNumOfDiffAlleles = (double) sumDiffAlleles.get() / lstat.size();
        double rounded = Math.round(avgNumOfDiffAlleles * 1000.0) / 1000.0;
        System.out.println("Average number of different alleles (Na): "+rounded);
    }

    public static void numOfEffAlleles(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicReference<Double> sumLociEffAlleles = new AtomicReference<>(0.0);
        lstat.forEach(locus->{
            int allAlleleAtLocus = Integer.parseInt(locus.split(":")[3]) + Integer.valueOf(locus.split(":")[4]);
            double refFreq = (double) Integer.parseInt(locus.split(":")[3]) / allAlleleAtLocus;
            double altFreq = (double) Integer.parseInt(locus.split(":")[4]) / allAlleleAtLocus;
            HashMap<String, Double> alleleFreq = new HashMap<>();
            alleleFreq.put("REF",refFreq);
            alleleFreq.put("ALT",altFreq);
            // h = 1 – sum(pi^2)
            // Ne = 1/(1-h)
            alleleFreq.forEach((key, value) -> {
                alleleFreq.put(key, Math.pow(value, 2));
            });
            double h = 1 - (alleleFreq.values().stream()
                    .mapToDouble(Double::doubleValue)
                    .sum());
            double Ne = 1/ (1 - h);
            sumLociEffAlleles.getAndUpdate(v -> v + Ne);
        });
        double avgnumOfEffAlleles = sumLociEffAlleles.get() / lstat.size();
        double rounded = Math.round(avgnumOfEffAlleles * 1000.0) / 1000.0;
        System.out.println("Average number of effective alleles (Ne): "+rounded);
    }

    public static void numOfHeterozygotes(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumHeterozygotes = new AtomicInteger(0);
        lstat.forEach(locus->{
            int heterozygotes = Integer.parseInt(locus.split(":")[1]);
            sumHeterozygotes.addAndGet(heterozygotes);
        });
        System.out.println("Number of heterozygotes: "+sumHeterozygotes.get());
    }

    public static void numOfHomozygotes(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumHomozygotes = new AtomicInteger(0);
        lstat.forEach(locus->{
            int homozygotes = Integer.parseInt(locus.split(":")[0]);
            sumHomozygotes.addAndGet(homozygotes);
        });
        System.out.println("Number of homozygotes: "+sumHomozygotes.get());
    }

    public static void calcObsHeterozygosity(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicReference<Double> sumObsHet = new AtomicReference<>(0.0);
        lstat.forEach(locus->{
            int numOfHomozygotes = Integer.parseInt(locus.split(":")[0]);
            int numOfHeterozygote = Integer.parseInt(locus.split(":")[1]);
            int numOfGenotypes = numOfHomozygotes + numOfHeterozygote;
            double obsHet = (double) numOfHeterozygote / numOfGenotypes;
            sumObsHet.updateAndGet(value -> value + obsHet);
        });
        double summedAvgObsHet = sumObsHet.get() / lstat.size();
        double rounded = Math.round(summedAvgObsHet * 1000.0) / 1000.0;
        System.out.println("Average Observed Heterozygosity (Ho): "+rounded);
    }

    public static void calcExpHeterozygosity(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicReference<Double> sumExpHet = new AtomicReference<>(0.0);
        lstat.forEach(locus->{
            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
            int numOfAlleles = numOfRefAllele + numOfAltAllele;
            // h = 1 – p^2 – q^2
            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
            double expHet = 1 - refAlleleFreq2 - altAlleleFreq2;
            sumExpHet.updateAndGet(value -> value + expHet);
        });
        double summedAvgExpHet = sumExpHet.get() / lstat.size();
        double rounded = Math.round(summedAvgExpHet * 1000.0) / 1000.0;
        System.out.println("Average Expected Heterozygosity (He): "+rounded);
    }

    public static void calcUnbiasedExpHeterozygosity(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicReference<Double> sumUnbiasedExpHet = new AtomicReference<>(0.0);
        lstat.forEach(locus->{
            int numOfHomGenotypes = Integer.parseInt(locus.split(":")[0]);
            int numOfHetGenotypes = Integer.parseInt(locus.split(":")[1]);
            int numOfGenotypes = numOfHomGenotypes + numOfHetGenotypes;
            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
            int numOfAlleles = numOfRefAllele + numOfAltAllele;
            // h = 1 – sum(pi^2)
            // uHe = Unbiased Expected Heterozygosity = (2N / (2N-1)) * h
            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
            double expHet = 1 - (refAlleleFreq2 + altAlleleFreq2);
            double twoN = 2 * numOfGenotypes;
            double uHe = (twoN / (twoN-1)) * expHet;
            sumUnbiasedExpHet.updateAndGet(value -> value + uHe);
        });
        double summedAvgUnbiasedExpHet = sumUnbiasedExpHet.get() / lstat.size();
        double rounded = Math.round(summedAvgUnbiasedExpHet * 1000.0) / 1000.0;
        System.out.println("Average Unbiased Expected Heterozygosity (uHe): "+rounded);
    }

    public static void calcShannonsI(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicReference<Double> sumShannon = new AtomicReference<>(0.0);
        lstat.forEach(locus->{
            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
            int numOfAlleles = numOfRefAllele + numOfAltAllele;
            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
            double LNrefAlleleFreq;
            if (refAlleleFreq != 0.0) {
                LNrefAlleleFreq = Math.log(refAlleleFreq);
            } else {
                LNrefAlleleFreq = 0.0;
            }
            double LNAlleleFreq;
            if (altAlleleFreq != 0.0) {
                LNAlleleFreq = Math.log(altAlleleFreq);
            } else {
                LNAlleleFreq = 0.0;
            }
            double refAllele = refAlleleFreq * LNrefAlleleFreq;
            double altAllele = altAlleleFreq * LNAlleleFreq;
            double ShannonI = refAllele + altAllele;
            sumShannon.updateAndGet(value -> value + ShannonI);
        });
        double summedAvgShannon = sumShannon.get() / lstat.size();
        double rounded = Math.round(summedAvgShannon * 1000.0) / 1000.0;
        System.out.println("Average Shannon's Information Index (SI): "+rounded);
    }

    public static void calcSimpsonsI(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicReference<Double> sumSimpson = new AtomicReference<>(0.0);
        lstat.forEach(locus->{
            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
            int numOfAlleles = numOfRefAllele + numOfAltAllele;
            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
            double sumAlleleFreqs = refAlleleFreq2 + altAlleleFreq2;
            double SimpsonD = 1 - sumAlleleFreqs;
            sumSimpson.updateAndGet(value -> value + SimpsonD);
        });
        double summedAvgSimpson = sumSimpson.get() / lstat.size();
        double rounded = Math.round(summedAvgSimpson * 1000.0) / 1000.0;
        System.out.println("Average Simpson's Diversity Index (D): "+rounded);
    }

    public static void calcFixationI(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicReference<Double> sumF = new AtomicReference<>(0.0);
        lstat.forEach(locus->{
            // Calculate He
            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
            int numOfAlleles = numOfRefAllele + numOfAltAllele;
            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
            double expHet = 1 - refAlleleFreq2 - altAlleleFreq2;
            // Calculate Ho
            int numOfHomozygotes = Integer.parseInt(locus.split(":")[0]);
            int numOfHeterozygote = Integer.parseInt(locus.split(":")[1]);
            int numOfGenotypes = numOfHomozygotes + numOfHeterozygote;
            double obsHet = (double) numOfHeterozygote / numOfGenotypes;
            // Calculate F: 1 - (Ho / He)
            double Findex;
            if(obsHet != 0 && expHet != 0){
                Findex = 1 - (obsHet / expHet);
            } else {
                Findex = 0.0;
            }
            sumF.updateAndGet(value -> value + Findex);
        });
        double summedFindex = sumF.get() / lstat.size();
        double rounded = Math.round(summedFindex * 1000.0) / 1000.0;
        System.out.println("Average Fixation Index (F): "+rounded);
    }

    public static void calcAR(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicReference<Double> sumAR = new AtomicReference<>(0.0);
        lstat.forEach(locus->{
            // Count 'n' the number of samples
            int numOfHomozygotes = Integer.parseInt(locus.split(":")[0]);
            int numOfHeterozygote = Integer.parseInt(locus.split(":")[1]);
            int numOfMissingGenotyes = Integer.parseInt(locus.split(":")[2]);
            int numSamples = numOfHomozygotes + numOfHeterozygote + numOfMissingGenotyes;
            // Count 'N' the number of "genes", number of alleles observed
            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
            int numAlleles = numOfRefAllele + numOfAltAllele;
            // Calculate based on Hurlbert (1971)
            int N = numAlleles - numOfRefAllele;
            int n = numSamples; // n must be between 0 and N
            double numCombinationsInNumeratorRef;
            if (n < 0 || n > N) {
                numCombinationsInNumeratorRef = 0.0;    // By convention (N choose n) = 0 when N < n
            } else {
                numCombinationsInNumeratorRef = (long) CombinatoricsUtils.binomialCoefficientDouble(N, n);
            }
            N = numAlleles - numOfAltAllele;
            n = numSamples; // n must be between 0 and N
            double numCombinationsInNumeratorAlt;
            if (n < 0 || n > N) {
                numCombinationsInNumeratorAlt = 0.0;    // By convention (N choose n) = 0 when N < n
            } else {
                numCombinationsInNumeratorAlt = (long) CombinatoricsUtils.binomialCoefficientDouble(N, n);
            }
            N = numAlleles;
            n = numSamples; // n must be between 0 and N
            double numCombinationsInDenominator;
            if (n < 0 || n > N) {
                numCombinationsInDenominator = 0.0;    // By convention (N choose n) = 0 when N < n
            } else {
                numCombinationsInDenominator = (long) CombinatoricsUtils.binomialCoefficientDouble(N, n);
            }
            double refAR;
            double altAR;
            if (numCombinationsInDenominator != 0.0){
                refAR = numCombinationsInNumeratorRef / numCombinationsInDenominator;
                altAR = numCombinationsInNumeratorAlt / numCombinationsInDenominator;
            } else {
                refAR = 0.0;
                altAR = 0.0;
            }
            double sumARForLocus = (1 - refAR) + (1 - altAR);
            sumAR.updateAndGet(value -> value + sumARForLocus);
        });
        double summedAR = sumAR.get() / lstat.size();
        double rounded = Math.round(summedAR * 1000.0) / 1000.0;
        System.out.println("Average Allelic Richness (Ar): "+rounded);
    }


}