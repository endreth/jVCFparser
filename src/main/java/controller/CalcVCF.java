package controller;

import org.apache.commons.math3.util.CombinatoricsUtils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

public class CalcVCF {

    public void numOfMissingGT(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumMissing = new AtomicInteger(0);
        lstat.forEach(locus->{
            int missing = Integer.parseInt(locus.split(":")[2]);
            sumMissing.addAndGet(missing);
        });
        System.out.println("Number of missing genotypes: "+sumMissing.get());
    }

    public void numOfREFAllele(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumRefAllele = new AtomicInteger(0);
        lstat.forEach(locus->{
            int refAllele = Integer.parseInt(locus.split(":")[3]);
            sumRefAllele.addAndGet(refAllele);
        });
        System.out.println("Number of REF alleles: "+sumRefAllele.get());
    }

    public void numOfALTAllele(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumAltAllele = new AtomicInteger(0);
        lstat.forEach(locus->{
            int altAllele = Integer.parseInt(locus.split(":")[4]);
            sumAltAllele.addAndGet(altAllele);
        });
        System.out.println("Number of ALT alleles: "+sumAltAllele.get());
    }

    public void numOfGenotypes(Map vcf){
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

    public void numOfDiffGenotypes(Map vcf){
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

    public void numOfDiffAlleles(Map vcf){
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

    public void numOfEffAlleles(Map vcf){
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

    public void numOfHeterozygotes(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumHeterozygotes = new AtomicInteger(0);
        lstat.forEach(locus->{
            int heterozygotes = Integer.parseInt(locus.split(":")[1]);
            sumHeterozygotes.addAndGet(heterozygotes);
        });
        System.out.println("Number of heterozygotes: "+sumHeterozygotes.get());
    }

    public void numOfHomozygotes(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        AtomicInteger sumHomozygotes = new AtomicInteger(0);
        lstat.forEach(locus->{
            int homozygotes = Integer.parseInt(locus.split(":")[0]);
            sumHomozygotes.addAndGet(homozygotes);
        });
        System.out.println("Number of homozygotes: "+sumHomozygotes.get());
    }

    public void calcObsHeterozygosity(Map vcf){
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

    public void calcExpHeterozygosity(Map vcf){
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

    public void calcUnbiasedExpHeterozygosity(Map vcf){
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

    public void calcShannonsI(Map vcf){
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
        System.out.println("Average Shannon's Information Index (H): "+rounded);
    }

    public void calcSimpsonsI(Map vcf){
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

    public void calcFixationI(Map vcf){
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

    public void calcAR(Map vcf){
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
