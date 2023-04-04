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
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nNumber of missing genotypes by locus:\n<locus> <no. missing genotypes>");
        AtomicInteger sumMissing = new AtomicInteger(0);
        for (int i = 0; i < lstat.size(); i++) {
            int missing = Integer.parseInt(lstat.get(i).split(":")[2]);
            sumMissing.addAndGet(missing);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            System.out.println(locus+"\t"+missing);
        }
        System.out.println("---");
        System.out.println("Total number of missing genotypes: "+sumMissing.get());
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicInteger sumMissing = new AtomicInteger(0);
//        lstat.forEach(locus->{
//            int missing = Integer.parseInt(locus.split(":")[2]);
//            sumMissing.addAndGet(missing);
//        });
//        System.out.println("Number of missing genotypes: "+sumMissing.get());
    }

    public void numOfREFAllele(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nNumber of REF alleles by locus:\n<locus> <no. REF alleles>");
        AtomicInteger sumRefAllele = new AtomicInteger(0);
        for (int i = 0; i < lstat.size(); i++) {
            int refAllele = Integer.parseInt(lstat.get(i).split(":")[3]);
            sumRefAllele.addAndGet(refAllele);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            System.out.println(locus+"\t"+refAllele);
        }
        System.out.println("---");
        System.out.println("Total number of REF alleles: "+sumRefAllele.get());
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicInteger sumRefAllele = new AtomicInteger(0);
//        lstat.forEach(locus->{
//            int refAllele = Integer.parseInt(locus.split(":")[3]);
//            sumRefAllele.addAndGet(refAllele);
//        });
//        System.out.println("Number of REF alleles: "+sumRefAllele.get());
    }

    public void numOfALTAllele(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nNumber of ALT alleles by locus:\n<locus> <no. ALT alleles>");
        AtomicInteger sumAltAllele = new AtomicInteger(0);
        for (int i = 0; i < lstat.size(); i++) {
            int altAllele = Integer.parseInt(lstat.get(i).split(":")[4]);
            sumAltAllele.addAndGet(altAllele);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            System.out.println(locus+"\t"+altAllele);
        }
        System.out.println("---");
        System.out.println("Total number of ALT alleles: "+sumAltAllele.get());
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicInteger sumAltAllele = new AtomicInteger(0);
//        lstat.forEach(locus->{
//            int altAllele = Integer.parseInt(locus.split(":")[4]);
//            sumAltAllele.addAndGet(altAllele);
//        });
//        System.out.println("Number of ALT alleles: "+sumAltAllele.get());
    }

    public void numOfGenotypes(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nNumber of genotypes by locus:\n<locus> <homozygotes> <heterozygotes>");
        AtomicInteger sumHomo = new AtomicInteger(0);
        AtomicInteger sumHet = new AtomicInteger(0);
        for (int i = 0; i < lstat.size(); i++) {
            int homozygote = Integer.parseInt(lstat.get(i).split(":")[0]);
            sumHomo.addAndGet(homozygote);
            int heterozygote = Integer.parseInt(lstat.get(i).split(":")[1]);
            sumHet.addAndGet(heterozygote);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            System.out.println(locus+"\t"+homozygote+"\t"+heterozygote);
        }
        int numOfGenotypes = sumHomo.get()+sumHet.get();
        System.out.println("---");
        System.out.println("Total number of genotypes: "+numOfGenotypes);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicInteger sumHomo = new AtomicInteger(0);
//        lstat.forEach(locus->{
//            int homozygote = Integer.parseInt(locus.split(":")[0]);
//            sumHomo.addAndGet(homozygote);
//        });
//        AtomicInteger sumHet = new AtomicInteger(0);
//        lstat.forEach(locus->{
//            int heterozygote = Integer.parseInt(locus.split(":")[1]);
//            sumHet.addAndGet(heterozygote);
//        });
//        int numOfGenotypes = sumHomo.get()+sumHet.get();
//        System.out.println("Number of genotypes: "+numOfGenotypes);
    }

    public void numOfDiffGenotypes(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nNumber of different genotypes (Ng) by locus:\n<locus> <no. diff. genotypes>");
        AtomicInteger sumUniqueGenotype = new AtomicInteger(0);
        for (int i = 0; i < lstat.size(); i++) {
            int uniqueGenotype = Integer.parseInt(lstat.get(i).split(":")[5]);
            sumUniqueGenotype.addAndGet(uniqueGenotype);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            System.out.println(locus+"\t"+uniqueGenotype);
        }
        double avgNumOfDiffGenotypes = (double) sumUniqueGenotype.get() / lstat.size();
        double rounded = Math.round(avgNumOfDiffGenotypes * 1000.0) / 1000.0;
        System.out.println("---");
        System.out.println("Average number of different genotypes (Ng): "+rounded);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicInteger sumUniqueGenotype = new AtomicInteger(0);
//        lstat.forEach(locus->{
//            int uniqueGenotype = Integer.parseInt(locus.split(":")[5]);
//            sumUniqueGenotype.addAndGet(uniqueGenotype);
//        });
//        double avgNumOfDiffGenotypes = (double) sumUniqueGenotype.get() / lstat.size();
//        double rounded = Math.round(avgNumOfDiffGenotypes * 1000.0) / 1000.0;
//        System.out.println("Average number of different genotypes (Ng): "+rounded);
    }

    public void numOfDiffAlleles(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nNumber of different alleles (Na) by locus:\n<locus> <no. diff. alleles>");
        AtomicInteger sumDiffAlleles = new AtomicInteger(0);
        for (int i = 0; i < lstat.size(); i++) {
            int diffAllele = Integer.parseInt(lstat.get(i).split(":")[6]);
            sumDiffAlleles.addAndGet(diffAllele);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            System.out.println(locus+"\t"+diffAllele);
        }
        double avgNumOfDiffAlleles = (double) sumDiffAlleles.get() / lstat.size();
        double rounded = Math.round(avgNumOfDiffAlleles * 1000.0) / 1000.0;
        System.out.println("---");
        System.out.println("Average number of different alleles (Na): "+rounded);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicInteger sumDiffAlleles = new AtomicInteger(0);
//        lstat.forEach(locus->{
//            int diffAllele = Integer.parseInt(locus.split(":")[6]);
//            sumDiffAlleles.addAndGet(diffAllele);
//        });
//        double avgNumOfDiffAlleles = (double) sumDiffAlleles.get() / lstat.size();
//        double rounded = Math.round(avgNumOfDiffAlleles * 1000.0) / 1000.0;
//        System.out.println("Average number of different alleles (Na): "+rounded);
    }

    public void numOfEffAlleles(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nNumber of effective alleles (Ne) by locus:\n<locus> <no. eff. alleles>");
        AtomicReference<Double> sumLociEffAlleles = new AtomicReference<>(0.0);
        for (int i = 0; i < lstat.size(); i++) {
            int allAlleleAtLocus = Integer.parseInt(lstat.get(i).split(":")[3]) + Integer.valueOf(lstat.get(i).split(":")[4]);
            double refFreq = (double) Integer.parseInt(lstat.get(i).split(":")[3]) / allAlleleAtLocus;
            double altFreq = (double) Integer.parseInt(lstat.get(i).split(":")[4]) / allAlleleAtLocus;
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
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            double roundedNe = Math.round(Ne * 1000.0) / 1000.0;
            System.out.println(locus+"\t"+roundedNe);
        }
        double avgnumOfEffAlleles = sumLociEffAlleles.get() / lstat.size();
        double rounded = Math.round(avgnumOfEffAlleles * 1000.0) / 1000.0;
        System.out.println("---");
        System.out.println("Average number of effective alleles (Ne): "+rounded);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicReference<Double> sumLociEffAlleles = new AtomicReference<>(0.0);
//        lstat.forEach(locus->{
//            int allAlleleAtLocus = Integer.parseInt(locus.split(":")[3]) + Integer.valueOf(locus.split(":")[4]);
//            double refFreq = (double) Integer.parseInt(locus.split(":")[3]) / allAlleleAtLocus;
//            double altFreq = (double) Integer.parseInt(locus.split(":")[4]) / allAlleleAtLocus;
//            HashMap<String, Double> alleleFreq = new HashMap<>();
//            alleleFreq.put("REF",refFreq);
//            alleleFreq.put("ALT",altFreq);
//            // h = 1 – sum(pi^2)
//            // Ne = 1/(1-h)
//            alleleFreq.forEach((key, value) -> {
//                alleleFreq.put(key, Math.pow(value, 2));
//            });
//            double h = 1 - (alleleFreq.values().stream()
//                    .mapToDouble(Double::doubleValue)
//                    .sum());
//            double Ne = 1/ (1 - h);
//            sumLociEffAlleles.getAndUpdate(v -> v + Ne);
//        });
//        double avgnumOfEffAlleles = sumLociEffAlleles.get() / lstat.size();
//        double rounded = Math.round(avgnumOfEffAlleles * 1000.0) / 1000.0;
//        System.out.println("Average number of effective alleles (Ne): "+rounded);
    }

    public void numOfHeterozygotes(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nNumber of heterozygotes by locus:\n<locus> <no. of heterozygotes>");
        AtomicInteger sumHeterozygotes = new AtomicInteger(0);
        for (int i = 0; i < lstat.size(); i++) {
            int heterozygotes = Integer.parseInt(lstat.get(i).split(":")[1]);
            sumHeterozygotes.addAndGet(heterozygotes);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            System.out.println(locus+"\t"+heterozygotes);
        }
        System.out.println("---");
        System.out.println("Total number of heterozygotes: "+sumHeterozygotes.get());
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicInteger sumHeterozygotes = new AtomicInteger(0);
//        lstat.forEach(locus->{
//            int heterozygotes = Integer.parseInt(locus.split(":")[1]);
//            sumHeterozygotes.addAndGet(heterozygotes);
//        });
//        System.out.println("Number of heterozygotes: "+sumHeterozygotes.get());
    }

    public void numOfHomozygotes(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nNumber of homozygotes by locus:\n<locus> <no. of homozygotes>");
        AtomicInteger sumHomozygotes = new AtomicInteger(0);
        for (int i = 0; i < lstat.size(); i++) {
            int homozygotes = Integer.parseInt(lstat.get(i).split(":")[0]);
            sumHomozygotes.addAndGet(homozygotes);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            System.out.println(locus+"\t"+homozygotes);
        }
        System.out.println("---");
        System.out.println("Total number of homozygotes: "+sumHomozygotes.get());
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicInteger sumHomozygotes = new AtomicInteger(0);
//        lstat.forEach(locus->{
//            int homozygotes = Integer.parseInt(locus.split(":")[0]);
//            sumHomozygotes.addAndGet(homozygotes);
//        });
//        System.out.println("Number of homozygotes: "+sumHomozygotes.get());
    }

    public void calcObsHeterozygosity(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nObserved Heterozygosity (Ho) by locus:\n<locus> <Ho>");
        AtomicReference<Double> sumObsHet = new AtomicReference<>(0.0);
        for (int i = 0; i < lstat.size(); i++) {
            int numOfHomozygotes = Integer.parseInt(lstat.get(i).split(":")[0]);
            int numOfHeterozygote = Integer.parseInt(lstat.get(i).split(":")[1]);
            int numOfGenotypes = numOfHomozygotes + numOfHeterozygote;
            double obsHet = (double) numOfHeterozygote / numOfGenotypes;
            sumObsHet.updateAndGet(value -> value + obsHet);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            double roundedHo = Math.round(obsHet * 1000.0) / 1000.0;
            System.out.println(locus+"\t"+roundedHo);
        }
        double summedAvgObsHet = sumObsHet.get() / lstat.size();
        double rounded = Math.round(summedAvgObsHet * 1000.0) / 1000.0;
        System.out.println("---");
        System.out.println("Average Observed Heterozygosity (Ho): "+rounded);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicReference<Double> sumObsHet = new AtomicReference<>(0.0);
//        lstat.forEach(locus->{
//            int numOfHomozygotes = Integer.parseInt(locus.split(":")[0]);
//            int numOfHeterozygote = Integer.parseInt(locus.split(":")[1]);
//            int numOfGenotypes = numOfHomozygotes + numOfHeterozygote;
//            double obsHet = (double) numOfHeterozygote / numOfGenotypes;
//            sumObsHet.updateAndGet(value -> value + obsHet);
//        });
//        double summedAvgObsHet = sumObsHet.get() / lstat.size();
//        double rounded = Math.round(summedAvgObsHet * 1000.0) / 1000.0;
//        System.out.println("Average Observed Heterozygosity (Ho): "+rounded);
    }

    public void calcExpHeterozygosity(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nExpected Heterozygosity (He) by locus:\n<locus> <He>");
        AtomicReference<Double> sumExpHet = new AtomicReference<>(0.0);
        for (int i = 0; i < lstat.size(); i++) {
            int numOfRefAllele = Integer.parseInt(lstat.get(i).split(":")[3]);
            int numOfAltAllele = Integer.parseInt(lstat.get(i).split(":")[4]);
            int numOfAlleles = numOfRefAllele + numOfAltAllele;
            // h = 1 – p^2 – q^2
            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
            double expHet = 1 - refAlleleFreq2 - altAlleleFreq2;
            sumExpHet.updateAndGet(value -> value + expHet);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            double roundedHe = Math.round(expHet * 1000.0) / 1000.0;
            System.out.println(locus+"\t"+roundedHe);
        }
        double summedAvgExpHet = sumExpHet.get() / lstat.size();
        double rounded = Math.round(summedAvgExpHet * 1000.0) / 1000.0;
        System.out.println("---");
        System.out.println("Average Expected Heterozygosity (He): "+rounded);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicReference<Double> sumExpHet = new AtomicReference<>(0.0);
//        lstat.forEach(locus->{
//            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
//            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
//            int numOfAlleles = numOfRefAllele + numOfAltAllele;
//            // h = 1 – p^2 – q^2
//            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
//            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
//            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
//            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
//            double expHet = 1 - refAlleleFreq2 - altAlleleFreq2;
//            sumExpHet.updateAndGet(value -> value + expHet);
//        });
//        double summedAvgExpHet = sumExpHet.get() / lstat.size();
//        double rounded = Math.round(summedAvgExpHet * 1000.0) / 1000.0;
//        System.out.println("Average Expected Heterozygosity (He): "+rounded);
    }

    public void calcUnbiasedExpHeterozygosity(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nUnbiased Expected Heterozygosity (uHe) by locus:\n<locus> <uHe>");
        AtomicReference<Double> sumUnbiasedExpHet = new AtomicReference<>(0.0);
        for (int i = 0; i < lstat.size(); i++) {
            int numOfHomGenotypes = Integer.parseInt(lstat.get(i).split(":")[0]);
            int numOfHetGenotypes = Integer.parseInt(lstat.get(i).split(":")[1]);
            int numOfGenotypes = numOfHomGenotypes + numOfHetGenotypes;
            int numOfRefAllele = Integer.parseInt(lstat.get(i).split(":")[3]);
            int numOfAltAllele = Integer.parseInt(lstat.get(i).split(":")[4]);
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
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            double roundeduHe = Math.round(uHe * 1000.0) / 1000.0;
            System.out.println(locus+"\t"+roundeduHe);
        }
        double summedAvgUnbiasedExpHet = sumUnbiasedExpHet.get() / lstat.size();
        double rounded = Math.round(summedAvgUnbiasedExpHet * 1000.0) / 1000.0;
        System.out.println("---");
        System.out.println("Average Unbiased Expected Heterozygosity (uHe): "+rounded);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicReference<Double> sumUnbiasedExpHet = new AtomicReference<>(0.0);
//        lstat.forEach(locus->{
//            int numOfHomGenotypes = Integer.parseInt(locus.split(":")[0]);
//            int numOfHetGenotypes = Integer.parseInt(locus.split(":")[1]);
//            int numOfGenotypes = numOfHomGenotypes + numOfHetGenotypes;
//            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
//            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
//            int numOfAlleles = numOfRefAllele + numOfAltAllele;
//            // h = 1 – sum(pi^2)
//            // uHe = Unbiased Expected Heterozygosity = (2N / (2N-1)) * h
//            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
//            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
//            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
//            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
//            double expHet = 1 - (refAlleleFreq2 + altAlleleFreq2);
//            double twoN = 2 * numOfGenotypes;
//            double uHe = (twoN / (twoN-1)) * expHet;
//            sumUnbiasedExpHet.updateAndGet(value -> value + uHe);
//        });
//        double summedAvgUnbiasedExpHet = sumUnbiasedExpHet.get() / lstat.size();
//        double rounded = Math.round(summedAvgUnbiasedExpHet * 1000.0) / 1000.0;
//        System.out.println("Average Unbiased Expected Heterozygosity (uHe): "+rounded);
    }

    public void calcShannonsI(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nShannon's Information Index (H) by locus:\n<locus> <H>");
        AtomicReference<Double> sumShannon = new AtomicReference<>(0.0);
        for (int i = 0; i < lstat.size(); i++) {
            int numOfRefAllele = Integer.parseInt(lstat.get(i).split(":")[3]);
            int numOfAltAllele = Integer.parseInt(lstat.get(i).split(":")[4]);
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
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            double roundedShannonI = Math.round(ShannonI * 1000.0) / 1000.0;
            System.out.println(locus+"\t"+roundedShannonI);
        }
        double summedAvgShannon = sumShannon.get() / lstat.size();
        double rounded = Math.round(summedAvgShannon * 1000.0) / 1000.0;
        System.out.println("---");
        System.out.println("Average Shannon's Information Index (H): "+rounded);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicReference<Double> sumShannon = new AtomicReference<>(0.0);
//        lstat.forEach(locus->{
//            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
//            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
//            int numOfAlleles = numOfRefAllele + numOfAltAllele;
//            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
//            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
//            double LNrefAlleleFreq;
//            if (refAlleleFreq != 0.0) {
//                LNrefAlleleFreq = Math.log(refAlleleFreq);
//            } else {
//                LNrefAlleleFreq = 0.0;
//            }
//            double LNAlleleFreq;
//            if (altAlleleFreq != 0.0) {
//                LNAlleleFreq = Math.log(altAlleleFreq);
//            } else {
//                LNAlleleFreq = 0.0;
//            }
//            double refAllele = refAlleleFreq * LNrefAlleleFreq;
//            double altAllele = altAlleleFreq * LNAlleleFreq;
//            double ShannonI = refAllele + altAllele;
//            sumShannon.updateAndGet(value -> value + ShannonI);
//        });
//        double summedAvgShannon = sumShannon.get() / lstat.size();
//        double rounded = Math.round(summedAvgShannon * 1000.0) / 1000.0;
//        System.out.println("Average Shannon's Information Index (H): "+rounded);
    }

    public void calcSimpsonsD(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nSimpson's Diversity Index (D) by locus:\n<locus> <D>");
        AtomicReference<Double> sumSimpson = new AtomicReference<>(0.0);
        for (int i = 0; i < lstat.size(); i++) {
            int numOfRefAllele = Integer.parseInt(lstat.get(i).split(":")[3]);
            int numOfAltAllele = Integer.parseInt(lstat.get(i).split(":")[4]);
            int numOfAlleles = numOfRefAllele + numOfAltAllele;
            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
            double sumAlleleFreqs = refAlleleFreq2 + altAlleleFreq2;
            double SimpsonD = 1 - sumAlleleFreqs;
            sumSimpson.updateAndGet(value -> value + SimpsonD);
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            double roundedSimpsonD = Math.round(SimpsonD * 1000.0) / 1000.0;
            System.out.println(locus+"\t"+roundedSimpsonD);
        }
        double summedAvgSimpson = sumSimpson.get() / lstat.size();
        double rounded = Math.round(summedAvgSimpson * 1000.0) / 1000.0;
        System.out.println("---");
        System.out.println("Average Simpson's Diversity Index (D): "+rounded);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicReference<Double> sumSimpson = new AtomicReference<>(0.0);
//        lstat.forEach(locus->{
//            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
//            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
//            int numOfAlleles = numOfRefAllele + numOfAltAllele;
//            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
//            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
//            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
//            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
//            double sumAlleleFreqs = refAlleleFreq2 + altAlleleFreq2;
//            double SimpsonD = 1 - sumAlleleFreqs;
//            sumSimpson.updateAndGet(value -> value + SimpsonD);
//        });
//        double summedAvgSimpson = sumSimpson.get() / lstat.size();
//        double rounded = Math.round(summedAvgSimpson * 1000.0) / 1000.0;
//        System.out.println("Average Simpson's Diversity Index (D): "+rounded);
    }

    public void calcFixationI(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nFixation Index (F) by locus:\n<locus> <F>");
        AtomicReference<Double> sumF = new AtomicReference<>(0.0);
        for (int i = 0; i < lstat.size(); i++) {
            // Calculate He
            int numOfRefAllele = Integer.parseInt(lstat.get(i).split(":")[3]);
            int numOfAltAllele = Integer.parseInt(lstat.get(i).split(":")[4]);
            int numOfAlleles = numOfRefAllele + numOfAltAllele;
            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
            double expHet = 1 - refAlleleFreq2 - altAlleleFreq2;
            // Calculate Ho
            int numOfHomozygotes = Integer.parseInt(lstat.get(i).split(":")[0]);
            int numOfHeterozygote = Integer.parseInt(lstat.get(i).split(":")[1]);
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
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            double roundedFindex = Math.round(Findex * 1000.0) / 1000.0;
            System.out.println(locus+"\t"+roundedFindex);
        }
        double summedFindex = sumF.get() / lstat.size();
        double rounded = Math.round(summedFindex * 1000.0) / 1000.0;
        System.out.println("---");
        System.out.println("Average Fixation Index (F): "+rounded);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicReference<Double> sumF = new AtomicReference<>(0.0);
//        lstat.forEach(locus->{
//            // Calculate He
//            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
//            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
//            int numOfAlleles = numOfRefAllele + numOfAltAllele;
//            double refAlleleFreq = (double) numOfRefAllele / numOfAlleles;
//            double altAlleleFreq = (double) numOfAltAllele / numOfAlleles;
//            double refAlleleFreq2 = Math.pow(refAlleleFreq, 2);
//            double altAlleleFreq2 = Math.pow(altAlleleFreq, 2);
//            double expHet = 1 - refAlleleFreq2 - altAlleleFreq2;
//            // Calculate Ho
//            int numOfHomozygotes = Integer.parseInt(locus.split(":")[0]);
//            int numOfHeterozygote = Integer.parseInt(locus.split(":")[1]);
//            int numOfGenotypes = numOfHomozygotes + numOfHeterozygote;
//            double obsHet = (double) numOfHeterozygote / numOfGenotypes;
//            // Calculate F: 1 - (Ho / He)
//            double Findex;
//            if(obsHet != 0 && expHet != 0){
//                Findex = 1 - (obsHet / expHet);
//            } else {
//                Findex = 0.0;
//            }
//            sumF.updateAndGet(value -> value + Findex);
//        });
//        double summedFindex = sumF.get() / lstat.size();
//        double rounded = Math.round(summedFindex * 1000.0) / 1000.0;
//        System.out.println("Average Fixation Index (F): "+rounded);
    }

    public void calcAR(Map vcf){
        var lstat = (List<String>) vcf.get("LOCUSSTAT");
        var lid = (List<String>) vcf.get("ID");
        System.out.println("\nAllelic Richness (Ar) by locus:\n<locus> <Ar>");
        AtomicReference<Double> sumAR = new AtomicReference<>(0.0);
        for (int i = 0; i < lstat.size(); i++) {
            // Count 'n' the number of samples
            int numOfHomozygotes = Integer.parseInt(lstat.get(i).split(":")[0]);
            int numOfHeterozygote = Integer.parseInt(lstat.get(i).split(":")[1]);
            int numOfMissingGenotyes = Integer.parseInt(lstat.get(i).split(":")[2]);
            int numSamples = numOfHomozygotes + numOfHeterozygote + numOfMissingGenotyes;
            // Count 'N' the number of "genes", number of alleles observed
            int numOfRefAllele = Integer.parseInt(lstat.get(i).split(":")[3]);
            int numOfAltAllele = Integer.parseInt(lstat.get(i).split(":")[4]);
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
            String locus = lid.get(i).split(":")[0] + "_"+lid.get(i).split(":")[1];
            double roundedsumARForLocus = Math.round(sumARForLocus * 1000.0) / 1000.0;
            System.out.println(locus+"\t"+roundedsumARForLocus);
        }
        double summedAR = sumAR.get() / lstat.size();
        double rounded = Math.round(summedAR * 1000.0) / 1000.0;
        System.out.println("---");
        System.out.println("Average Allelic Richness (Ar): "+rounded);
//        var lstat = (List<String>) vcf.get("LOCUSSTAT");
//        AtomicReference<Double> sumAR = new AtomicReference<>(0.0);
//        lstat.forEach(locus->{
//            // Count 'n' the number of samples
//            int numOfHomozygotes = Integer.parseInt(locus.split(":")[0]);
//            int numOfHeterozygote = Integer.parseInt(locus.split(":")[1]);
//            int numOfMissingGenotyes = Integer.parseInt(locus.split(":")[2]);
//            int numSamples = numOfHomozygotes + numOfHeterozygote + numOfMissingGenotyes;
//            // Count 'N' the number of "genes", number of alleles observed
//            int numOfRefAllele = Integer.parseInt(locus.split(":")[3]);
//            int numOfAltAllele = Integer.parseInt(locus.split(":")[4]);
//            int numAlleles = numOfRefAllele + numOfAltAllele;
//            // Calculate based on Hurlbert (1971)
//            int N = numAlleles - numOfRefAllele;
//            int n = numSamples; // n must be between 0 and N
//            double numCombinationsInNumeratorRef;
//            if (n < 0 || n > N) {
//                numCombinationsInNumeratorRef = 0.0;    // By convention (N choose n) = 0 when N < n
//            } else {
//                numCombinationsInNumeratorRef = (long) CombinatoricsUtils.binomialCoefficientDouble(N, n);
//            }
//            N = numAlleles - numOfAltAllele;
//            n = numSamples; // n must be between 0 and N
//            double numCombinationsInNumeratorAlt;
//            if (n < 0 || n > N) {
//                numCombinationsInNumeratorAlt = 0.0;    // By convention (N choose n) = 0 when N < n
//            } else {
//                numCombinationsInNumeratorAlt = (long) CombinatoricsUtils.binomialCoefficientDouble(N, n);
//            }
//            N = numAlleles;
//            n = numSamples; // n must be between 0 and N
//            double numCombinationsInDenominator;
//            if (n < 0 || n > N) {
//                numCombinationsInDenominator = 0.0;    // By convention (N choose n) = 0 when N < n
//            } else {
//                numCombinationsInDenominator = (long) CombinatoricsUtils.binomialCoefficientDouble(N, n);
//            }
//            double refAR;
//            double altAR;
//            if (numCombinationsInDenominator != 0.0){
//                refAR = numCombinationsInNumeratorRef / numCombinationsInDenominator;
//                altAR = numCombinationsInNumeratorAlt / numCombinationsInDenominator;
//            } else {
//                refAR = 0.0;
//                altAR = 0.0;
//            }
//            double sumARForLocus = (1 - refAR) + (1 - altAR);
//            sumAR.updateAndGet(value -> value + sumARForLocus);
//        });
//        double summedAR = sumAR.get() / lstat.size();
//        double rounded = Math.round(summedAR * 1000.0) / 1000.0;
//        System.out.println("Average Allelic Richness (Ar): "+rounded);
    }


}
