package controller;

import model.VCF;
import org.apache.commons.cli.*;

public class Main {
    public static void main(String[] args) {

        //TODO:
        // Percentage of missing alleles
        // Percentage of missing genotypes
        // PIC (Polymorphism Information Content)
        // Sample (Individual) subsetting
        // Loci (SNP) subsetting
        // Command line parser with parameters
        // Lot other things..

        Options options = new Options();

        Option source = new Option("f", "file", true, "VCF file path");
        source.setRequired(true);
        options.addOption(source);

        options.addOption(new Option("agc", "agcounts", false, "Allele-Genotype counts"));
        options.addOption(new Option("oh", "obshet", false, "Average Observed Heterozygosity (Ho)"));
        options.addOption(new Option("eh", "exphet", false, "Average Expected Heterozygosity (He)"));
        options.addOption(new Option("ueh", "uexphet", false, "Average Unbiased Expected Heterozygosity (uHe)"));
        options.addOption(new Option("sh", "shann", false, "Average Shannon's Information Index (H)"));
        options.addOption(new Option("si", "simp", false, "Average Simpson's Diversity Index (D)"));
        options.addOption(new Option("fx", "fix", false, "Average Fixation Index (F)"));
        options.addOption(new Option("ar", "arich", false, "Average Allelic Richness (Ar)"));

        CommandLine cmd;
        CommandLineParser parser = new BasicParser();
        HelpFormatter helper = new HelpFormatter();

        try {
            cmd = parser.parse(options, args);
            CalcVCF clc = new CalcVCF();
            ParseVCF parse = new ParseVCF();
            VCF vcf = new VCF();

            if (cmd.hasOption("f")) {
                String opt_config = cmd.getOptionValue("file");
                System.out.println("VCF file path set to " + opt_config);
                vcf.setVcfMap(parse.readVCF(opt_config));
            }

            if(cmd.hasOption("agc")) {
                clc.numOfMissingGT(vcf.getVcfMap());
                clc.numOfREFAllele(vcf.getVcfMap());
                clc.numOfALTAllele(vcf.getVcfMap());
                clc.numOfGenotypes(vcf.getVcfMap());
                clc.numOfDiffGenotypes(vcf.getVcfMap());
                clc.numOfDiffAlleles(vcf.getVcfMap());
                clc.numOfEffAlleles(vcf.getVcfMap());
                clc.numOfHeterozygotes(vcf.getVcfMap());
                clc.numOfHomozygotes(vcf.getVcfMap());
            }

            if(cmd.hasOption("oh")) {
                clc.calcObsHeterozygosity(vcf.getVcfMap());
            }
            if(cmd.hasOption("eh")) {
                clc.calcExpHeterozygosity(vcf.getVcfMap());
            }
            if(cmd.hasOption("ueh")) {
                clc.calcUnbiasedExpHeterozygosity(vcf.getVcfMap());
            }
            if(cmd.hasOption("sh")) {
                clc.calcShannonsI(vcf.getVcfMap());
            }
            if(cmd.hasOption("si")) {
                clc.calcSimpsonsI(vcf.getVcfMap());
            }
            if(cmd.hasOption("fx")) {
                clc.calcFixationI(vcf.getVcfMap());
            }
            if(cmd.hasOption("a")) {
                clc.calcAR(vcf.getVcfMap());
            }

        } catch (ParseException e) {
            System.out.println(e.getMessage());
            helper.printHelp("Usage:", options);
            System.exit(0);
        }

    }


}