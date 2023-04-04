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

        options.addOption(new Option("mg", "missg", false, "Missing genotype counts"));
        options.addOption(new Option("ra", "refa", false, "REF allele counts"));
        options.addOption(new Option("aa", "alta", false, "ALT allele counts"));
        options.addOption(new Option("gc", "gcounts", false, "Genotype counts"));
        options.addOption(new Option("dgc", "diffgcounts", false, "Different genotype counts"));
        options.addOption(new Option("da", "dacounts", false, "Different allele counts"));
        options.addOption(new Option("ea", "eacounts", false, "Effective allele counts"));
        options.addOption(new Option("het", "hetcounts", false, "Heterozygote counts"));
        options.addOption(new Option("hom", "homcounts", false, "Homozygote counts"));

        options.addOption(new Option("oh", "obshet", false, "Average Observed Heterozygosity (Ho)"));
        options.addOption(new Option("eh", "exphet", false, "Average Expected Heterozygosity (He)"));
        options.addOption(new Option("ueh", "uexphet", false, "Average Unbiased Expected Heterozygosity (uHe)"));
        options.addOption(new Option("sh", "shann", false, "Average Shannon's Information Index (H)"));
        options.addOption(new Option("sd", "simp", false, "Average Simpson's Diversity Index (D)"));
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

            if(cmd.hasOption("mg")) {
                clc.numOfMissingGT(vcf.getVcfMap());
            }
            if(cmd.hasOption("ra")) {
                clc.numOfREFAllele(vcf.getVcfMap());
            }
            if(cmd.hasOption("aa")) {
                clc.numOfALTAllele(vcf.getVcfMap());
            }
            if(cmd.hasOption("gc")) {
                clc.numOfGenotypes(vcf.getVcfMap());
            }
            if(cmd.hasOption("dgc")) {
                clc.numOfDiffGenotypes(vcf.getVcfMap());
            }
            if(cmd.hasOption("da")) {
                clc.numOfDiffAlleles(vcf.getVcfMap());
            }
            if(cmd.hasOption("ea")) {
                clc.numOfEffAlleles(vcf.getVcfMap());
            }
            if(cmd.hasOption("het")) {
                clc.numOfHeterozygotes(vcf.getVcfMap());
            }
            if(cmd.hasOption("hom")) {
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
            if(cmd.hasOption("sd")) {
                clc.calcSimpsonsD(vcf.getVcfMap());
            }
            if(cmd.hasOption("fx")) {
                clc.calcFixationI(vcf.getVcfMap());
            }
            if(cmd.hasOption("ar")) {
                clc.calcAR(vcf.getVcfMap());
            }

        } catch (ParseException e) {
            System.out.println(e.getMessage());
            helper.printHelp("Usage:", options);
            System.exit(0);
        }

    }


}