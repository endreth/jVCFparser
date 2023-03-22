package controller;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class ParseVCF {

    public Map readVCF(String filePath) {

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

}
