import org.apache.commons.cli.*;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

import java.io.*;
import java.util.*;

public class Runner {

    public static void main(String[] args) throws IOException {
        long startingtime = System.currentTimeMillis();

        HashMap<Double, Double> hg_adjusted = new HashMap<>();
        HashMap<Double, Double> fej_adjusted = new HashMap<>();
        HashMap<Double, Double> ks_adjusted = new HashMap<>();

        //--------------------------------------------------------------------------------------------------------------
        //CLP
        Options options = new Options();
        //obo
        Option obo = new Option("obo", true, "obo file");
        options.addOption(obo);
        //root
        Option root = new Option("root", true, "root");
        options.addOption(root);
        //mapping
        Option mapping = new Option("mapping", true, "mapping");
        options.addOption(mapping);
        //mappinngtype
        Option mappingtype = new Option("mappingtype", true, "mappingtype");
        options.addOption(mappingtype);
        //enrich
        Option enrich = new Option("enrich", true, "enrich");
        options.addOption(enrich);
        //minsize
        Option minsize = new Option("minsize", true, "minsize");
        options.addOption(minsize);
        //maxsize
        Option maxsize = new Option("maxsize", true, "maxsize");
        options.addOption(maxsize);
        //outputoverlap
        Option overlapout = new Option("overlapout", true, "overlapout");
        options.addOption(overlapout);
        //output
        Option o = new Option("o", true, "output tsv");
        options.addOption(o);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("utility-name", options);
            System.exit(1);
        }
        //usage info
        if (cmd.getOptions().length == 0 || cmd.hasOption("help")) {
            System.out.println("add usage info");
        }
        //--------------------------------------------------------------------------------------------------------------
        //output
        File outputfile = new File(cmd.getOptionValue("o"));
        FileWriter fw = new FileWriter(outputfile);
        fw.write("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true\n");

        oboReader oboReader = new oboReader(cmd.getOptionValue("obo"), cmd.getOptionValue("root"));

        //File outputfile = new File(cmd.getOptionValue("o"));

        //______________________________________________________________________________________________________
        // TODO: 10.02.2022 -mappingtype go

        MappingReader_new mappingReader_new = new MappingReader_new(cmd.getOptionValue("mapping"), oboReader.getAll_DAGNodes(), cmd.getOptionValue("mappingtype"));


        //______________________________________________________________________________________________________


        //read input of mappingfile
        //MappingReader mappingReader = new MappingReader(cmd.getOptionValue("mapping"), oboReader.getAll_DAGNodes(), cmd.getOptionValue("mappingtype"));

        //System.out.println(mapping1.getGenes().size());

        EnrichReader enrichReader = new EnrichReader(cmd.getOptionValue("enrich"), oboReader.getAll_DAGNodes(), mappingReader_new.getGenes());


        //size
        oboReader.size();

        //shortest path to a true
        int minsizeInt = Integer.parseInt(cmd.getOptionValue("minsize"));
        int maxsizeInt = Integer.parseInt(cmd.getOptionValue("maxsize"));

        oboReader.get_shortest_path(enrichReader.getEnriched_nodes(), minsizeInt, maxsizeInt);


        //pvalues
        int numAllGenes = oboReader.numAllGenes();  //all significant genes
        System.out.println("All significant genes: "+numAllGenes);
        HashSet<Gene> allGenesInTree = oboReader.allGenesInTree();  //all genes that occur in tree


        //--------------------------------------------------------------------------------------------------------------
        //go over all DAG nodes and get pvalues for Benjamini Hochberg
        ArrayList<Double> all_hg_pval = new ArrayList<>();
        ArrayList<Double> all_fej_pval = new ArrayList<>();
        ArrayList<Double> all_ks_pval = new ArrayList<>();
        for (DAGNode d : oboReader.getAll_DAGNodes().values()) {

            if (d.hasCorrectSize(minsizeInt, maxsizeInt)) {
                //System.out.println("d has correct size");
                //------------------------------------------------------------------------------------------------------
                //hg_pval
                HypergeometricDistribution hg = new HypergeometricDistribution(numAllGenes, oboReader.getNumSigGenes(), d.getGenes().size());
                double hg_pval = hg.upperCumulativeProbability(d.getnoverlap());
                all_hg_pval.add(hg_pval);
                d.setHg_pval(hg_pval);
                //------------------------------------------------------------------------------------------------------
                //fej_pval  (wie hg aber fuer alle eingesetzten werte -1)
                double fej_pval = 0.0;
                HypergeometricDistribution fish = new HypergeometricDistribution(numAllGenes - 1, oboReader.getNumSigGenes() - 1, d.getGenes().size() - 1);
                fej_pval = fish.upperCumulativeProbability(d.getnoverlap() - 1);
                all_fej_pval.add(fej_pval);
                d.setFej_pval(fej_pval);
                //------------------------------------------------------------------------------------------------------
                //ks_pval and ks_stat
                KolmogorovSmirnovTest ks = new KolmogorovSmirnovTest();
                ArrayList<Double> in_set = new ArrayList<>();
                ArrayList<Double> bg = new ArrayList<>();
                int count = 0;

                //fill bg distribution & in_set
                count = 0;
                for (Gene g : allGenesInTree) {
                    if (d.getGenes().contains(g)) { //in_set
                        in_set.add(g.getFc());
                    } else {  //bg
                        bg.add(g.getFc());
                    }
                }

                //arraylist to arry
                double[] in_set_distrib = new double[in_set.size()];
                for (Double x : in_set) {
                    in_set_distrib[count] = x;
                    count++;
                }
                count = 0;
                double[] bg_distrib = new double[bg.size()];
                for (Double x : bg) {
                    bg_distrib[count] = x;
                    count++;
                }

                //ks
                double ks_pval = ks.kolmogorovSmirnovTest(in_set_distrib, bg_distrib);
                double ks_stat = ks.kolmogorovSmirnovStatistic(in_set_distrib, bg_distrib);
                all_ks_pval.add(ks_pval);
                d.setKs_pval(ks_pval);
                d.setKs_stat(ks_stat);
            }

        }

        //Benjamini Hochberg
        calculate_BenjaminiHochberg_FDR(all_hg_pval, hg_adjusted);
        calculate_BenjaminiHochberg_FDR(all_fej_pval, fej_adjusted);
        calculate_BenjaminiHochberg_FDR(all_ks_pval, ks_adjusted);


        //write output
        for (DAGNode d : oboReader.getAll_DAGNodes().values()) {

            if (d.hasCorrectSize(Integer.parseInt(cmd.getOptionValue("minsize")), Integer.parseInt(cmd.getOptionValue("maxsize")))) {
                Output out = new Output();

                out.setHg_pval(d.getHg_pval());
                out.setFej_pval(d.getFej_pval());
                out.setKs_pval(d.getKs_pval());
                out.setKs_stat(d.getKs_stat());

                out.setHg_fdr(hg_adjusted.get(d.getHg_pval()));
                out.setFej_fdr(fej_adjusted.get(d.getFej_pval()));
                out.setKs_fdr(ks_adjusted.get(d.getKs_pval()));

                out.setTerm(d.getId());
                out.setName(d.getName());
                out.setSize(d.getSize());
                out.setIs_true(d.isEnriched());
                out.setNoverlap(d.getnoverlap());
                out.setShortest_path_to_a_true(d.getShortest_path_to_a_true());
                fw.write(out.toString());
            }

        }


        fw.close();
        //--------------------------------------------------------------------------------------------------------------

        int counter = 0;
        long startingtime_secondOut = 0;
        if (cmd.hasOption("overlapout")) {
            System.out.println("second out ...");
            startingtime_secondOut = System.currentTimeMillis();


            //todo second output
            File outputfile2 = new File(cmd.getOptionValue("overlapout"));
            FileWriter fw2 = new FileWriter(outputfile2);
            fw2.write("term1\tterm2\tis_relative\tpath_length\tnum_overlapping\tmax_ov_percent\n");

            //go through all genes and get info on DAG entries with shared mapped genes
            HashMap<String, DAGNode> all_DAGNodes = oboReader.getAll_DAGNodes(); // all DAG nodes
            boolean is_relative;
            int path_length;
            int num_overlapping;
            double max_ov_percent;

            HashSet<HashSet<String>> already_seen = new HashSet<>(); //save all valid pairs that were seen before (im output muss jedes Paar genau einmal stehen)


            for (Gene gene : mappingReader_new.getGenes().values()) {
                //System.out.println(gene.getGene_name());

                HashSet<HashSet<String>> validPairs = gene.getAllValidGoPairs(minsizeInt, maxsizeInt, oboReader.getAll_DAGNodes());

                for (HashSet<String> pair : validPairs) { //todo computations for one pair

                    if (already_seen.contains(pair)) { // if a pair was seen before: continue
                        continue;
                    }

                    String[] pair_as_array = new String[2]; // pair hashset to array (damit man ueber index drauf zugreifen kann)
                    pair.toArray(pair_as_array);

                    //is_relative
                     if (all_DAGNodes.get(pair_as_array[0]).getAllParents().contains(pair_as_array[1]) ||
                            all_DAGNodes.get(pair_as_array[1]).getAllParents().contains(pair_as_array[0])) { 

                        is_relative = true;
                    } else {
                        is_relative = false;
                    }


                    //path length
                    path_length = all_DAGNodes.get(pair_as_array[0]).getShortestPath_to_partner(pair_as_array[1], all_DAGNodes);


                    //(a)num_overlapping & (b)max_ov_percent
                    // (a) number of gene ids associated to both DAG entries
                    // (b) maximum percentage (a float value between 0.0 and 100.0) of the
                    //     shared gene ids to all associated gene ids to term1 or term2
                    double[] foo = all_DAGNodes.get(pair_as_array[0]).getNumOverlapping(pair_as_array[1], all_DAGNodes);
                    num_overlapping = (int) foo[0];
                    max_ov_percent = foo[1] * 100.0;


                    Output_second output_second = new Output_second();

                    output_second.setTerm1(pair_as_array[0]);
                    output_second.setTerm2(pair_as_array[1]);
                    output_second.setIs_relative(is_relative);
                    output_second.setPath_length(path_length);
                    output_second.setNum_overlapping(num_overlapping);
                    output_second.setMax_ov_percent(max_ov_percent);

                    counter++;


                    fw2.write(output_second.toString());

                    already_seen.add(pair);

                }


            }

            fw2.close();
        }


        long endingtime = System.currentTimeMillis();
        System.out.println("time total: " + ((endingtime - startingtime) / 1000));
        System.out.println("time second out: " + ((endingtime - startingtime_secondOut) / 1000));



    }


    public static Double[] calculate_BenjaminiHochberg_FDR(ArrayList<Double> pvalues_list, HashMap<Double, Double> adjusted_hashmap) {

        Double[] pvalues = pvalues_list.toArray(new Double[pvalues_list.size()]);

        Double[] adjustedPvalues = new Double[pvalues.length];

        // order the pvalues.
        Arrays.sort(pvalues);

        // iterate through all p-values:  largest to smallest
        for (int i = pvalues.length - 1; i >= 0; i--) {
            if (i == pvalues.length - 1) {
                adjustedPvalues[i] = pvalues[i];
            } else {
                double unadjustedPvalue = pvalues[i];
                int divideByM = i + 1;
                double left = adjustedPvalues[i + 1];
                double right = (pvalues.length / (double) divideByM) * unadjustedPvalue;
                adjustedPvalues[i] = Math.min(left, right);
            }

            adjusted_hashmap.put(pvalues[i], adjustedPvalues[i]);

        }

        return adjustedPvalues;
    }

}
