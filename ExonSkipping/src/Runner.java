import org.apache.commons.cli.*;
import java.io.*;
import java.util.*;

public class Runner {

    //make output file
    static HashSet<ESoutput> outputLine = new HashSet<>();

    public static void main(String [] args) {
        Options options = new Options();
        //-gtf
        Option gtf = new Option("gtf", true, "gtf input");
        options.addOption(gtf);
        //-o
        Option outputDir = new Option("o", true, "output path");
        options.addOption(outputDir);

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
        //print usage info
        if (cmd.getOptions().length < 2 || cmd.hasOption("help")) {
            System.out.println("Please add the following parameters:\n -gtf <gtf_file>\n" +
                    " -o <output_tsv>\n");
            return;
        }

        //--------------------------------------------------------------------------------------------------------------
        long startTime = System.currentTimeMillis();
        File outputfile = new File(cmd.getOptionValue("o"));

        try {
            //System.out.println(cmd.getOptionValue("o"));
            FileWriter fw = new FileWriter(cmd.getOptionValue("o"));
            //BufferedReader out = new BufferedReader(new FileReader(outputfile));
            //Headers
            fw.write("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tSV_prots\tWT_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases\n");
            fw.close();

            //create genome
            Genome ge = new Genome();

            //GTFParser
            GTFParser gtfParser = new GTFParser(cmd.getOptionValue("gtf"), ge);

            //----------------------------------------------------------------------------------------------------------
            //exon skipping logic


            //for every chr in genome
            for (var chrEntry : ge.getGenome().entrySet()){
                //String logicchr = chrEntry.getKey();
                //for every gene in chr
                for (var geneEntry : chrEntry.getValue().getChromosome().entrySet()) {
                    //String logicgene = geneEntry.getKey();

                    HashSet<Region> all_SV =new HashSet<>();
                    HashMap<String, RegionVector> all_regionVectors_in_transcript = geneEntry.getValue().getRegionVectors();
                    RegionVector all_introns_in_gene = new RegionVector();
                    RegionVector all_exons_in_gene = new RegionVector();

                    //for every transcript in gene
                    all_regionVectors_in_transcript.forEach((transID, regionVector) -> {
                    //for (var transID : all_regionVectors_in_transcript.entrySet()){
                        // all exons of this transcript
                        RegionVector all_cds_in_transcript = geneEntry.getValue().CDS.get(transID);

                        //add to all_introns_in_gene
                        //System.out.println(all_exons_in_transcript);
                        for (Region r : all_cds_in_transcript.getRegions()) {
                            all_exons_in_gene.addRegion(r);
                        }

                        RegionVector all_introns_in_transcript = regionVector.generateIntrons();

                        //System.out.println(all_introns_in_gene.getRegions().size());
                        for (Region r : all_introns_in_transcript.getRegions()) {
                            all_introns_in_gene.addRegion(r);
                            //System.out.println("x");
                        }
                        //System.out.println(all_introns_in_gene.getRegions().size());
                    });

                    //System.out.println(all_introns_in_gene.getRegions());

                    RegionVector SV = all_introns_in_gene.get_SV();

                    //now loop through all found SV

                    //System.out.println(SV.getRegions().size());
                    for (Region sv_intron : SV.getRegions()) {
                        if (all_SV.contains(sv_intron)) {
                            continue;
                        }
                        all_SV.add(sv_intron);

                        //Treesets to store output vars
                        TreeSet<Integer> num_skipped_exons = new TreeSet<>();
                        TreeSet<Integer> len_skipped_exons = new TreeSet<>();

                        //all transID for SV
                        HashSet<String> sv_transcript_ids = new HashSet<>();
                        //all transID for WT
                        HashSet<String> wt_transcript_ids = new HashSet<>();


                        //loop all transcripts to check if SV/WT transcripts
                        for(var transcriptID : all_regionVectors_in_transcript.entrySet()) {

                            RegionVector intron_RV = transcriptID.getValue().generateIntrons();

                            //if transcript has SV introns - its SV
                            if (intron_RV.getRegions().contains(sv_intron)) {
                                sv_transcript_ids.add(transcriptID.getKey());
                            }
                            //if it has WT its WT
                            else if (intron_RV.has_wt_introns(sv_intron)) {
                                wt_transcript_ids.add(transcriptID.getKey());
                            }

                        }
                        //System.out.println(sv_transcript_ids.size());
                        //System.out.println(wt_transcript_ids.size());
                        if (wt_transcript_ids.isEmpty()) {
                            continue;
                        }
                        //----------------------------------------------------------------------------------------------
                        //SV_protIDs
                        TreeSet<String> sv_prots = new TreeSet<>();
                        //loop SVtransID
                        for (String sv_trans_id : sv_transcript_ids) {
                            try {
                                for (String prot_id : geneEntry.getValue().getAllTranscripts().get(sv_trans_id).getProtIds()) {
                                    //System.out.println("here");
                                    sv_prots.add(prot_id);
                                    //System.out.println(prot_id);
                                }
                            } catch (Exception e) {
                                System.out.println("SV prots error");
                            }
                        }
                        //----------------------------------------------------------------------------------------------
                        //WT
                        TreeSet<String> wt_prots = new TreeSet<>();
                        RegionVector wt_introns = new RegionVector();

                        //loop WTtransIDs
                        for (String wt_trans_id : wt_transcript_ids) {
                            //protIDs for WT trans
                            try {
                                //System.out.println("test")
                                //System.out.println("WT prots " + geneEntry.getValue().getAllTranscripts().get(wt_transcript_id).getProtIds());
                                for (String prot_id : geneEntry.getValue().getAllTranscripts().get(wt_trans_id).getProtIds()) {
                                    wt_prots.add(prot_id);

                                }
                            } catch (Exception e) {
                                System.out.println("WT prots error");
                            }

                            //
                            int num_WT_introns = 0;
                            int len_WT_introns = 0;
                            TreeSet<Region> WT_trans_introns = geneEntry.getValue().CDS.get(wt_trans_id).generateIntrons().getRegions();

                            //get introns from transID
                            for (Region reg : WT_trans_introns) {
                                if (sv_intron.getStart() <= reg.getStart() && sv_intron.getEnd() >= reg.getEnd()) {
                                    //This is a WT intron
                                    wt_introns.addRegion(reg);
                                    num_WT_introns++;
                                    len_WT_introns += reg.getLength();
                                }

                            }

                            //if more than 2 introns we Skipped an Exon
                            if (num_WT_introns >= 2) {
                                //subtract 1 cause cause less exons in this case
                                int num_skipped_exons_in_trans = num_WT_introns - 1;
                                int len_skipped_exons_in_trans = sv_intron.getLength() - len_WT_introns;
                                //add to output
                                num_skipped_exons.add(num_skipped_exons_in_trans);
                                len_skipped_exons.add(len_skipped_exons_in_trans);
                            }

                        }
                        //System.out.println(wt_prots);


                        //----------------------------------------------------------------------------------------------
                        //output

                        //GeneID, Name, Chr, Strand, nProt, nTrans
                        Gene output_Gene = geneEntry.getValue();
                        if (output_Gene != null) {
                            //output
                            ESoutput esoutput = new ESoutput();
                            esoutput.setGene_id(output_Gene.getId());
                            esoutput.setGene_symbol(output_Gene.getName());
                            esoutput.setChr(chrEntry.getKey());
                            esoutput.setStrand(output_Gene.getStrand());
                            esoutput.setNprots(output_Gene.getCDS().size());
                            esoutput.setNtrans(output_Gene.getAllTranscripts().size());
                            esoutput.setSV_intron(sv_intron.toString());
                            esoutput.setWT_introns(wt_introns.toString());
                            esoutput.setSV_prots(String.join("|", sv_prots));
                            esoutput.setWT_prots(String.join("|",wt_prots));
                            esoutput.setMin_skippped_exons(Collections.min(num_skipped_exons));
                            esoutput.setMax_skipped_exons(Collections.max(num_skipped_exons));
                            esoutput.setMin_skipped_bases(Collections.min(len_skipped_exons));
                            esoutput.setMax_skipped_bases(Collections.max(len_skipped_exons));

//                            //todo testing why 2 same objects
//                            if (output_Gene.getId().equals("ENSG00000113761")){
//                                System.out.println(outputLine.size());
//                                System.out.println(esoutput);
//                                ESoutput test = new ESoutput();
//                                test.setGene_id(output_Gene.getId());
//                                test.setGene_symbol(output_Gene.getName());
//                                test.setChr(chrEntry.getKey());
//                                test.setStrand(output_Gene.getStrand());
//                                test.setNprots(output_Gene.getCDS().size());
//                                test.setNtrans(output_Gene.getAllTranscripts().size());
//                                test.setSV_intron(sv_intron.toString());
//                                test.setWT_introns(wt_introns.toString());
//                                test.setSV_prots(String.join("|", sv_prots));
//                                test.setWT_prots(String.join("|",wt_prots));
//                                test.setMin_skippped_exons(Collections.min(num_skipped_exons));
//                                test.setMax_skipped_exons(Collections.max(num_skipped_exons));
//                                test.setMin_skipped_bases(Collections.min(len_skipped_exons));
//                                test.setMax_skipped_bases(Collections.max(len_skipped_exons));
//                                System.out.println(test);
//                                outputLine.add(test);
//                                System.out.println(outputLine.size());
//                            }


                            outputLine.add(esoutput);

                            //old line by line output
                            //replaced to be able to remove duplicates
//                            output.add(output_Gene.getId() + "\t" + output_Gene.getName() +
//                                    "\t" + chrEntry.getValue().getChrID() + "\t" + output_Gene.getStrand() + "\t" +
//                                    output_Gene.getCDS().size() + "\t" + output_Gene.getAllTranscripts().size() + "\t" +
//                                    sv_intron + "\t" + wt_introns +
//                                    "\t" + String.join("|", sv_prots) + "\t" + String.join("|",wt_prots) + "\t" +
//                                    "\t" + Collections.min(num_skipped_exons) + "\t" + Collections.max(num_skipped_exons) +
//                                    "\t" + Collections.min(len_skipped_exons) + "\t" + Collections.max(len_skipped_exons) + "\n");
                        }

                    }

                }

            }
        }catch (IOException e) {
            e.printStackTrace();
        }
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(outputfile, true));
            for(ESoutput lines: outputLine){
                bw.write(lines + "\n");
            }
            bw.close();
        }
        catch (IOException e) {
            System.out.println("unexpected writing error");
            e.printStackTrace();
        }

        long EndTime = System.currentTimeMillis();

        //System.out.println("Time: " + (EndTime - startTime));

    }

}
