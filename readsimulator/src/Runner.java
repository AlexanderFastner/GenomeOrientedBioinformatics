import org.apache.commons.cli.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Random;

public class Runner {

    public static void main(String [] args) throws IOException {


        //vars
        //--------------------------------------------------------------------------------------------------------------
        File fastaFile;
        File fidxFile;
        File readcountsFile;
        int readLength;
        int nd;
        int increment = 0;
        boolean foundTransStart;
        boolean foundTransEnd;
        int skippedExonPosition;
        int remainingSequence;
        long genomicStart;
        long genomicEnd;
        double mutationRate;
        double meanFrLength;
        double standardDeviation;
        String fragmentSequence;
        String fw;
        String rw;
        StringBuilder transcriptSequence = new StringBuilder();
        StringBuilder fwMuts = new StringBuilder();
        StringBuilder rwMuts = new StringBuilder();
        ArrayList<Integer> drawnFragments = new ArrayList<>();
        ArrayList<Integer> fwmutationLocation = new ArrayList<>();
        ArrayList<Integer> rwmutationLocation = new ArrayList<>();
        RegionVector fwgenomicRV = new RegionVector();
        RegionVector rwgenomicRV = new RegionVector();

        //--------------------------------------------------------------------------------------------------------------
        //CLP
        Options options = new Options();
        //length
        Option length = new Option("length", true, "input length");
        options.addOption(length);
        //frlength fragment length distribution
        Option frlength = new Option("frlength", true, "mean fragment length");
        options.addOption(frlength);
        //SD standard deviation
        Option SD = new Option("SD", true, "standard deviation");
        options.addOption(SD);
        //mutationrate
        Option mutationrate = new Option("mutationrate", true, "mutationrate");
        options.addOption(mutationrate);
        //readcounts table input
        Option readcounts = new Option("readcounts", true, "read counts table");
        options.addOption(readcounts);
        //fasta
        Option fasta = new Option("fasta", true, "fasta input");
        options.addOption(fasta);
        //fidx fasta index
        Option fidx = new Option("fidx", true, "fasta ids");
        options.addOption(fidx);
        //gtf
        Option gtf = new Option("gtf", true, "gtf");
        options.addOption(gtf);
        //od
        Option outputDir = new Option("od", true, "output dir");
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
        //usage info
        if (cmd.getOptions().length == 0 || cmd.hasOption("help")) {
            System.out.println("""
                    -length <int> : read length
                    • -frlength <int> -SD <int>: fragment length distribution. The input is the mean and standard deviation of a normal distribution. This must be used to draw fragment lengths (you
                    have to round the drawn numbers). For simplicity re-draw the fragment length if the value is
                    smaller than length or larger than the transcript length.
                    1
                    • -readcounts: table of gene_id, transcript_id, count tuples: For every transcript simulate
                    count many reads. For each read pair, select a fragment length (see frlength), and draw a
                    start position from a uniform distribution of all possible start positions on the transcript.
                    • -mutationrate: mutation rate in percent (between 0.0 and 100.0): percent of the simulated
                    bases to be changed from the original.
                    • -fasta: genome FASTA file
                    • -fidx: genome FASTA file index
                    • -gtf: annotation file for the transcript locations
                    • -od: output directory (where the files fw.fastq, rw.fastq, read.mappinginfo will be written)""");
        }
        //--------------------------------------------------------------------------------------------------------------

        fastaFile = new File(cmd.getOptionValue("fasta"));
        fidxFile = new File(cmd.getOptionValue("fidx"));
        readcountsFile = new File(cmd.getOptionValue("readcounts"));
        readLength = Integer.parseInt(cmd.getOptionValue("length"));
        meanFrLength = Integer.parseInt(cmd.getOptionValue("frlength"));
        standardDeviation = Integer.parseInt(cmd.getOptionValue("SD"));
        mutationRate = Double.parseDouble(cmd.getOptionValue("mutationrate"));
//        System.out.println(fastaFile);
//        System.out.println(fidxFile);
//        System.out.println(gtfFile);
//        System.out.println(readcountsFile);
//        System.out.println(readLength);
//        System.out.println(meanFrLength);
//        System.out.println(standardDeviation);
//        System.out.println(mutationRate);

        //create genome
        Genome ge = new Genome();
        //parse
        GTFParser gtfParser = new GTFParser(cmd.getOptionValue("gtf"), ge);
        //raf file
        GenomicSequenceExtractor gse = new GenomicSequenceExtractor(fastaFile, fidxFile);
        //--------------------------------------------------------------------------------------------------------------

        try {
            //----------------------------------------------------------------------------------------------------------
            //output
            File fwfastq = new File(cmd.getOptionValue("od") + "/fw.fastq");
            File rwfastq = new File(cmd.getOptionValue("od") + "/rw.fastq");
            File mappinginfo = new File(cmd.getOptionValue("od") + "/read.mappinginfo");
            FileWriter fwfileWriter = new FileWriter(fwfastq);
            FileWriter rwfileWriter = new FileWriter(rwfastq);
            FileWriter mappinginfofileWriter = new FileWriter(mappinginfo);

            //add header
            mappinginfofileWriter.append("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut").append("\n");
            //----------------------------------------------------------------------------------------------------------

            BufferedReader brReadCounts = new BufferedReader((new FileReader(readcountsFile)));
            //skip header line
            brReadCounts.readLine();
            String line = brReadCounts.readLine();


            //for every entry in readcounts
            while (line != null) {
                //System.out.println(line);
                String[] splitLine = line.split("\t");

                //get start end and chr for transcript

                String chr = ge.findChromosome(splitLine[0]);

                //loop exons and get sequence from fasta
                for (Region exon : ge.getGenome().get(chr).getChromosome().get(splitLine[0]).getExons().get(splitLine[1]).getRegions()) {
//                    if (splitLine[0].equals("ENST00000361361")) {
//                        System.out.println(exon.getStart() + " " + exon.getEnd());
//                    }
                    //combine exons into transcript sequence
                    transcriptSequence.append(gse.getSequence(chr,  exon.getStart(), exon.getEnd()));
                }
                //get reverse complement after all exons read
                if ((ge.getGenome().get(chr).getChromosome().get(splitLine[0]).getStrand().equals("-"))) {
                    transcriptSequence = new StringBuilder(GenomicUtils.reverseComplement(transcriptSequence.toString()));
                }
//                if (splitLine[1].equals("ENST00000361361")) {
//                    System.out.println(splitLine[1]);
//                    String test= transcriptSequence.toString();
//                    int count = 0;
//                    for(char c: test.toCharArray()){
//                        System.out.print(c);
//                        count++;
//                        if (count == 60){
//                            System.out.print("\n");
//                            count = 0;
//                        }
//                    }
//                    System.out.println();
//                }

                //------------------------------------------------------------------------------------------------------

                //ne random
                Random random = new Random();

                //draw n many fragments
                for (int i = 0; i < Integer.parseInt(splitLine[2]); i++) {
                    int FL = 0;
                    //redraw if fragment length is smaller than (input)length or larger than transcript length
                    while (FL < readLength || FL > transcriptSequence.length()) {
                        nd = (int) (random.nextGaussian() * standardDeviation + meanFrLength);
                        FL = Math.max(readLength, nd);
                    }
                    //list of drawn fragment lengths
                    //System.out.println(readLength + " " + FL);
                    drawnFragments.add(FL);
                }
                //------------------------------------------------------------------------------------------------------

                //for each entry in drawnFragments
                //System.out.println("drawnfragmens size " + drawnFragments.size());
                for (int frag : drawnFragments) {

                    //get random start pos
                    int randomPos = (int) (Math.random() * ((transcriptSequence.length() - frag)));
                    //get fragment sequence
                    //System.out.println("random " + randomPos + " :    frag :" + frag + "     length: " + transcriptSequence.length());

                    fragmentSequence = transcriptSequence.substring(randomPos, randomPos + frag);
                    //fragmentSequence = fragmentSequence.replaceAll("\n", "");
//                    System.out.print(fragmentSequence);
//                    System.out.println("F");
                    //System.out.println("fragment sequence: " + fragmentSequence);

                    //get read sequences of length readlength (second read is reverse complement)
                    fw = fragmentSequence.substring(0, readLength);
                    //System.out.println("fw " + fw);
                    rw = fragmentSequence.substring((fragmentSequence.length() - readLength));
                    //System.out.println("rw " + rw);
                    //reverse complement of read
                    rw = GenomicUtils.reverseComplement(rw);
                    //System.out.println("rc of rw " + rw);
                    //System.out.println("**");

                    //--------------------------------------------------------------------------------------------------
                    //simulate mutations
                    Random r = new Random();
                    StringBuilder fwBuiltString = new StringBuilder(fw);
                    StringBuilder rwBuiltString = new StringBuilder(rw);


                    //for each char in fw
                    for (int j = 0; j < fw.length(); j++) {
                        //draw random number from Math.Random() * 100.
                        float probability = r.nextFloat() * 100;
                        //System.out.println(probability);
                        //System.out.println(mutationRate);
                        //if <= mutationRate
                        if (probability <= mutationRate) {
                            //replace that bp with something else
                            fwBuiltString.setCharAt(j, GenomicUtils.mutate(fw.charAt(j)));
                            fwmutationLocation.add(j);
                        }
                    }

                    //same for rw
                    for (int k = 0; k < rw.length(); k++) {
                        //draw random number from Math.Random() * 100.
                        float probability = r.nextFloat() * 100;
                        //System.out.println(probability);
                        //System.out.println(mutationRate);
                        //if <= mutationRate
                        if (probability <= mutationRate) {
                            //replace that bp with something else
                            rwBuiltString.setCharAt(k, GenomicUtils.mutate(rw.charAt(k)));
                            rwmutationLocation.add(k);
                        }
                    }
                    //replace fw with mutated
                    fw = fwBuiltString.toString();
                    //replace rw with mutated
                    rw = rwBuiltString.toString();

                    //convert mutation arrays into comma separated String
                    for (int entry : fwmutationLocation) {
                        fwMuts.append(entry);
                        fwMuts.append(',');
                    }
                    if (fwmutationLocation.size() > 0) {
                        fwMuts.deleteCharAt(fwMuts.length() - 1);
                    }
                    for (int entry : rwmutationLocation) {
                        rwMuts.append(entry);
                        rwMuts.append(',');
                    }
                    if (rwmutationLocation.size() > 0) {
                        rwMuts.deleteCharAt(rwMuts.length() - 1);
                    }
                    //--------------------------------------------------------------------------------------------------
                    //transcript region
                    Region t_fw_regvec = new Region(randomPos, (randomPos + readLength));
                    Region t_rw_regvec = new Region((randomPos + frag - readLength), (randomPos + frag));
                    //System.out.println(t_fw_regvec[0] + " " + t_fw_regvec[1]);
                    //System.out.println(t_rw_regvec[0] + " " + t_rw_regvec[1]);
                    //--------------------------------------------------------------------------------------------------
                    //convert transcript to genomic region

                    //if minus strand change t regions
                    Region tempfw = new Region(t_fw_regvec.getStart(), t_fw_regvec.getEnd());
                    Region temprw = new Region(t_rw_regvec.getStart(), t_rw_regvec.getEnd());
                    if (ge.getGenome().get(chr).getChromosome().get(splitLine[0]).getStrand().equals("-")) {
                        //fw
                        long oldStart = t_fw_regvec.getStart();
                        t_fw_regvec.setStart(transcriptSequence.length() - t_fw_regvec.getEnd());
                        t_fw_regvec.setEnd(transcriptSequence.length() - oldStart);

                        //rw
                        oldStart = t_rw_regvec.getStart();
                        t_rw_regvec.setStart(transcriptSequence.length() - t_rw_regvec.getEnd());
                        t_rw_regvec.setEnd(transcriptSequence.length() - oldStart);
                    }

                    //--------------------------------------------------------------------------------------------------

                    foundTransStart = false;
                    foundTransEnd = false;
                    skippedExonPosition = 0;
                    remainingSequence = (int) t_fw_regvec.getLength() - 1;

//                    if(splitLine[1].equals("ENST00000427231")) {
//                        System.out.println(ge.getGenome().get(chr).getChromosome().get(splitLine[0]).getExons().get(splitLine[1]).getRegions());
//                        System.out.println(ge.getGenome().get(chr).getChromosome().get(splitLine[0]).getStrand());
//                    }

                    //fw
                    //for exons in regions
                    for (Region exons : ge.getGenome().get(chr).getChromosome().get(splitLine[0]).getExons().get(splitLine[1]).getRegions()) {
                        //check start pos
                        if(!foundTransStart) {
                            //if not found
                            //System.out.println(splitLine[1]);
//                            if(splitLine[1].equals("ENST00000427231")) {
//                                System.out.println("exons\t" + exons.getLength() + "\t" + exons);
//                                System.out.println("t_fw_regvec\t" + t_fw_regvec.getStart());
//                                System.out.println("skipped exons\t" + skippedExonPosition);
//                            }
                            if((exons.getStart() + t_fw_regvec.getStart() - skippedExonPosition) <= exons.getEnd()){
                                genomicStart = exons.getStart() + t_fw_regvec.getStart() - skippedExonPosition;
                                foundTransStart = true;
                            }
                            else {
                                //transcript region start isnt on this exon
                                skippedExonPosition += exons.getLength();
                                continue;
                            }
                        }
                        //we found a start but no end yet
                        else{
                            genomicStart = exons.getStart();
                        }
                        //System.out.println("genomic Start\t" + genomicStart);
                        //System.out.println("len" + exons.getLength());

                        //check end pos
                        if(exons.getEnd() >= (genomicStart + remainingSequence - 1)){
                            genomicEnd = genomicStart + remainingSequence;
                            foundTransEnd = true;
                        }
                        //possible off by 1
                        else {
                            genomicEnd = (exons.getEnd() + 1);
                            remainingSequence = (int) (genomicStart + remainingSequence - exons.getEnd() -1);
                        }
                        Region re = new Region(genomicStart, genomicEnd);
                        fwgenomicRV.addRegion(re);
                        //if found both break
                        if(foundTransStart && foundTransEnd){
                            break;
                        }
                    }

                    //System.out.println(fwgenomicRV.getRegions().size());

                    foundTransStart = false;
                    foundTransEnd = false;
                    remainingSequence = (int)t_rw_regvec.getLength() -1;
                    skippedExonPosition = 0;

                    //rw
                    //for exons in regions
                    for (Region exons : ge.getGenome().get(chr).getChromosome().get(splitLine[0]).getExons().get(splitLine[1]).getRegions()) {
                        //check start pos
                        if(!foundTransStart) {
                            //if not found
                            //System.out.println("exons\t" + exons.getLength() + "\t" + exons);
                            if((exons.getStart() + t_rw_regvec.getStart() - skippedExonPosition) <= exons.getEnd()){
                                genomicStart = exons.getStart() + t_rw_regvec.getStart() - skippedExonPosition;
                                foundTransStart = true;
                            }
                            else {
                                //transcript region start isnt on this exon
                                skippedExonPosition += exons.getLength();
                                continue;
                            }
                        }
                        //we found a start but no end yet
                        else{
                            genomicStart = exons.getStart();
                        }
                        //System.out.println("len" + exons.getLength());
                        //check end pos
                        if(exons.getEnd() >= (genomicStart + remainingSequence - 1)){
                            genomicEnd = genomicStart + remainingSequence;
                            foundTransEnd = true;
                        }
                        else {
                            genomicEnd = exons.getEnd() + 1;
                            remainingSequence = (int) (genomicStart + remainingSequence - exons.getEnd() -1);
                        }
                        Region re = new Region(genomicStart, genomicEnd);
                        rwgenomicRV.addRegion(re);
                        //if found both break
                        if(foundTransStart && foundTransEnd){
                            break;
                        }

                    }

                    //--------------------------------------------------------------------------------------------------
                    //fw
                    StringBuilder fwGout = new StringBuilder();
                    //format output
                    for(Region fwregion : fwgenomicRV.getRegions()){
                        fwGout.append(fwregion.getStart()).append("-").append(fwregion.getEnd()).append("|");
                    }
                    //remove extra |
                    fwGout = new StringBuilder(fwGout.substring(0, fwGout.length() - 1));

                    //rw
                    StringBuilder rwGout = new StringBuilder();
                    for(Region rwregion: rwgenomicRV.getRegions()){
                        rwGout.append(rwregion.getStart()).append("-").append(rwregion.getEnd()).append("|");
                    }
                    rwGout = new StringBuilder(rwGout.substring(0, rwGout.length() - 1));
                    //System.out.println(increment + " " + fw.length() + " : " + fw);
                    //System.out.println(increment + " " + rw.length() + " : " + rw);
                    //System.out.println(increment);
                    //--------------------------------------------------------------------------------------------------
                    //output
                    //fw
                    fwfileWriter.append("@").append(String.valueOf(increment)).append("\n");
                    fwfileWriter.append(fw).append("\n");
                    fwfileWriter.append("+").append(String.valueOf(increment)).append("\n");
                    fwfileWriter.append("I".repeat(fw.length())).append("\n");
                    //rw
                    rwfileWriter.append("@").append(String.valueOf(increment)).append("\n");
                    rwfileWriter.append(rw).append("\n");
                    rwfileWriter.append("+").append(String.valueOf(increment)).append("\n");
                    rwfileWriter.append("I".repeat(rw.length())).append("\n");

                    //mappinginfo
                    //readid	chr	gene	transcript	t_fw_regvec	t_rw_regvec	fw_regvec	rw_regvec	fw_mut	rw_mut
                    mappinginfofileWriter.append(String.valueOf(increment)).append("\t").append(chr).append("\t")
                        .append(splitLine[0]).append("\t").append(splitLine[1]).append("\t")
                        .append(String.valueOf(tempfw.getStart())).append("-").append(String.valueOf(tempfw.getEnd())).append("\t")
                        .append(String.valueOf(temprw.getStart())).append("-").append(String.valueOf(temprw.getEnd())).append("\t")
                        .append(fwGout.toString()).append("\t")
                        .append(rwGout.toString()).append("\t")
                        .append(fwMuts.toString()).append("\t").append(rwMuts.toString())
                        .append("\n");

                    //increase increment counter
                    increment++;
                    //clear stuff
                    fwMuts.setLength(0);
                    rwMuts.setLength(0);
                    fwmutationLocation.clear();
                    rwmutationLocation.clear();
                    fwgenomicRV = new RegionVector();
                    rwgenomicRV = new RegionVector();
                }

                //clear drawn fragments
                drawnFragments.clear();

                line = brReadCounts.readLine();
                //reset transcript sequence for next line
                transcriptSequence.setLength(0);
            }
            fwfileWriter.close();
            rwfileWriter.close();
            mappinginfofileWriter.close();

        }catch (Exception e){
            e.printStackTrace();
            System.out.println("not good");
        }
    }
}
