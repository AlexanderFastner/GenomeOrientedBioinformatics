import augmentedTree.IntervalTree;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.apache.commons.cli.*;
import java.io.*;
import java.util.*;

public class Runner {

    public static void main(String [] args) throws IOException, OutOfMemoryError {

        long startingtime = System.currentTimeMillis();

        //--------------------------------------------------------------------------------------------------------------
        //CLP
        Options options = new Options();
        //gtf
        Option gtf = new Option("gtf", true, "gtf file");
        options.addOption(gtf);
        //bam
        Option bam = new Option("bam", true, "bam file input");
        options.addOption(bam);
        //o
        Option outputtsv = new Option("o", true, "output tsv");
        options.addOption(outputtsv);
        //frstrand
        Option frstrand = new Option("frstrand", true, "frstrand true/false");
        options.addOption(frstrand);

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
            System.out.println("-gtf: enter input gtf file\n-bam: enter the input bam file\n-o: name of output tsv\n-frstrand: strandedness of experiment");
        }
        //--------------------------------------------------------------------------------------------------------------

        File bamfile = new File(cmd.getOptionValue("bam"));
        ArrayList<bamOutput> outputs = new ArrayList<bamOutput>();
        HashMap<String, SAMRecord> lookup = new HashMap<String, SAMRecord>();

        Genome ge = new Genome();
        //gtf
        GTFParser gtfParser = new GTFParser(cmd.getOptionValue("gtf"), ge);

        //read in bam file
        SAMFileReader samreader = new SAMFileReader(bamfile, false);
        samreader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        Iterator<SAMRecord> it = samreader.iterator();

        //output
        File outputFile = new File(cmd.getOptionValue("o"));
        FileWriter fw = new FileWriter(outputFile);

        File rpkMall = new File("RPKMall");
        FileWriter rpkm = new FileWriter(rpkMall);

        File rpkmpcr = new File("RPKMPCRO");
        FileWriter rpkmpcrO = new FileWriter(rpkmpcr);

        //iterate through chromosomes to get all the trees
        HashMap<String, IntervalTree<Region>> intervaltreeP = new HashMap<>();
        HashMap<String, IntervalTree<Region>> intervaltreeN = new HashMap<>();

        HashMap<String, Integer> tmicount = new HashMap<>();
        HashMap<String, Integer> pcrzero = new HashMap<>();

        //fill hashmap with interval trees
        for (Chromosome c : ge.getGenome().values()) {
            intervaltreeP.put(c.getChrID(), c.getitP());
            intervaltreeN.put(c.getChrID(), c.getitN());
//            System.out.println("chr: " + c.getChrID() + "\tInterval tree P " + intervaltreeP);
//            System.out.println("chr: " + c.getChrID() + "\tInterval tree N " + intervaltreeN);
//            System.out.println();
        }
//        for(Chromosome ce: ge.getGenome().values()){
//            System.out.println(intervaltreeP.get(ce.getChrID()));
//        }
        //System.out.println(intervaltreeP);

        //vars
        int clipping = 0;
        int starting = 0;
        int ending = 0;
        String currentChr = "";
        SAMRecord fwread;
        //hashmap of all merged RVs
        HashMap<String, Integer> mergedRVs = new HashMap<>();
        ArrayList <String> cleanup = new ArrayList<String>();
        int mergedRVcleanup = 30000;

        //iterate all bam entries
        while (it.hasNext()) {
            SAMRecord sr = it.next();
            //clear lookup if chr changes
            if (!currentChr.equals(sr.getReferenceName())) {
                System.out.println("Chromosome " + sr.getReferenceName());
                System.out.println(System.currentTimeMillis() - startingtime);
//                System.out.println(tmicount.size());
//                System.out.println(pcrzero.size());
                System.out.println();
                //clear lookup
                lookup.clear();
                //clear hashmap of all merged RVs
                mergedRVs.clear();
            }
            //get otherseen from lookup
            SAMRecord otherseen = lookup.get(sr.getReadName());

            //testing
//            System.out.println("ReadName:" + sr.getReadName());
//            System.out.println("NotPrimaryAlignmentFlag:" +sr.getNotPrimaryAlignmentFlag());
//            System.out.println("ReferenceName:" + sr.getReferenceName());
//            System.out.println("ReadPairedFlag:" + sr.getReadPairedFlag());
//            System.out.println("ReadNegativeStrandFlag:" + sr.getReadNegativeStrandFlag());
//            System.out.println("MateUnmappedFlag:" + sr.getMateUnmappedFlag());
//            System.out.println("FirstOfPairFlag:" + sr.getFirstOfPairFlag());

            //loop till both reads of pair have been found
            if (otherseen != null) {
                RegionVector mergedReads = new RegionVector();
                RegionVector mergedIntronReads = new RegionVector();
                int pcrindex = 0;

                if (sr.getFirstOfPairFlag()) {
                    fwread = sr;
                } else {
                    fwread = otherseen;
                }

                //System.out.println(otherseen.getReadName());

                //output
                bamOutput bamO = new bamOutput(Integer.parseInt(otherseen.getReadName()));

                //handle read pair
                //------------------------------------------------------------------------------------------------------
                //precheck
                //get smaller start and bigger end from both genes
                starting = Math.min(sr.getAlignmentStart(), otherseen.getAlignmentStart());
                ending = Math.max(sr.getAlignmentEnd(), otherseen.getAlignmentEnd());

//                if(sr.getReadName().equals("50078")){
//                    System.out.println(sr.getAlignmentStart());
//                    System.out.println(otherseen.getAlignmentStart());
//                    System.out.println(starting);
//                    System.out.println(sr.getAlignmentEnd());
//                    System.out.println(otherseen.getAlignmentEnd());
//                    System.out.println(ending);
//                }

                IntervalTree<Region> correctTree;
                IntervalTree<Region> otherTree_for_antisense = null;

                if (cmd.getOptionValue("frstrand") == null) {
                    correctTree = intervaltreeP.get(sr.getReferenceName());
                    correctTree.addAll(intervaltreeN.get(sr.getReferenceName()));
                } else if ((cmd.getOptionValue("frstrand").equals("false") && fwread.getReadNegativeStrandFlag()) || (cmd.getOptionValue("frstrand").equals("true") && !fwread.getReadNegativeStrandFlag())) {
                    correctTree = intervaltreeP.get(sr.getReferenceName());
                    //for antisnese
                    otherTree_for_antisense = intervaltreeN.get(sr.getReferenceName());
                } else {
                    correctTree = intervaltreeN.get(sr.getReferenceName());
                    //for antisense
                    otherTree_for_antisense = intervaltreeP.get(sr.getReferenceName());
                }

                //------------------------------------------------------------------------------------------------------
                //antisense
                boolean antisense = false;
                if (!(otherTree_for_antisense == null)) {
                    ArrayList<Region> cgenes_other = otherTree_for_antisense.getIntervalsSpanning(starting, ending, new ArrayList<Region>());

                    if (cgenes_other.size() > 0) {
                        antisense = true;
                    }
                    cgenes_other.clear();
                }
                bamO.setAntisense(antisense);

                //------------------------------------------------------------------------------------------------------
//                if(sr.getReadName().equals("50078")) {
//                    System.out.println("start and end:" + starting + " " + ending);
//                }

                ArrayList<Region> cgenes = correctTree.getIntervalsSpanning(starting, ending, new ArrayList<Region>());
                ArrayList<Region> igenes = correctTree.getIntervalsSpannedBy(starting, ending, new ArrayList<Region>());

//                System.out.println("igenes " + igenes);
//                System.out.println("cgenes " + cgenes);
                //cgenes - ein teil vom read ist gene
                //getIntervalsSpanning
                //start stop from start of first to end of second
                //get smaller start and bigger end from both genes
                //igenes - der read includes 1 whole gene
                //getIntervalsSpannedBy

                //pre check (useful to remove from memory early in program)
                //if |igenes| > 0 && |cgenes| == 0
                if (igenes.size() > 0 && cgenes.size() == 0) {
                    //System.out.println("this is one we can skip + " + otherseen.getReadName());
                    continue;
                }
                //------------------------------------------------------------------------------------------------------
                //mm
                Integer nm = (Integer) sr.getAttribute("NM");
                nm = (nm != null) ? nm : (Integer) sr.getAttribute("nM");
                nm = (nm != null) ? nm : (Integer) sr.getAttribute("XM");
                Integer nmo = (Integer) otherseen.getAttribute("NM");
                nmo = (nmo != null) ? nmo : (Integer) otherseen.getAttribute("nM");
                nmo = (nmo != null) ? nmo : (Integer) otherseen.getAttribute("XM");
                bamO.setMm(nm + nmo);
                //------------------------------------------------------------------------------------------------------
                //clipping size
                //get difference between getallignmentX and getUnclippedX

                clipping = (((sr.getAlignmentStart() - sr.getUnclippedStart()) + (sr.getUnclippedEnd() - sr.getAlignmentEnd()))
                        + ((otherseen.getAlignmentStart() - otherseen.getUnclippedStart()) + (otherseen.getUnclippedEnd() - otherseen.getAlignmentEnd())));

                bamO.setClippingSize(clipping);

                //------------------------------------------------------------------------------------------------------
                //split count
                int splitcount = 0;
                boolean splitinconsistency = false;

                //get exons from alignment blocks
                RegionVector exonssr = exonsFromAB(sr.getAlignmentBlocks()).remove_overlap();
                RegionVector exonsotherseen = exonsFromAB(otherseen.getAlignmentBlocks()).remove_overlap();
                //convert to introns
                RegionVector intronssr = exonssr.generateIntrons();
                RegionVector intronsotherseen = exonsotherseen.generateIntrons();

                boolean foundsc = false;
                //if at least 1 read has implied introns
                if (exonssr.regions.size() > 1 || exonsotherseen.regions.size() > 1) {

                    int smallerend = Math.min(exonssr.regions.last().getStop(), exonsotherseen.regions.last().getStop());
                    int biggerstart = Math.max(exonssr.regions.first().getStart(), exonsotherseen.regions.first().getStart());

                    //No overlap
                    if (smallerend < biggerstart) {
                        splitcount = intronssr.regions.size() + intronsotherseen.regions.size();
                        bamO.setNsplit(splitcount);
                        foundsc = true;

                    }
                    //Overlap
                    if (!foundsc) {
                        int startoverlap = biggerstart;
                        int endoverlap = smallerend;
                        Region overlap = new Region(startoverlap, endoverlap);
                        RegionVector intronsOverlapSr = intronsinOverlap(overlap, intronssr);
                        RegionVector intronsOverlapOtherseen = intronsinOverlap(overlap, intronsotherseen);

                        if (intronsOverlapSr.regions.equals(intronsOverlapOtherseen.regions)) {
                            splitcount = intronssr.regions.size() + intronsotherseen.regions.size() - intronsOverlapSr.regions.size();
                            bamO.setNsplit(splitcount);
                            foundsc = true;
                        }
                        else {
                            //System.out.println("splitinconsistency");
                            splitinconsistency = true;
                            bamO.setNsplit(-1);
                        }
                    }
                }
                //------------------------------------------------------------------------------------------------------
                //GDIST
                int gdist = 0;
                boolean is_intergenic = false;
                if(!splitinconsistency){
                    //if intergenic
                    if(cgenes.size() == 0){
                        is_intergenic = true;
                        bamO.setIntergenic(is_intergenic);
                        ArrayList<Region> leftneighbor = new ArrayList<>();
                        ArrayList<Region> rightneighbor = new ArrayList<>();
                        leftneighbor = correctTree.getIntervalsLeftNeighbor(starting, ending, new ArrayList<Region>());
                        rightneighbor = correctTree.getIntervalsRightNeighbor(starting, ending, new ArrayList<Region>());

                        boolean neighborsequal = false;
                        //check if left and right neighbor are equal, if yes gdist = 0
                        if (leftneighbor.equals(rightneighbor)) {
                            gdist = 0;
                            neighborsequal = true;
                            //continue;
                        }
                        if (!neighborsequal) {
                            int distLeft;
                            if (leftneighbor.size() == 0) {
                                distLeft = Integer.MAX_VALUE;
                            }
                            //get distance of start and end of left neighbor to fragStart and fragEnd and select min
                            else {
                                int leftend = leftneighbor.get(0).getStop();
                                int distLeft_a = Math.min(Math.abs(starting - leftend), Math.abs(ending - leftend));

                                leftend = leftneighbor.get(0).getStart();
                                int distLeft_b = Math.min(Math.abs(starting - leftend), Math.abs(ending - leftend));

                                distLeft = Math.min(distLeft_a, distLeft_b);

                            }
                            int distRight;
                            if (rightneighbor.size() == 0) {
                                distRight = Integer.MAX_VALUE;
                            } else {    //get distance of start and end of right neighbor to fragStart and fragEnd and select min
                                int right_start = rightneighbor.get(0).getStart();
                                int distRight_a = Math.min(Math.abs(starting - right_start), Math.abs(ending - right_start));
                                right_start = rightneighbor.get(0).getStop();
                                int distRight_b = Math.min(Math.abs(starting - right_start), Math.abs(ending - right_start));

                                distRight = Math.min(distRight_a, distRight_b);
                            }
                            // get min of left and right neighbor
                            gdist = Math.min(distLeft, distRight) - 1;
                        }
                    }
                }
                bamO.setGdist(gdist);

                //------------------------------------------------------------------------------------------------------
                //transcriptomic/merged transcriptomic/intronic

                StringBuilder transstring = new StringBuilder();
                boolean transcriptomicfound = false;
                boolean mergedfound = false;
                //for each transcript check if transcriptomic

                //key = geneID, value = transcriptIDs if transcriptomic || MERGED || INTRONIC
                HashMap<String, ArrayList<String>> transcriptsbyGene = new HashMap<>();

                //merge transcripts from sr and otherseen
                for (Region r : exonssr.regions) {
                    mergedReads.addRegion(new Region(r.getStart(), r.getStop()));
                }
                for (Region r : exonsotherseen.regions) {
                    mergedReads.addRegion(new Region(r.getStart(), r.getStop()));
                }
                //merged introns
                for (Region r : intronssr.regions) {
                    mergedIntronReads.addRegion(new Region(r.getStart(), r.getStop()));
                }
                for (Region r : intronsotherseen.regions) {
                    mergedIntronReads.addRegion(new Region(r.getStart(), r.getStop()));
                }
                mergedReads = mergedReads.remove_overlap();

                //------------------------------------------------------------------------------------------------------
                //pcr index
                //for pcrindex

                String combined = mergedReads.toString();
//                System.out.println(mergedReads);

                if(cmd.getOptionValue("frstrand") != null) {
                    if (fwread.getReadNegativeStrandFlag()) {
                        combined += '-';
                    }
                    else {
                        combined += '+';
                    }
                }

//                if(sr.getReadName().equals("3202142")){
//                    System.out.println(combined);
//                    System.out.println(mergedRVs);
//                }

//                System.out.println("MergedRVS " + mergedRVs.size());
                //calculate PCRindex
                //for all RV in hashmap
                if(mergedRVs.containsKey(combined)){
                    pcrindex = mergedRVs.get(combined);
                    mergedRVs.put(combined, pcrindex + 1);
                }
                else{
                    mergedRVs.put(combined, 1);
                }
                if (splitinconsistency){
                    mergedRVs.put(combined, pcrindex);
                }
                bamO.setPcrindex(pcrindex);
//                System.out.println(mergedRVs.size());

                //for all cgenes
                for(Region geneSpan : cgenes){
                    //RPKM
                    if(tmicount.containsKey(geneSpan.getG().getId())){
                        tmicount.put(geneSpan.getG().getId(), tmicount.get(geneSpan.getG().getId())+1);
                    }
                    else {
                        tmicount.put(geneSpan.getG().getId(), 1);
                    }

                    if(pcrzero.containsKey(geneSpan.getG().getId())) {
                        if (pcrindex == 0) {
                            pcrzero.put(geneSpan.getG().getId(), pcrzero.get(geneSpan.getG().getId()) + 1);
                        }
                    }
                    else{
                        pcrzero.put(geneSpan.getG().getId(),1);
                    }


                    //cleanup merged RV
                    if(mergedRVcleanup < mergedRVs.size()){
                        //loop Merged RV and check if ends are < current start
                        for(String s: mergedRVs.keySet()){
                            int start = mergedReads.getStart();
//                            System.out.println(s);
                            int mergedend = Integer.parseInt(s.substring(s.lastIndexOf(":") + 1, s.length()-1));
//                            System.out.println(mergedend);
                            if(start > mergedend) {
//                                System.out.println(start);
//                                System.out.println(mergedend);
                                cleanup.add(s);
                            }
                        }
//                        System.out.println(mergedRVs.size());
                        for(String s: cleanup){
                            mergedRVs.remove(s);
                        }
//                        System.out.println(mergedRVs.size());
                        cleanup.clear();
                    }

                    //--------------------------------------------------------------------------------------------------


                    //get transcripts from this gene
                    //System.out.println(geneSpan);
//                    System.out.println(correctTree);
                    HashMap<String, RegionVector> transcripts = null;
                    ArrayList<Region> foundGenes = correctTree.getIntervalsEqual(geneSpan.getStart(), geneSpan.getStop(), new ArrayList<Region>());
                    for(Region region: foundGenes) {
                        transcriptomicfound = false;
                        mergedfound = false;

                        Gene curGene = region.getG();
                        transcripts = curGene.getExons();

                        //for each transcript check if transcriptomic
                        for (Map.Entry<String, RegionVector> regionVectorEntry : transcripts.entrySet()) {
                            //check if transcriptomic

                            boolean srfound = false;
                            boolean otherseenfound = false;

                            exonssr = exonssr.remove_overlap();
                            exonsotherseen = exonsotherseen.remove_overlap();

                            if (regionVectorEntry.getValue().subRV(exonssr)) {
                                srfound = true;
                            }
                            if (regionVectorEntry.getValue().subRV(exonsotherseen)) {
                                otherseenfound = true;
                            }
                            if (srfound && otherseenfound) {
                                transcriptomicfound = true;
                                try {
                                    //check if exists
                                    transcriptsbyGene.get(curGene.getId());
                                    transcriptsbyGene.get(curGene.getId()).add(regionVectorEntry.getKey());

                                } catch (Exception e) {
                                    ArrayList<String> transIdslist = new ArrayList<>();
                                    transIdslist.add(regionVectorEntry.getKey());
                                    transcriptsbyGene.put(curGene.getId(), transIdslist);
                                }

                            }

                        }

                        if (!transcriptomicfound) {
                            //------------------------------------------------------------------------------------------------------
                            RegionVector mergedTranscript = new RegionVector();
                            for (RegionVector rv : transcripts.values()) {
                                for (Region r : rv.regions) {
                                    mergedTranscript.addRegion(new Region(r.getStart(), r.getStop()));
                                }
                            }
                            mergedTranscript = mergedTranscript.remove_overlap();

                            //------------------------------------------------------------------------------------------------------
                            //compare to merged transcript
                            if (mergedTranscript.contained(mergedReads)) {
                                //if match stop and skip intronic
                                mergedfound = true;
                                try {
                                    transcriptsbyGene.get(curGene.getId());
                                    transcriptsbyGene.get(curGene.getId()).add("MERGED");
                                } catch (Exception e) {
                                    ArrayList<String> transIdslist = new ArrayList<>();
                                    transIdslist.add("MERGED");
                                    transcriptsbyGene.put(curGene.getId(), transIdslist);
                                }
                            }

                            //else its intronic
                            if (!mergedfound && !transcriptomicfound) {
                                try {
                                    transcriptsbyGene.get(curGene.getId());
                                    transcriptsbyGene.get(curGene.getId()).add("INTRON");
                                } catch (Exception e) {
                                    ArrayList<String> transIdslist = new ArrayList<>();
                                    transIdslist.add("INTRON");
                                    transcriptsbyGene.put(curGene.getId(), transIdslist);
                                }
                            }
                        }
                    }


                    //output
                    if(!splitinconsistency) {

                        boolean is_transcriptomic = false;
                        boolean is_merged = false;

                        //only print highest level of tmi
                        for(Map.Entry<String, ArrayList<String>> entry: transcriptsbyGene.entrySet()){
                            ArrayList <String> undeterminedTMI = entry.getValue();
                            if (!undeterminedTMI.get(0).equals("MERGED") && !undeterminedTMI.get(0).equals("INTRON")) {
                                is_transcriptomic = true;
                            } else if (undeterminedTMI.get(0).equals("MERGED")) {
                                is_merged = true;
                            }
                        }

                        int gcount = 0;
                        StringBuilder tmi = new StringBuilder();
                        for(Map.Entry<String, ArrayList<String>> entry: transcriptsbyGene.entrySet()){
                            ArrayList <String> undeterminedTMI = entry.getValue();
                            if(is_transcriptomic){
                                if (!undeterminedTMI.get(0).equals("MERGED") && !undeterminedTMI.get(0).equals("INTRON")) {
                                    gcount++;
                                    tmi.append(entry.getKey()).append(",").append(ge.getGenome().get(currentChr).getChromosome().get(entry.getKey()).getType()).append(":");
                                    for (String t_id : undeterminedTMI) {
                                        tmi.append(t_id).append(",");
                                    }
                                }
                                //remove extra comma
                                if (tmi.toString().endsWith(",")) {
                                    tmi = new StringBuilder(tmi.substring(0, tmi.length() - 1));
                                }
                                tmi.append("|");
                            }
                            else if (is_merged) {
                                if (undeterminedTMI.get(0).equals("MERGED")) {
                                    gcount++;
                                    tmi.append(entry.getKey()).append(",").append(ge.getGenome().get(currentChr).getChromosome().get(entry.getKey()).getType()).append(":").append(undeterminedTMI.get(0)).append("|");
                                }

                            }
                            else{
                                gcount++;
                                tmi.append(entry.getKey()).append(",").append(ge.getGenome().get(currentChr).getChromosome().get(entry.getKey()).getType()).append(":").append(undeterminedTMI.get(0)).append("|");
                            }

                        }
                        bamO.setGcount(gcount);

                        //remove last |
                        //stupid terrible bad hack
                        if(tmi.length() > 2) {
                            while (tmi.charAt(tmi.length() - 1) == '|') {
                                tmi = new StringBuilder(tmi.substring(0, tmi.length() - 1));
                            }
                            while (tmi.charAt(0) == '|'){
                                tmi = new StringBuilder(tmi.substring(1, tmi.length()));
                            }
                        }

//                            System.out.println(tmi);
//                        bamO.setTmi(tmi.toString());
                        //todo replace
                        tmi = null;


                    }

                }

                //------------------------------------------------------------------------------------------------------

//                fw.append(bamO.toString());

                //todo clear this
                bamO = null;
                continue;
            }

            //check what reads we can ignore
            if (bamOutput.canIgnore(sr)) {
                continue;
            }
            //add to lookup
            lookup.put(sr.getReadName(), sr);
            //if chromosome changes clear lookup
            currentChr = sr.getReferenceName();
        }

//        System.out.println("Finished Bam file");

        //hashmap Gene tmiscore
        double RPKMall = 0;
        double RPKMPCR0 = 0;
        double totalreads = 0;
        totalreads = tmicount.values().size();
        for(String s: tmicount.keySet()) {
            long genelen = 0;
            Gene temp = ge.getGenome().get(ge.findChromosome(s)).getChromosome().get(s);
            for(RegionVector rv : temp.getExons().values()) {
                genelen+= rv.getlength();
            }

            //1000000000 * (reads per gene / total reads) / gene length
            RPKMall = 1000000000 * ((double)tmicount.get(temp.getId()) / totalreads) / genelen;

            //get reads with pcrindex = 0
            //new hashmap where pcr = 0
            RPKMPCR0 = 1000000000 * ((double)pcrzero.get(temp.getId()) / totalreads ) / genelen;
//            System.out.println(RPKMall + " " +  RPKMPCR0);

            rpkm.write(RPKMall + "\n");
            rpkmpcrO.write(RPKMPCR0 + "\n");

        }


        //close
        fw.close();
        rpkm.close();
        rpkmpcrO.close();
    }

    private static RegionVector intronsinOverlap(Region overlap, RegionVector allintrons) {
        RegionVector intronsinOverlap = new RegionVector();
        //loop for all introns
        for (Region intron : allintrons.regions) {
            if (intron.getEnd() >= overlap.getStart() && intron.getStart() <= overlap.getStart()) {
                Region newr = new Region(overlap.getStart(), intron.getEnd());
                intronsinOverlap.addRegion(newr);
                continue;
            }
            if (overlap.getStart() <= intron.getStart() && intron.getStart() <= overlap.getEnd()) {
                if (intron.getEnd() <= overlap.getEnd()) {
                    intronsinOverlap.addRegion(intron);
                } else {
                    Region newr = new Region(intron.getStart(), overlap.getEnd());
                    intronsinOverlap.addRegion(newr);
                }
            }
        }
        return intronsinOverlap;
    }

    private static RegionVector exonsFromAB(List<AlignmentBlock> alignmentBlocks) {
        // iterate over alignment blocks, add exons to list
        RegionVector exons = new RegionVector();
        Region prevregion = null;

        for (AlignmentBlock ab : alignmentBlocks) {
            int refstart = ab.getReferenceStart();
            int readstart = ab.getReadStart();
            int refend = refstart + ab.getLength() - 1;
            int readend = readstart + ab.getLength();

            //check if == to refstart
            if (prevregion != null && prevregion.getEnd() + 1 == refstart) {
                prevregion.setEnd(refend);
                continue;
            }
            Region newexon = new Region(refstart, refend);
            exons.addRegion(newexon);
            exons.remove_overlap();
            prevregion = newexon;
        }
        return exons;
    }

//    //search through intervall tree and find a gene with given start and stop
//    public Gene findGene(Region searchregion, IntervalTree <Region> correcttree){
//
//        for(Region r: correcttree){
//            if(searchregion == r){
//                return r.getG();
//            }
//        }
//
//        System.out.println("No gene found");
//        return null;
//    }




}
