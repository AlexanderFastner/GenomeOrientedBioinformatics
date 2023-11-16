import augmentedTree.IntervalTree;

import java.io.*;
import java.util.Objects;
import java.util.stream.StreamSupport;

public class GTFParser {

    public int nTrans = 0;


    public GTFParser(String s, Genome ge) {
        try {
        BufferedReader myreader = new BufferedReader(new FileReader(s));
        String line = myreader.readLine();
        String[] splitLine;

        //variables
        //--------------------------------------------------------------------------------------------------------------

        int nProts = 0;
        int exon_start = 0;
        int exon_end = 0;
        String chr = "";
        String gene_id = "";
        String gene_name = "";
        String gene_biotype = "";
        String gene_strand = "";
        String transcript_ID = "";
        String exon_id = "";
        String [] annotations = {};
        String [] split_annotations = {};
        String oldGene = "";
        String oldChr = "";
        //--------------------------------------------------------------------------------------------------------------
        //loop
        while (line != null) {
            line = myreader.readLine();
            if (line == null) {
                break;
            }
            splitLine = line.split("\t");
            if (splitLine.length < 8) {
                continue;
            }

            chr = splitLine[0];
            gene_strand = splitLine[6];
            annotations = (splitLine[8].replace("\"", "").replace("'", "").replace(" gene_id", "gene_id").split("; "));
            for (int i = 0; i < annotations.length; i++) {

                split_annotations = annotations[i].split(" ");
                //System.out.println(split_annotations[0]);
                if (split_annotations[0].equals("gene_id")) {
                    gene_id = split_annotations[1];
                }
                if (split_annotations[0].equals("gene_name")) {
                    gene_name = split_annotations[1];
                }
                if (split_annotations[0].equals("gene_biotype")) {
                    gene_biotype = split_annotations[1];
                }
                if (split_annotations[0].equals("transcript_id")) {
                    transcript_ID = split_annotations[1];
                }
            }
            //----------------------------------------------------------------------------------------------------------
            if (splitLine[2].equals("exon")) {
                nProts++;
                exon_start = Integer.parseInt(splitLine[3]);
                exon_end = Integer.parseInt(splitLine[4]);

                //System.out.println(annotations[0]);

                //make region
                Region r = new Region(exon_start, exon_end);

                //fill into data structure

                //find chromosome
                if (ge.getGenome().get(chr) != null) {
                    //check if gene Exists
                    if(ge.getGenome().get(chr).getChromosome().get(gene_id) != null) {
                        //check if this RegionVector exists
                        if(ge.getGenome().get(chr).getChromosome().get(gene_id).getRegionVectors().get(transcript_ID) != null) {
                            //if yes add region here
                            ge.getGenome().get(chr).getChromosome().get(gene_id).getRegionVectors().get(transcript_ID).addRegion(r);
                            //add to exons
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addExon(transcript_ID, exon_start, exon_end);
                            //add transcript
                            Transcripts t = new Transcripts(transcript_ID);
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addTranscript(transcript_ID, t);
                        }
                        //else make new RegionVector
                        else {
                            RegionVector new_regionVector = new RegionVector();
                            //add region to regionVector
                            new_regionVector.addRegion(r);
                            //add to CDS
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addExon(transcript_ID, exon_start, exon_end);
                            //add to gene
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addRegionVector(transcript_ID, new_regionVector);
                            //add transcript
                            Transcripts t = new Transcripts(transcript_ID);
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addTranscript(transcript_ID, t);
                        }
                    }
                    //else make new gene
                    else {
                        Gene new_gene = new Gene(gene_id, gene_strand, gene_biotype, gene_name);
                        //make regionVector
                        RegionVector new_regionVector = new RegionVector();
                        //add region to regionVector
                        new_regionVector.addRegion(r);
                        //add regionVector to gene
                        new_gene.addRegionVector(transcript_ID, new_regionVector);
                        //add gene to chromosome
                        ge.getGenome().get(chr).addGene(gene_id, new_gene);
                        //add to CDS
                        ge.getGenome().get(chr).getChromosome().get(gene_id).addExon(transcript_ID, exon_start, exon_end);
                        //add transcript
                        Transcripts t = new Transcripts(transcript_ID);
                        ge.getGenome().get(chr).getChromosome().get(gene_id).addTranscript(transcript_ID, t);
                    }
                }
                //else make new chromosome

                else{
                    //make chromosome
                    //System.out.println(chr);
                    Chromosome new_chromosome = new Chromosome(chr);
                    //make gene
                    Gene new_gene = new Gene(gene_id, gene_strand, gene_biotype, gene_name);
                    //make regionVector
                    RegionVector new_regionVector = new RegionVector();
                    //add region to regionVector
                    new_regionVector.addRegion(r);
                    //add regionVector to gene
                    new_gene.addRegionVector(transcript_ID, new_regionVector);
                    //add gene to chromosome
                    new_chromosome.addGene(gene_id, new_gene);
                    //add chromosome to genome
                    ge.addChromosome(chr, new_chromosome);
                    //add to CDS
                    ge.getGenome().get(chr).getChromosome().get(gene_id).addExon(transcript_ID, exon_start, exon_end);
                    //add transcript
                    Transcripts t = new Transcripts(transcript_ID);
                    ge.getGenome().get(chr).getChromosome().get(gene_id).addTranscript(transcript_ID, t);
                }
                oldGene = gene_id;
            }
            //----------------------------------------------------------------------------------------------------------
            //ntrans
            else if(splitLine[2].equals("transcript") ) {
                //if no existing chr
                if(ge.getGenome().get(splitLine[0]) == null) {
                    //if no make chr
                    Chromosome new_chromosome = new Chromosome(chr);
                    //make gene
                    Gene new_gene = new Gene(gene_id, gene_strand, gene_biotype, gene_name);
                    //make trans
                    Transcripts t = new Transcripts(transcript_ID);
                    //add chr to genome
                    ge.addChromosome(chr, new_chromosome);
                    //add gene to chr
                    new_chromosome.addGene(gene_id, new_gene);
                    //add trans to gene
                    new_gene.getAllTranscripts().put(transcript_ID, t);
                }
                else {
                    //if yes check that chr for geneID
                    if (ge.getGenome().get(chr).getChromosome().get(gene_id) == null) {
                        //make new gene
                        Gene new_gene = new Gene(gene_id, gene_strand, gene_biotype, gene_name);
                        //add gene to chr
                        ge.getGenome().get(chr).addGene(gene_id, new_gene);
                        //make trans
                        Transcripts t = new Transcripts(transcript_ID);
                        //add trans
                        new_gene.getAllTranscripts().put(transcript_ID, t);
                    }
                    //gene already exists
                    else {
                        //then add transcript to this gene
                        //make trans
                        Transcripts t = new Transcripts(transcript_ID);
                        //add trans
                        ge.getGenome().get(chr).getChromosome().get(gene_id).getAllTranscripts().put(transcript_ID, t);
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        //todo build Intervaltrees
            //holds start end info on genes
        //for each chromosome
            //make Interval trees + and -
            //add (start, end) of each gene into tree

//        System.out.println(ge.getGenome().size());

        for (Chromosome c: ge.getGenome().values()){
//            System.out.println(c.getChrID() + " " + c.getChromosome().size());
            //fill IntervalTree
            for(Gene g: ge.getGenome().get(c.getChrID()).getChromosome().values()){
//                System.out.println("gene" + g.exons.keySet());
                //get all exons of each gene

                long minstart = Integer.MAX_VALUE;
                long maxstop = 0;
                for(RegionVector rv : g.getExons().values()) {
                    if (rv.getStart() < minstart) {
                        minstart = rv.getStart();
                    }
                    if (rv.getEnd() > maxstop) {
                        maxstop = rv.getEnd();
                    }
                }
                Region r = new Region(minstart, maxstop, g);
//                    System.out.println(g.getStrand());
                if(Objects.equals(g.getStrand(), "-")){
                    c.additN(r);
                }
                else{
                    c.additP(r);
                }


            }
        }


        }catch (FileNotFoundException e){
            System.out.println("No GTF found");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public int getnTrans() {
        return nTrans;
    }

    public void setnTrans(int nTrans) {
        this.nTrans = nTrans;
    }
}
