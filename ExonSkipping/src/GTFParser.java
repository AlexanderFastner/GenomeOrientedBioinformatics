import java.io.*;
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
        String protein_ID = "";
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
                if (split_annotations[0].equals("protein_id")) {
                    protein_ID = split_annotations[1].replace(";", "");
                }
            }
            //----------------------------------------------------------------------------------------------------------
            if (splitLine[2].equals("CDS")) {
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
                            //add to CDS
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addCDS(transcript_ID, exon_start, exon_end);
                            //add transcript
                            Transcripts t = new Transcripts(transcript_ID);
                            //add protID
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addTranscript(transcript_ID, t);
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addProtID_to_transcript(transcript_ID, protein_ID);
                        }
                        //else make new RegionVector
                        else {
                            RegionVector new_regionVector = new RegionVector();
                            //add region to regionVector
                            new_regionVector.addRegion(r);
                            //add to CDS
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addCDS(transcript_ID, exon_start, exon_end);
                            //add to gene
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addRegionVector(transcript_ID, new_regionVector);
                            //add transcript
                            Transcripts t = new Transcripts(transcript_ID);
                            //add protID
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addTranscript(transcript_ID, t);
                            ge.getGenome().get(chr).getChromosome().get(gene_id).addProtID_to_transcript(transcript_ID, protein_ID);
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
                        ge.getGenome().get(chr).getChromosome().get(gene_id).addCDS(transcript_ID, exon_start, exon_end);
                        //add transcript
                        Transcripts t = new Transcripts(transcript_ID);
                        //add protID
                        ge.getGenome().get(chr).getChromosome().get(gene_id).addTranscript(transcript_ID, t);
                        ge.getGenome().get(chr).getChromosome().get(gene_id).addProtID_to_transcript(transcript_ID, protein_ID);
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
                    ge.getGenome().get(chr).getChromosome().get(gene_id).addCDS(transcript_ID, exon_start, exon_end);
                    //add transcript
                    Transcripts t = new Transcripts(transcript_ID);
                    //add protID
                    ge.getGenome().get(chr).getChromosome().get(gene_id).addTranscript(transcript_ID, t);
                    ge.getGenome().get(chr).getChromosome().get(gene_id).addProtID_to_transcript(transcript_ID, protein_ID);
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

        //old introns
        //generate introns
//        for (var c : ge.getGenome().entrySet()){
//            for (var g : c.getValue().getChromosome().entrySet()) {
//                for(var rv : g.getValue().getRegionVectors().entrySet()){
//                    RegionVector new_intron = new RegionVector(rv.getKey());
//                    new_intron.generateIntrons(rv.getValue());
//                    g.getValue().addAllIntrons(new_intron.getTranscriptID(), new_intron);
//                }
//            }
//        }


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
