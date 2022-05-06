import java.util.HashMap;

public class Genome {
    //a genome has many chromosomes

    //A Genome has a hashmap of all the chromosomes
    private HashMap<String, Chromosome> genome = new HashMap();

    public Genome () {

    }

    public String findChromosome(String geneID){
        //given a gene ID find what chromosome its in
        //search each chromosome in genome to find one with gene ID
        for(String chr :genome.keySet()){
            if (genome.get(chr).getChromosome().get(geneID) !=null){
                return chr;
            }
        }
        System.out.println("error in find chromosome");
        return null;
    }

    public void addChromosome(String chrNumber, Chromosome chr){
        genome.put(chrNumber, chr);
    }

    public void setGenome(HashMap<String, Chromosome> genome) {
        this.genome = genome;
    }

    public HashMap<String, Chromosome> getGenome() {
        return genome;
    }
}
