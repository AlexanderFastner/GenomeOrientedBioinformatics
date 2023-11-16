import java.util.HashMap;

public class Chromosome {

    //each chromosome has a string identifier
    private String chrID;
    //Each chromosome has many Genes
    private HashMap<String, Gene> chromosome = new HashMap<>();

    public Chromosome (String chrID){
        this.chrID = chrID;
    }

    public void addGene(String geneID, Gene g){
        chromosome.put(geneID, g);
    }

    public HashMap<String, Gene> getChromosome() {
        return chromosome;
    }

    public void setChromosome(HashMap<String, Gene> chromosome) {
        this.chromosome = chromosome;
    }

    public String getChrID() {
        return chrID;
    }

    public void setChrID(String chrID) {
        this.chrID = chrID;
    }
}
