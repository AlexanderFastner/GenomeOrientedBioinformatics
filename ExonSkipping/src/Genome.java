import java.util.HashMap;

public class Genome {
    //a genome has many chromosomes

    //A Genome has a hashmap of all the chromosomes
    private HashMap<String, Chromosome> genome = new HashMap();

    public Genome () {

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
