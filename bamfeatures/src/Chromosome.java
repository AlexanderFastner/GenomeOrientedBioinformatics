import augmentedTree.IntervalTree;

import java.util.ArrayList;
import java.util.HashMap;

public class Chromosome {

    //each chromosome has a string identifier
    private String chrID;
    //Each chromosome has many Genes
    private HashMap<String, Gene> chromosome = new HashMap<>();

    //interval trees to keep start and end positions of genes
    private IntervalTree<Region> itP = new IntervalTree<>();
    private IntervalTree<Region> itN = new IntervalTree<>();

    public IntervalTree<Region> getitP() {
        return itP;
    }

    public void setitP(IntervalTree<Region> it) {
        this.itP = it;
    }

    public void additP(Region region){
        itP.add(region);
    }

    public IntervalTree<Region> getitN() {
        return itN;
    }

    public void setitN(IntervalTree<Region> it) {
        this.itN = it;
    }

    public void additN(Region region){
        itN.add(region);
    }

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
