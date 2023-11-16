import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class DAGNode {

    private String id;
    private HashSet<String> parents = new HashSet<>();    //direct parents that can be found under that id in obo file
    private HashSet<String> allParents = new HashSet<>();      //all parents (including implicit parents)
    private HashSet<Gene> genes = new HashSet<>();
    private boolean enriched = false;
    private String name;
    private int size;   //for output size

    private Double hg_pval = 0.0;
    private Double fej_pval = 0.0;
    private Double ks_pval = 0.0;
    private Double ks_stat = 0.0;

    private HashSet<ArrayList<String>> paths = new HashSet<>(); //non corrected paths

    private HashSet<String> correct_paths = new HashSet<>(); //corrected paths

    private String shortest_path_to_a_true;


    public Double getKs_stat() {
        return ks_stat;
    }

    public void setKs_stat(Double ks_stat) {
        this.ks_stat = ks_stat;
    }

    public Double getHg_pval() {
        return hg_pval;
    }

    public void setHg_pval(Double hg_pval) {
        this.hg_pval = hg_pval;
    }

    public Double getFej_pval() {
        return fej_pval;
    }

    public void setFej_pval(Double fej_pval) {
        this.fej_pval = fej_pval;
    }

    public Double getKs_pval() {
        return ks_pval;
    }

    public void setKs_pval(Double ks_pval) {
        this.ks_pval = ks_pval;
    }

    public int getSize() {
        return size;
    } //for output size

    public void setSize(int size) { //for output size
        this.size = size;
    }

    public void setGenes(HashSet<Gene> genes) {
        this.genes = genes;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public boolean isEnriched() {
        return enriched;
    }

    public void setEnriched(boolean enriched) {
        this.enriched = enriched;
    }

    public DAGNode(String id) {
        this.id = id;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public HashSet<String> getParents() {
        return parents;
    }

    public void addDirectParent(String parent) {
        this.parents.add(parent);
        this.allParents.add(parent);
    }

    public void addParent(String parent) {
        this.allParents.add(parent);
    }

    public void setParents(HashSet<String> parents) {
        this.parents = parents;
    }

    public HashSet<String> getAllParents() {
        return allParents;
    }

    public void setAllParents(HashSet<String> allParents) {
        this.allParents = allParents;
    }


    public void addGene(Gene g) {
        this.genes.add(g);
    }

    public HashSet<Gene> getGenes() {
        return genes;
    }

    public boolean hasCorrectSize(int minsize, int maxsize) {

        if (size >= minsize && size <= maxsize) {
            return true;
        }

        return false;


    }

    public HashSet<ArrayList<String>> getPaths() {
        return paths;
    }

    public void setPaths(HashSet<ArrayList<String>> paths) {
        this.paths = paths;
    }


    public HashSet<String> getCorrect_paths() {
        return correct_paths;
    }



    public void setCorrect_paths(HashSet<String> correct_paths) {
        this.correct_paths = correct_paths;
    }


    public String getShortest_path_to_a_true() {
        return shortest_path_to_a_true;
    }

    public void setShortest_path_to_a_true(String shortest_path_to_a_true) {
        this.shortest_path_to_a_true = shortest_path_to_a_true;
    }

    //todo noverlap
    public int getnoverlap(){
        int count = 0;
        for (Gene g: this.genes){
            if(g.isSignif()){
                count++;
            }
        }
        return count;
    }


    public void addPath(ArrayList<String> path) {

        paths.add(path);

    }


    //____________________________________________________________________________________________
    //shortest path
    public int getShortestPath_to_partner(String partner, HashMap<String, DAGNode> all_DAGNodes) {


        HashSet<String> commonAncestors = getCommonAncestors(partner, all_DAGNodes);


        ArrayList<String>[] shortest_path = new ArrayList[2];


        int shortest_path_size = Integer.MAX_VALUE;



        /*if (shortest_path[0] != null) {
            shortest_path_size = shortest_path[0].size() + shortest_path[1].size();
        }*/


        //get all paths to common ancestors
        for (String commonAncestor : commonAncestors) {

            ArrayList<String> path_ca = getPathToCommonAncestor(all_DAGNodes.get(id), commonAncestor);
            ArrayList<String> path_ca_partner = getPathToCommonAncestor(all_DAGNodes.get(partner), commonAncestor);

            //check if current path is shorter
            int len_new = path_ca.size() + path_ca_partner.size();

            if (len_new < shortest_path_size) {
                shortest_path[0] = saveArrayList(path_ca);
                shortest_path[1] = saveArrayList(path_ca_partner);
                shortest_path_size = shortest_path[0].size() + shortest_path[1].size();
            }

        }


        return shortest_path_size - 2;


    }

    private HashSet<String> getCommonAncestors(String go_id, HashMap<String, DAGNode> all_DAGNodes) {


        //d.getAllParents().add(d.getId()); //wichtig: die eigene id zu allParents hinzufuegen (das wurde schon gemacht weil shortest path von oboReader schon aufgerufen wurde?!)
        //enriched_node.getAllParents().add(enriched_node.getId()); //wichtig: die eigene id zu allParents hinzufuegen (das wurde schon gemacht weil shortest path von oboReader schon aufgerufen wurde?!)

        HashSet<String> allParents_partner = all_DAGNodes.get(go_id).getAllParents();    //all parents of enriched node


        HashSet<String> common_ancestors = new HashSet<>(); //hashset with intersections = common ancestors
        for (String s : allParents) {
            common_ancestors.add(s);
        }

        common_ancestors.retainAll(allParents_partner); //save intersection in common_ancestors


        return common_ancestors;

    }


    private ArrayList<String> saveArrayList(ArrayList<String> my_list) {

        ArrayList<String> result = new ArrayList<>();

        for (String s : my_list) {
            result.add(s);
        }

        return result;

    }

    public ArrayList<String> getPathToCommonAncestor(DAGNode node, String commonAncestor) {

        HashSet<String> all_paths_of_node = node.getCorrect_paths(); //all paths of original node

        //String min_path = "";
        ArrayList<String> min_path = new ArrayList<>(); //todo stimmt init ueber isEmpty??
        ArrayList<String> building_path = new ArrayList<>();


        loop1: for (String path : all_paths_of_node) { //go over all paths of node until ancestor is found


            building_path.clear();
            String[] ids_of_path = path.split("\\|");


            for (String id : ids_of_path) {

                building_path.add(id);

                if (id.equals(commonAncestor)) {

                    if (building_path.size() < min_path.size() || min_path.isEmpty()) { //current building_path is shorter

                        min_path.clear();
                        for (String s : building_path) { //add current path to min_path
                            min_path.add(s);
                        }
                        continue loop1;


                    }

                }

            }


        }

        return min_path;

    }

    //____________________________________________________________________________________________



    public double[] getNumOverlapping(String partner, HashMap<String,DAGNode> all_DAGNodes) {

        double[] result = new double[2];

        double ov_percent_first;
        double ov_percent_second;
        double max_ov_percent;

        HashSet<Gene> overlapping = new HashSet<>();
        for (Gene g: genes) {   //fill overlapping hashset with all genes (for retainAll)
            overlapping.add(g);
        }

        HashSet<Gene> genes_partner = all_DAGNodes.get(partner).getGenes(); //get all genes of partner go id

        overlapping.retainAll(genes_partner); //via retainAll get intersection of genes

        //max_ov_percent
        ov_percent_first = (double) overlapping.size() / (double) genes.size();
        ov_percent_second = (double) overlapping.size() / (double) genes_partner.size();

        max_ov_percent = Math.max(ov_percent_first, ov_percent_second);

        //result
        result[0] = (double) overlapping.size();
        result[1] = max_ov_percent;


        return result;
    }
}
