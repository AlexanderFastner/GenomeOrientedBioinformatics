import java.util.HashMap;
import java.util.HashSet;

public class Gene {

    String gene_id;
    String gene_name;
    private HashSet<String> go_ids_of_gene;
    private boolean isenriched;
    double fc;
    boolean signif;

    public Gene(String gene_id, String gene_name, HashSet<String> go_ids_of_gene) {
        this.gene_id = gene_id;
        this.gene_name = gene_name;
        this.go_ids_of_gene = go_ids_of_gene;
    }


    @Override
    public String toString() {
        String result;

        result = gene_id + "\t" + gene_name + "\t";

        for (String s : go_ids_of_gene) {
            result += s + "|";
        }
        if (result.endsWith("|")) {
            result = result.substring(0, result.length() - 1);
        }

        result += " fc=" + fc + " s=" + signif;

        return result;
    }

    public String getGene_name() {
        return gene_name;
    }

    public void setGene_name(String gene_name) {
        this.gene_name = gene_name;
    }

    public boolean isIsenriched() {
        return isenriched;
    }

    public void setIsenriched(boolean isenriched) {
        this.isenriched = isenriched;
    }

    public HashSet<String> getGo_ids_of_gene() {
        return go_ids_of_gene;
    }

    public void setGo_ids_of_gene(HashSet<String> go_ids_of_gene) {
        this.go_ids_of_gene = go_ids_of_gene;
    }

    public double getFc() {
        return fc;
    }

    public boolean isSignif() {
        return signif;
    }

    public void setFc(double fc) {
        this.fc = fc;
    }

    public void setSignif(boolean signif) {
        this.signif = signif;
    }


    public HashSet<HashSet<String>> getAllValidGoPairs(int minsizeInt, int maxsizeInt, HashMap<String, DAGNode> all_DAGnodes) {

        removeInvalidGoIds(minsizeInt, maxsizeInt, all_DAGnodes); //only keep go ids with size in [minsize, maxsize]

        HashSet<HashSet<String>> all_validGoPairs = new HashSet<>();


        // get all pairs (via arrays)
        String[] ids_in_array = new String[go_ids_of_gene.size()];
        go_ids_of_gene.toArray(ids_in_array);

        //for (int i = 0; i < one.size() - 1; i++) {
        //            for (int j = 1; j < one.size(); j++) {

        for (int i = 0; i < ids_in_array.length; i++) {
            for (int j = i + 1; j < ids_in_array.length; j++) {

                HashSet<String> validPair = new HashSet<>();
                validPair.add(ids_in_array[i]);
                validPair.add(ids_in_array[j]);

                all_validGoPairs.add(validPair);


            }


        }


        return all_validGoPairs;
    }

    private void removeInvalidGoIds(int minsizeInt, int maxsizeInt, HashMap<String, DAGNode> all_DAGnodes) {

        HashSet<String> invalid_ids = new HashSet<>();

        for (String go_id : go_ids_of_gene) {
            if (!all_DAGnodes.containsKey(go_id)) { // check if DAGNode tree contains go id
                invalid_ids.add(go_id);
                continue;
            }

            if (!all_DAGnodes.get(go_id).hasCorrectSize(minsizeInt, maxsizeInt)) {
                invalid_ids.add(go_id);
            }
        }

        for (String s : invalid_ids) {
            go_ids_of_gene.remove(s);
        }

    }
}
