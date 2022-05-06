import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class oboReader {
    private HashMap<String, DAGNode> all_DAGNodes;
    private int numSigGenes;

    public oboReader(String obo, String root) {

        boolean idfound = false;

        all_DAGNodes = new HashMap<>(); //all terms/DAGNodes that were already found in the obo file
        DAGNode dag = null;
        try {
            BufferedReader myreader = new BufferedReader(new FileReader(obo));
            String line = myreader.readLine();
            boolean first_entry = true;

            while (line != null) {
                line = myreader.readLine();

                if (line == null) {
                    break;
                }

                if (line.startsWith("id: GO:")) { //get go id

                    if (!first_entry && idfound) {
                        all_DAGNodes.put(dag.getId(), dag); //add (previous) DAGNode to all_DAGNodes
                    }

                    //dag = new DAGNode(line.substring(line.length() - 7));
                    dag = new DAGNode(line.split("\\s+")[1]);   //save new GO id as DAGNode

                    idfound = true;
                    first_entry = false;
                    continue;
                }

                //if id not found keep going until an id is found
                if (!idfound) {
                    continue;
                }

                //if root is wrong do not save dag
                if (line.startsWith("namespace") && !line.contains(root)) {
                    //System.out.println(line);
                    idfound = false;
                    continue;
                }

                //also ignore if obselete
                if (line.startsWith("is_obsolete: true")) {
                    //System.out.println(line);
                    idfound = false;
                    continue;
                }

                //add name to Dagnode
                if (line.startsWith("name: ")) {
                    //System.out.println(line.substring(6));
                    dag.setName(line.substring(6));
                }

                if (line.contains("is_a:") && idfound) {    //we found a direct parent
                    String parent = line.split("\\s+")[1];  //get GO id of parent, currently in format "GO:2001316"
                    dag.addDirectParent(parent);
                }

                //System.out.println(line);

            }


            //add last entry to all_DAGNodes
            all_DAGNodes.put(dag.getId(), dag); //add (previous) DAGNode to all_DAGNodes


        } catch (Exception e) {
            System.out.println(e);
        }

        //non direct parents and path to shortest true
        for (String d : all_DAGNodes.keySet()) {    //go through all dag nodes

            ArrayList<String> all_paths = getParents(d); //get non corrected paths

            StringBuilder prev_path = new StringBuilder();  //last path that was added (we need this to get entire path for multiple parents)
            StringBuilder building_path = new StringBuilder(); //this path is being built
            boolean found_end = false;
            for (String id : all_paths) { //go through not corrected paths

                if (id.equals(d)) { //same id as original node (=starting point)
                    building_path = new StringBuilder(d);   //found original node (=starting point)
                } else if (!found_end) {

                    building_path.append("|").append(id); //add new id to corrected path


                } else if (found_end) { //-> jetzt new original node (=starting point)

                    String p = prev_path.toString();
                    String[] p_split = p.split("\\|");

                    for (String pp : p_split) { //go through previous path until you find LCA (splitting point)

                        building_path = building_path.append(pp);
                        if (pp.equals(id)) {    //we found LCA (splitting point)
                            break;
                        } else {
                            building_path.append("|");
                        }

                    }

                }


                if (all_DAGNodes.get(id).getParents().isEmpty()) { //found end
                    found_end = true;
                    all_DAGNodes.get(d).getCorrect_paths().add(building_path.toString()); //add correct path to dagnode
                    prev_path = building_path;

                    building_path = new StringBuilder(); //empty building path
                } else {
                    found_end = false;
                }


            }


            //System.out.println(all_DAGNodes.get(d).getCorrect_paths());

            /*if (d.equals("GO:0046209")) {
                System.out.println(all_DAGNodes.get(d).getAllParents());
                break;
            }*/


        }

    }

    public HashMap<String, DAGNode> getAll_DAGNodes() {
        return all_DAGNodes;
    }

    //___________________________________________________________________________________________________________
    //non direct parents
    public ArrayList<String> getParents(String go_id) {
        HashSet<String> parents = new HashSet<>();
        ArrayList<String> result = new ArrayList<>();   //wir brauchen jetzt eine ArrayList statt HashSet, um uns die Reihenfolge bis zu einem end point zu merken
        if (all_DAGNodes.containsKey(go_id)) {
            parents = all_DAGNodes.get(go_id).getParents(); //get direct parents
        }
        if (!parents.isEmpty()) {
            for (String p : parents) {

                all_DAGNodes.get(go_id).getAllParents().add(p); //add parents of direct parents (= non direct parents)

                result.add(go_id);
                //result.add(p);

                //recursion: get non direct parents bis man zu einem Endknoten kommt
                for (String x : getParents(p)) {
                    all_DAGNodes.get(go_id).getAllParents().add(x);
                    result.add(x);
                }

            }

        } else { //we found root bzw. end point
            result.add(go_id);
        }

        all_DAGNodes.get(go_id).addPath(result);
        return result;
    }


    //___________________________________________________________________________________________________________


    public void size() {

        int sizecounter;
        //System.out.println(all_DAGNodes.size());
        for (DAGNode d : all_DAGNodes.values()) {
            //System.out.println(d.getGenes());
            sizecounter = 0;
            for (Gene g : d.getGenes()) {
                //System.out.println("f");
                if (g.isIsenriched()) {
                    //System.out.println("out");
                    sizecounter++;
                }
            }
            d.setSize(sizecounter); //save size in each DAGNode

            //System.out.println(d.getName() + " " + sizecounter);
        }
    }

    public int numAllGenes() {
        HashSet<Gene> geneSet = allGenesInTree();

        //count signif
        int counter = 0;
        for (Gene g : geneSet) {
            if (g.isSignif()) {
                counter++;
            }
        }
        numSigGenes = counter;

        return geneSet.size();
    }


    public HashSet<Gene> allGenesInTree() {
        HashSet<Gene> geneSet = new HashSet<>();
        for (DAGNode d : all_DAGNodes.values()) {
            for (Gene g : d.getGenes()) {
                geneSet.add(g);
            }
        }
        return geneSet;
    }


    public int getNumSigGenes() {
        return numSigGenes;
    }

    public void get_shortest_path(HashSet<DAGNode> enriched_DAGNodes, int minsize, int maxsize) {


        loop1:
        for (DAGNode d : all_DAGNodes.values()) {//go over all existing dagnodes

            if (!d.hasCorrectSize(minsize, maxsize)) {
                continue;
            }


            ArrayList<String>[] shortest_path = new ArrayList[2];
            /*ArrayList<String> foo = new ArrayList<>(); // avoid nullpointerexception
            foo.add("");
            shortest_path[0] = foo;
            shortest_path[1] = foo;*/

            if (d.isEnriched()) {   //der node selbst ist enriched -> nichts bei shortest path reinschreiben
                d.setShortest_path_to_a_true("");
                continue;
            }

            d.getAllParents().add(d.getId()); //wichtig: die eigene id zu allParents hinzufuegen
            HashSet<String> original_AllParents = d.getAllParents();    //all parents of original node

            for (DAGNode enriched_node : enriched_DAGNodes) { //compare to all enriched nodes

                enriched_node.getAllParents().add(enriched_node.getId()); //wichtig: die eigene id zu allParents hinzufuegen
                HashSet<String> enriched_AllParents = enriched_node.getAllParents();    //all parents of enriched node


                HashSet<String> common_ancestors = new HashSet<>(); //hashset with intersections = common ancestors
                for (String s : original_AllParents) {
                    common_ancestors.add(s);
                }

                common_ancestors.retainAll(enriched_AllParents); //save intersection in common_ancestors

                int shortest_path_size = Integer.MAX_VALUE;

                if (shortest_path[0] != null) {
                    shortest_path_size = shortest_path[0].size() + shortest_path[1].size();
                }


                //get all paths to common ancestors
                for (String commonAncestor : common_ancestors) {

                    ArrayList<String> path_ca_original = getPathToCommonAncestor(d, commonAncestor);
                    ArrayList<String> path_ca_enriched = getPathToCommonAncestor(enriched_node, commonAncestor);

                    //check if current path is shorter
                    int len_new = path_ca_original.size() + path_ca_enriched.size();

                    if (len_new < shortest_path_size) {
                        shortest_path[0] = saveArrayList(path_ca_original);
                        shortest_path[1] = saveArrayList(path_ca_enriched);
                        shortest_path_size = shortest_path[0].size() + shortest_path[1].size();
                    }

                }

            }
            //result ids
            //System.out.println("\tID " + d.getId());
            StringBuilder result = new StringBuilder();
            for (String s : shortest_path[0]) {
                result.append("|").append(s);
            }
            result.append(" * ");
            int i = shortest_path[1].size() - 2;
            while (i >= 0) {
                result.append("|").append(shortest_path[1].get(i));
                i--;
            }
            result = new StringBuilder(result.substring(1));
            //System.out.println(result);

            //result for out
            result = new StringBuilder();
            for (String s : shortest_path[0]) {
                //if (all_DAGNodes.containsKey(s)) {
                result.append("|").append(all_DAGNodes.get(s).getName());
                //}
            }
            result.append(" * ");
            i = shortest_path[1].size() - 2;
            while (i >= 0) {
                result.append("|").append(all_DAGNodes.get(shortest_path[1].get(i)).getName());
                i--;
            }
            result = new StringBuilder(result.substring(1));
            d.setShortest_path_to_a_true(result.toString());
            //System.out.println(result);


            /*if (d.getId().equals("GO:0051046")) {
                break;
            }*/

        }


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
        ArrayList<String> min_path = new ArrayList<>();
        ArrayList<String> building_path = new ArrayList<>();


        loop1:
        for (String path : all_paths_of_node) { //go over all paths of node until ancestor is found


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


}
