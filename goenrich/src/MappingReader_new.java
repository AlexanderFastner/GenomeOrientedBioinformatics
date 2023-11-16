import java.io.*;
import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

public class MappingReader_new {

    public HashMap<String, Gene> getGenes() {
        return genes;
    }

    //key=gene_name, value=Gene(ensemblID, gene_name, gos)
    private HashMap<String, Gene> genes = new HashMap<>();


    public MappingReader_new(String path_to_mappingFile, HashMap<String, DAGNode> all_DAGNodes, String mappingtype) {


        if (mappingtype.equals("go")) {

            try {
                InputStream in = new GZIPInputStream(new FileInputStream(path_to_mappingFile));
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String line = "";
                String qualifier = "";
                String go_id = "";
                String gene_name = "";

                while (line != null) {
                    line = br.readLine();

                    if (line == null) {
                        break;
                    }

                    if(line.contains("GO:1900371")){
                        System.out.println("found");
                    }


                    if (line.startsWith("!")) {
                        continue;
                    }

                    String[] line_split = line.split("\t");

                    gene_name = line_split[2];
                    qualifier = line_split[3];
                    go_id = line_split[4];

                    if (!qualifier.equals("") || gene_name.equals("")) {    //Use only the associations without any association qualifier modifier
                        continue;
                    }



                    // handle gene
                    if (!genes.containsKey(gene_name)) { //gene not in genes HashMap yet
                        HashSet<String> go_ids_of_gene = new HashSet<>();
                        go_ids_of_gene.add(go_id);
                        Gene g = new Gene("", gene_name, go_ids_of_gene);
                        genes.put(gene_name, g);
                    } else { // gene already in HashMap
                        genes.get(gene_name).getGo_ids_of_gene().add(go_id);
                    }


                }


            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }


            //add genes to DAGNode tree
            //System.out.println(genes.size());
            boolean gene_in_tree;
            HashSet<String> remove_these_genes = new HashSet<>();
            for (Gene g : genes.values()) {
                gene_in_tree = false;

                HashSet<String> associated_goIds = g.getGo_ids_of_gene();

                HashSet<String> add_these_ids = new HashSet<>();

                for (String go_id : associated_goIds) {
                    if (all_DAGNodes.containsKey(go_id)) {
                        all_DAGNodes.get(go_id).getGenes().add(g);

                        //add gene to all parents of go_id
                        for (String parent : all_DAGNodes.get(go_id).getAllParents()) {
                            if (all_DAGNodes.containsKey(parent)) {
                                all_DAGNodes.get(parent).addGene(g); //add gene to DAGNode

                                add_these_ids.add(parent);

                            }

                        }





                        gene_in_tree = true;
                    }


                }
                //add go ids of parents to gene
                for (String s : add_these_ids) {
                    g.getGo_ids_of_gene().add(s);
                }
                if (!gene_in_tree) {
                    remove_these_genes.add(g.gene_name);
                }

            }

            //remove all Genes from genes that do not appear in dag node tree
            for (String g : remove_these_genes) {
                genes.remove(g);
            }


        }

        //__________________________________________________________________________________

        else {
            try {
                File input_mappingFile = new File(path_to_mappingFile);
                Scanner reader_mappingFile = new Scanner(input_mappingFile);

                //HashMap<String, gene> all_genes = new HashMap<String, gene>();

                String gene_id;
                String gene_name;
                String[] all_gos;

                //read mappingInfo file line by line
                boolean header = true;
                while (reader_mappingFile.hasNextLine()) {
                    String line = reader_mappingFile.nextLine(); //get next line
                    if (!header) {

                    /*if (line.startsWith("#")) {
                        continue;
                    }*/

                        String[] line_split = line.split("\t"); //split line

                        gene_id = line_split[0];
                        gene_name = line_split[1];

                        if (gene_name.equals("")) { //es kann sein, dass es keinen gene_name gibt (die Zeile soll dann ignoriert werden)
                            continue;
                        }

                        all_gos = line_split[2].split("\\|"); //get all gos
                        HashSet<String> gos_hashset = new HashSet<>();
                        for (String s : all_gos) {   //write gos to HashSet
                            gos_hashset.add(s);
                        }

                        Gene g;
                        if (!genes.containsKey(gene_name)) { //gene has not been found yet --> new gene
                            g = new Gene(gene_id, gene_name, gos_hashset); //new Gene
                        } else {
                            g = genes.get(gene_name); //gene was found before, only add new go_ids
                            for (String go_id : gos_hashset) {
                                g.getGo_ids_of_gene().add(go_id);
                            }
                        }

                        HashSet<String> add_these_parentIds_to_gene = new HashSet<>();

                        for (String go_id : gos_hashset) {
                            if (all_DAGNodes.containsKey(go_id)) {
                                all_DAGNodes.get(go_id).addGene(g); //add gene to DAGNode

                                //add gene to all parents of go_id
                                for (String parent_go : all_DAGNodes.get(go_id).getAllParents()) {
                                    if (all_DAGNodes.containsKey(parent_go)) {
                                        all_DAGNodes.get(parent_go).addGene(g);


                                        add_these_parentIds_to_gene.add(parent_go); // die parents go ids muessen auch zu dem Gen dazugehoeren!
                                    }

                                }
                            }

                        }

                        for (String id: add_these_parentIds_to_gene) {// die parent go ids muessen auch zu dem Gen dazugehoeren!
                            g.getGo_ids_of_gene().add(id);
                        }
                        genes.put(gene_name, g);

                    }
                    header = false;

                }

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }


        }
    }


}
