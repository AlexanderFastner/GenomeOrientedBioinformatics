import java.io.*;
import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

public class MappingReader {

    public HashMap<String, Gene> getGenes() {
        return genes;
    }

    //key=gene_name, value=Gene(ensemblID, gene_name, gos)
    private HashMap<String, Gene> genes = new HashMap<String, Gene>();


    public MappingReader(String path_to_mappingFile, HashMap<String, DAGNode> all_DAGNodes , String mappingtype) {


        if (mappingtype.equals("go")) {

            try {

                InputStream in = new GZIPInputStream(new FileInputStream(path_to_mappingFile));
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                br.readLine();
                while (br.readLine() != null) {
                    String line = br.readLine();
                    String[] line_split = line.split("\t"); //split line

                    if (line == null) {
                        break;
                    }

                    if (line.startsWith("!")) {
                        continue;
                    }

                    //geneid
                    String geneid = line_split[2];
                    //qualifier
                    String qualifier = line_split[3];
                    //GO
                    String GOid = line_split[4];

                    if (!qualifier.equals("")) {
                        continue;
                    }

                    if (geneid.equals("")) {
                        continue;
                    }

//                    System.out.println(geneid);
//                    System.out.println(qualifier);
//                    System.out.println(GOid);

                    //System.out.println(line);

                    //put gene names into dagnodes
                    //if same name then add Go ids

                    DAGNode d = new DAGNode(GOid);
                    //System.out.println(all_DAGNodes.size());
                    if(all_DAGNodes.containsKey(GOid)){
                        //then add new gene name to that GOid
                        HashSet<String> s = new HashSet<>();
                        s.add(GOid);
                        Gene g = new Gene(null, geneid, s);
                        all_DAGNodes.get(GOid).addGene(g);
                        //System.out.println(g);
                    }
                    else{
                        //make new Dagnode
                        HashSet<String> s = new HashSet<>();
                        s.add(GOid);
                        Gene g = new Gene(null, geneid, s);
                        d.addGene(g);
                        all_DAGNodes.put(GOid, d);
                    }



                    //fill the Mappingreader genes Hashmap
                    Gene g;
                    //if genes doesnt have that geneName make a new one
                    if (!genes.containsKey(geneid)) {
                        HashSet<String> s = new HashSet<>();
                        s.add(GOid);
                        g = new Gene(null, geneid, s);
                        genes.put(geneid, g);
                    } else {
                        //add GOid to go ids of gene
                        genes.get(geneid).getGo_ids_of_gene().add(GOid);
                    }

                    //System.out.println(genes.size());




                }

                in.close();
                br.close();

            } catch (IOException io) {
                System.out.println("problem with file");
            }
        }


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
                        for (String go_id: gos_hashset) {
                            g.getGo_ids_of_gene().add(go_id);
                        }
                    }

                        for (String go_id : gos_hashset) {
                            if (all_DAGNodes.containsKey(go_id)) {
                                all_DAGNodes.get(go_id).addGene(g); //add gene to DAGNode

                                //add gene to all parents of go_id
                                for (String parent_go : all_DAGNodes.get(go_id).getAllParents()) {
                                    if (all_DAGNodes.containsKey(parent_go)) {
                                        all_DAGNodes.get(parent_go).addGene(g);
                                    }
                                }

                            }

                        }

                        genes.put(gene_name, g);     //eigentlich ist diese hashmap auch unnoetig

                    }
                    header = false;

                }

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }


        }
    }






}
