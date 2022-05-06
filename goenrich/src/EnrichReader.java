import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

public class EnrichReader {

    private HashSet<DAGNode> enriched_nodes = new HashSet<>();


    public EnrichReader(String path_to_enrichfile, HashMap<String, DAGNode> all_DAGNodes, HashMap<String, Gene> genes) {

        try {
            BufferedReader myreader = new BufferedReader(new FileReader(path_to_enrichfile));
            String line = "";


            String gene_name;
            double fc;
            boolean signif;


            while (line != null) {
                line = myreader.readLine();

                if (line == null) {
                    break;
                }

                if (line.startsWith("#")) {     //we found an enriched GO id

                    if (all_DAGNodes.containsKey(line.substring(1))) {
                        enriched_nodes.add(all_DAGNodes.get(line.substring(1)));  //add enriched node to hashset todo vllt brauchen wir das doch nicht??
                        all_DAGNodes.get(line.substring(1)).setEnriched(true);//set DAGNode enriched=true
                    }

                } else if (!line.startsWith("id")) {    //save fc and signif for each gene
                    //System.out.println(line);

                    gene_name = line.split("\\s+")[0];
                    //System.out.println(gene_name);
                    fc = Double.valueOf(line.split("\\s+")[1]);
                    signif = Boolean.parseBoolean(line.split("\\s+")[2]);


                    if (genes.containsKey(gene_name)) {
                        genes.get(gene_name).setFc(fc);
                        genes.get(gene_name).setSignif(signif);
                        genes.get(gene_name).setIsenriched(true);
                        //System.out.println(genes.get(gene_name));
                    }


                }


            }


        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }


    public HashSet<DAGNode> getEnriched_nodes() {
        return enriched_nodes;
    }

    public void setEnriched_nodes(HashSet<DAGNode> enriched_nodes) {
        this.enriched_nodes = enriched_nodes;
    }
}

