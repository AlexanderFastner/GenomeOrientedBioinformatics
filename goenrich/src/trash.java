import java.util.HashSet;

public class trash {


    /*
    while (i2 <= i && i2 < ids_of_one_enriched_path.length) {
                                if (ids_of_original_path[i].equals(ids_of_one_enriched_path[i2])) { //diese id mit id von enriched vergleichen und den darauffolgenden Eebenen
                                    //todo match
                                    String result = "";
                                    int counter = 0;
                                    for (String string : ids_of_original_path) {

                                        String name_for_result = all_DAGNodes.get(string).getName(); //get name of dagnode

                                        result += "|" + name_for_result;

                                        if (counter == i) {
                                            result += "*";
                                            break;
                                        }
                                        counter++;
                                    }
                                    counter = 0;
                                    for (String string : ids_of_one_enriched_path) {

                                        String name_for_result = all_DAGNodes.get(string).getName(); //get name of dagnode

                                        result += "|" + name_for_result;

                                        if (counter == i2) {    //todo das muss wahrscheinlich ueber result +=, Stringbuilder?!
                                            break;
                                        }
                                        counter++;
                                    }

                                    result = result.substring(1);
                                    d.setShortest_path_to_a_true(result);
                                    continue loop1;

                                }
                                i2++;
                            }
     */


    /*public void get_shortest_path(HashSet<DAGNode> enriched_DAGNodes) {

        loop1:
        for (DAGNode d : all_DAGNodes.values()) { //go over all existing dagnodes

            HashSet<String> paths_of_original_node = d.getCorrect_paths(); //all paths of original node

            //todo fuer diese for loop ein i einfuegen damit erst breite dann Tiefe
            int i = 0;
            int max = 1;
            while (i < max) {
                loop2:
                for (String path_original_node : paths_of_original_node) { //go through all paths of original dagnode
                    //for (int j = 0; j<paths_of_original_node.size(); j++ ) { //go through all paths of original dagnode
                    //String path_original_node = paths_of_original_node[j];

                    String[] ids_of_original_path = path_original_node.split("\\|"); //get each go term (=id) of original path

                    //determine max to identify stop
                    if (ids_of_original_path.length > max) {
                        max = ids_of_original_path.length;
                    }
                    //check if current path length is <=i
                    if (ids_of_original_path.length <= i - 1) { //todo stimmt <?
                        continue;
                    }

                    //for (int i = 0; i < ids_of_original_path.length; i++) {  // go over each go term (=id) of path and check if it equals enriched node


                    if (all_DAGNodes.get(ids_of_original_path[i]).isEnriched()) { //schauen ob der original node selbst enriched ist

                        //todo match bei isEnriched wird nichts hingeschrieben!?
                        /*String result = "";
                        int counter = 0;
                        for (String string : ids_of_original_path) {

                            String name_for_result = all_DAGNodes.get(string).getName();

                            result += "|" + name_for_result;

                            if (counter == i) {
                                result += " * ";
                                break;
                            }
                            counter++;
                        }


                        result = result.substring(1);
                        d.setShortest_path_to_a_true(result);
                        continue loop1;
                    } else if (i == 0) { // wenn nein eine Ebene hoch bei original node
                        continue;
                    }


                    int i2 = 1; //zaehlt die Ebene von enriched node
                    while (i2 <= i) {

                        for (DAGNode enriched_dag : enriched_DAGNodes) {    //todo go over all enriched dagnodes Achtung: pro Ebene!!

                            HashSet<String> paths_of_enriched_node = enriched_dag.getCorrect_paths(); //all paths of enriched node
                            //todo erst eine ebene komplett durchgehen vor level up
                            for (String path_enriched_node : paths_of_enriched_node) {  //go through all paths of enriched node

                                String[] ids_of_one_enriched_path = path_enriched_node.split("\\|");   //get each go term (=id) of enriched path

                                if (ids_of_one_enriched_path.length <= i2) { //checken, dass der enriched path lang genug ist
                                    continue;
                                }

                                String compare_this_enriched_id = ids_of_one_enriched_path[i2]; // this is the enriched id that is now being compared to original

                                if (ids_of_original_path[i].equals(ids_of_one_enriched_path[i2])) { //found match
                                    //todo match
                                    String result = "";
                                    int counter = 0;
                                    for (String string : ids_of_original_path) {

                                        String name_for_result = all_DAGNodes.get(string).getName(); //get name of dagnode

                                        result += "|" + name_for_result;

                                        if (counter == i) {
                                            result += " * ";
                                            break;
                                        }
                                        counter++;
                                    }
                                    counter = 0;
                                    for (String string : ids_of_one_enriched_path) {

                                        String name_for_result = all_DAGNodes.get(string).getName(); //get name of dagnode

                                        result += "|" + name_for_result;

                                        if (counter == i2 - 1) {    //todo das muss wahrscheinlich ueber result +=, Stringbuilder?!
                                            break;
                                        }
                                        counter++;
                                    }

                                    result = result.substring(1);
                                    d.setShortest_path_to_a_true(result);
                                    continue loop1;
                                }


                            }

                        }


                        i2++;


                    }
                    //continue loop2;

                    //}


                }
                i++;
            }


        }


    }*/


}
