import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenomicSequenceExtractor {

    private RandomAccessFile raf;
    private File idx;


    //constructor
    //gets sequence from file using the index and access in raf
    public GenomicSequenceExtractor(File fasta, File idx) throws FileNotFoundException {
        //initialize raf file on fasta
        raf = new RandomAccessFile(fasta, "r");
        //save index file n variable
        this.idx = idx;
    }


    //method getSequence
    //finds where start is in file using index
    //todo this is broken
    public String getSequence(String chr, long start, long end){
        try {
            BufferedReader br = new BufferedReader(new FileReader(idx));
            String line = br.readLine();
            String[] splitLine;
            //iterate through index file
            while(line != null) {
                splitLine = line.split("\t");
                //found right chromosome
                if(splitLine[0].equals(chr)){

                    //seek to starting pos
                    //start chr + start exon -1
                    //account for extra bytes due to \n
                    //get start in fa file
                    int numNewLines = (int)(((start - 1) / Integer.parseInt(splitLine[3])) * ((Integer.parseInt(splitLine[4]) - Integer.parseInt(splitLine[3]))));
                    long fileStart = Long.parseLong(splitLine[2]) + (start - 1) + numNewLines;
                    raf.seek(fileStart);

                    //calculate num of lines for bytes
                    long numLines = (long)Math.floor(((end-start))/(Integer.parseInt(splitLine[3])));
                    long numBytes = (end - start + 1) + (numLines * ((Integer.parseInt(splitLine[4]) - Integer.parseInt(splitLine[3]))));

                    byte[] b = new byte[(int)numBytes];
                    //read all at once
                    raf.readFully(b);
                    // byte[] to string
                    String s = new String(b, StandardCharsets.UTF_8);
                    //strip out \n
                    s = s.replaceAll("\n", "").replaceAll("\r", "");
                    //vergleichen mit wielange der sein sollte + 1
                    //check length of s in case skipped 1

                    if ((end - start + 1) != s.length()) {
                        //System.out.println("used");
                        //add one extra byte
                        raf.seek(fileStart + numBytes);
                        byte[] newb = new byte[1];
                        raf.readFully(newb);
                        String newstring = new String(newb, StandardCharsets.UTF_8);
                        s += newstring;
                    }

                    return s;

                    //old code
//                    long counter = 0;
//                    long exonLength = (end - start);
//                    List<Byte> b = new ArrayList<Byte>();
//                    //loop through length many char, ignoring \n
//                    while(counter < exonLength) {
//                        if(raf.readChar() == '\n'){
//                            counter = counter -1;
//                            break;
//                        }
//                        b.add((byte)raf.readChar());
//                        counter++;
//                    }
//                    //arraylist to array
//                    byte[] result = new byte[b.size()];
//                    for(int i = 0; i < b.size(); i++) {
//                        result[i] = b.get(i);
//                    }
//                    // byte[] to string
//                    String s = new String(result, StandardCharsets.UTF_8);
//                    //System.out.println(s);
//                    return s.toString();
                }
                line = br.readLine();
                if (line == null) {
                    break;
                }
            }
        }
        catch (FileNotFoundException e) {
            System.out.println("index File not found");
            e.printStackTrace();
        }
        catch (IOException i) {
            System.out.println("io exception");
            i.printStackTrace();
        }
        System.out.println("something went wrong");
        return null;
    }


}
