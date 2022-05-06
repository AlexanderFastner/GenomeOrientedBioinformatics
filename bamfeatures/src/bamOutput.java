import net.sf.samtools.SAMRecord;

public class bamOutput {

    int readid;
    int mm;
    int clippingSize;
    int nsplit;
    int gcount;
    String tmi;
    int gdist;
    boolean antisense;
    int pcrindex;
    boolean intergenic;


    //todo add decider to only output the highest classification of transcriptomic merged intronic


    @Override
    public String toString() {
        if(nsplit == -1){
            return readid + "\tsplit-inconsistent:true\n";
        }

        if(tmi != null){
            return readid +
                    "\tmm:" + mm +
                    "\tclipping:" + clippingSize +
                    "\tnsplit:" + nsplit +
                    "\tgcount:" + gcount +
                    "\t" + tmi +
                    "\tpcrindex: " + pcrindex + "\n";
        }
        if(intergenic) {
            return
                    readid +
                            "\tmm:" + mm +
                            "\tclipping:" + clippingSize +
                            "\tnsplit:" + nsplit +
                            "\tgcount:" + gcount +
                            "\tgdist:" + gdist +
                            "\tantisense:" + antisense +
                            "\tpcrindex: " + pcrindex + "\n";
        }


        else return null;
    }

    public boolean isIntergenic() {
        return intergenic;
    }

    public void setIntergenic(boolean intergenic) {
        this.intergenic = intergenic;
    }

    public bamOutput(int readid) {
        this.readid = readid;
    }

    public int getReadid() {
        return readid;
    }

    public void setReadid(int readid) {
        this.readid = readid;
    }

    public int getMm() {
        return mm;
    }

    public void setMm(int mm) {
        this.mm = mm;
    }

    public int getClippingSize() {
        return clippingSize;
    }

    public void setClippingSize(int clippingSize) {
        this.clippingSize = clippingSize;
    }

    public int getNsplit() {
        return nsplit;
    }

    public void setNsplit(int nsplit) {
        this.nsplit = nsplit;
    }

    public int getGcount() {
        return gcount;
    }

    public void setGcount(int gcount) {
        this.gcount = gcount;
    }

    public String getTmi() {
        return tmi;
    }

    public void setTmi(String tmi) {
        this.tmi = tmi;
    }

    public int getGdist() {
        return gdist;
    }

    public void setGdist(int gdist) {
        this.gdist = gdist;
    }

    public boolean isAntisense() {
        return antisense;
    }

    public void setAntisense(boolean antisense) {
        this.antisense = antisense;
    }

    public int getPcrindex() {
        return pcrindex;
    }

    public void setPcrindex(int pcrindex) {
        this.pcrindex = pcrindex;
    }


    public static boolean canIgnore(SAMRecord sr){
        boolean test = false;
        if(sr.getNotPrimaryAlignmentFlag()){
            test = true;
        }
        if(sr.getReadUnmappedFlag()){
            test = true;
        }
        if(sr.getMateUnmappedFlag()) {
            test = true;
        }
        if(sr.getReferenceName().equals(sr.getMateReferenceName())){
            if((sr.getReadNegativeStrandFlag() == sr.getMateNegativeStrandFlag())){
                test = true;
            }
        }

        return test;
    }
}
