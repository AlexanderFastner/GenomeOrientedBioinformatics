public class cds {
    //like region but with prot id cause im lazy
    int start;
    int end;
    String protein_id;

    public cds(String protein_id, int start, int end) {
        this.start = start;
        this.end = end;
        this.protein_id = protein_id;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public String getProtein_id() {
        return protein_id;
    }

    public void setProtein_id(String protein_id) {
        this.protein_id = protein_id;
    }
}
