public class Output {

    private String term;
    private String name;
    private int size;
    private boolean is_true;
    private int noverlap;
    private double hg_pval;
    private double hg_fdr = 0.0;
    private double fej_fdr = 0.0;
    private double fej_pval;
    private double ks_stat;
    private double ks_pval;
    private double ks_fdr = 0.0;    //todo macht es sinn hier default values anzugeben? sonst nullpointerexception (?)
    private String shortest_path_to_a_true;


    @Override
    public String toString() {
        return
                term + '\t' +
                name + '\t' +
                size + '\t' +
                is_true + '\t' +
                noverlap + '\t' +
                hg_pval + '\t' +
                hg_fdr + '\t' +
                fej_pval + '\t' +
                fej_fdr + '\t' +
                ks_stat + '\t' +
                ks_pval + '\t' +
                ks_fdr + '\t' +
                shortest_path_to_a_true + '\n';
    }

    public double getFej_fdr() {
        return fej_fdr;
    }

    public void setFej_fdr(double fej_fdr) {
        this.fej_fdr = fej_fdr;
    }

    public double getKs_fdr() {
        return ks_fdr;
    }

    public void setKs_fdr(double ks_fdr) {
        this.ks_fdr = ks_fdr;
    }

    public String getTerm() {
        return term;
    }

    public void setTerm(String term) {
        this.term = term;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public int getSize() {
        return size;
    }

    public void setSize(int size) {
        this.size = size;
    }

    public boolean isIs_true() {
        return is_true;
    }

    public void setIs_true(boolean is_true) {
        this.is_true = is_true;
    }

    public int getNoverlap() {
        return noverlap;
    }

    public void setNoverlap(int noverlap) {
        this.noverlap = noverlap;
    }

    public double getHg_pval() {
        return hg_pval;
    }

    public void setHg_pval(double hg_pval) {
        this.hg_pval = hg_pval;
    }

    public double getHg_fdr() {
        return hg_fdr;
    }

    public void setHg_fdr(double hg_fdr) {
        this.hg_fdr = hg_fdr;
    }

    public double getFej_pval() {
        return fej_pval;
    }

    public void setFej_pval(double fej_pval) {
        this.fej_pval = fej_pval;
    }

    public double getKs_stat() {
        return ks_stat;
    }

    public void setKs_stat(double ks_stat) {
        this.ks_stat = ks_stat;
    }

    public double getKs_pval() {
        return ks_pval;
    }

    public void setKs_pval(double ks_pval) {
        this.ks_pval = ks_pval;
    }

    public String getShortest_path_to_a_true() {
        return shortest_path_to_a_true;
    }

    public void setShortest_path_to_a_true(String shortest_path_to_a_true) {
        this.shortest_path_to_a_true = shortest_path_to_a_true;
    }
}
