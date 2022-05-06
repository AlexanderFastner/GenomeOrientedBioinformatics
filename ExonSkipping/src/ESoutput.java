import java.util.Objects;
import java.util.TreeSet;
public class ESoutput {

    //------------------------------------------------------------------------------------------------------------------
    //outputs in order
    String gene_id;
    String gene_symbol;
    String chr;
    String strand;
    int nprots;
    int ntrans;
    String SV_intron;
    String WT_introns;
    String SV_prots;
    String WT_prots;
    int min_skippped_exons;
    int max_skipped_exons;
    int min_skipped_bases;
    int max_skipped_bases;

    //------------------------------------------------------------------------------------------------------------------
    //setters
    public void setGene_id(String gene_id) {
        this.gene_id = gene_id;
    }

    public void setGene_symbol(String gene_symbol) {
        this.gene_symbol = gene_symbol;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public void setNprots(int nprots) {
        this.nprots = nprots;
    }

    public void setNtrans(int ntrans) {
        this.ntrans = ntrans;
    }

    public void setSV_intron(String SV_intron) {
        this.SV_intron = SV_intron;
    }

    public void setWT_introns(String WT_introns) {
        this.WT_introns = WT_introns;
    }

    public void setSV_prots(String SV_prots) {
        this.SV_prots = SV_prots;
    }

    public void setWT_prots(String WT_prots) {
        this.WT_prots = WT_prots;
    }

    public void setMin_skippped_exons(int min_skippped_exons) {
        this.min_skippped_exons = min_skippped_exons;
    }

    public void setMax_skipped_exons(int max_skipped_exons) {
        this.max_skipped_exons = max_skipped_exons;
    }

    public void setMin_skipped_bases(int min_skipped_bases) {
        this.min_skipped_bases = min_skipped_bases;
    }

    public void setMax_skipped_bases(int max_skipped_bases) {
        this.max_skipped_bases = max_skipped_bases;
    }

    public String getGene_id() {
        return gene_id;
    }

    public String getGene_symbol() {
        return gene_symbol;
    }

    public String getChr() {
        return chr;
    }

    public String getStrand() {
        return strand;
    }

    public int getNprots() {
        return nprots;
    }

    public int getNtrans() {
        return ntrans;
    }

    public String getSV_intron() {
        return SV_intron;
    }

    public String getWT_introns() {
        return WT_introns;
    }

    public String getSV_prots() {
        return SV_prots;
    }

    public String getWT_prots() {
        return WT_prots;
    }

    public int getMin_skippped_exons() {
        return min_skippped_exons;
    }

    public int getMax_skipped_exons() {
        return max_skipped_exons;
    }

    public int getMin_skipped_bases() {
        return min_skipped_bases;
    }

    public int getMax_skipped_bases() {
        return max_skipped_bases;
    }

    @Override
    public String toString() {
        return gene_id + "\t" + gene_symbol + "\t" + chr + "\t" + strand + "\t" + nprots + "\t" +
                ntrans + "\t" + SV_intron + "\t" + WT_introns + "\t" + SV_prots + "\t" + WT_prots + "\t" +
                min_skippped_exons + "\t" + max_skipped_exons + "\t" + min_skipped_bases + "\t" + max_skipped_bases;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ESoutput eSoutput = (ESoutput) o;
        return getNprots() == eSoutput.getNprots() && getNtrans() == eSoutput.getNtrans() && getMin_skippped_exons() == eSoutput.getMin_skippped_exons() && getMax_skipped_exons() == eSoutput.getMax_skipped_exons() && getMin_skipped_bases() == eSoutput.getMin_skipped_bases() && getMax_skipped_bases() == eSoutput.getMax_skipped_bases() && Objects.equals(getGene_id(), eSoutput.getGene_id()) && Objects.equals(getGene_symbol(), eSoutput.getGene_symbol()) && Objects.equals(getChr(), eSoutput.getChr()) && Objects.equals(getStrand(), eSoutput.getStrand()) && Objects.equals(getSV_intron(), eSoutput.getSV_intron()) && Objects.equals(getWT_introns(), eSoutput.getWT_introns()) && Objects.equals(getSV_prots(), eSoutput.getSV_prots()) && Objects.equals(getWT_prots(), eSoutput.getWT_prots());
    }

    @Override
    public int hashCode() {
        return Objects.hash(getGene_id(), getGene_symbol(), getChr(), getStrand(), getNprots(), getNtrans(), getSV_intron(), getWT_introns(), getSV_prots(), getWT_prots(), getMin_skippped_exons(), getMax_skipped_exons(), getMin_skipped_bases(), getMax_skipped_bases());
    }
}
