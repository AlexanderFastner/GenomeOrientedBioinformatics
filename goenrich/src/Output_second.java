public class Output_second {


    private String term1;
    private String term2;
    private boolean is_relative;
    private int path_length;
    private int num_overlapping;
    private double max_ov_percent;


    @Override
    public String toString() {
        return term1 + "\t" +
                term2 + "\t" +
                is_relative + "\t" +
                path_length + "\t" +
                num_overlapping + "\t" +
                max_ov_percent + "\n";
    }


    public String getTerm1() {
        return term1;
    }

    public void setTerm1(String term1) {
        this.term1 = term1;
    }

    public String getTerm2() {
        return term2;
    }

    public void setTerm2(String term2) {
        this.term2 = term2;
    }

    public boolean isIs_relative() {
        return is_relative;
    }

    public void setIs_relative(boolean is_relative) {
        this.is_relative = is_relative;
    }

    public int getPath_length() {
        return path_length;
    }

    public void setPath_length(int path_length) {
        this.path_length = path_length;
    }

    public int getNum_overlapping() {
        return num_overlapping;
    }

    public void setNum_overlapping(int num_overlapping) {
        this.num_overlapping = num_overlapping;
    }

    public double getMax_ov_percent() {
        return max_ov_percent;
    }

    public void setMax_ov_percent(double max_ov_percent) {
        this.max_ov_percent = max_ov_percent;
    }

}