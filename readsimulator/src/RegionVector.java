import java.util.Comparator;
import java.util.TreeSet;

public class RegionVector {

    //A region vector has a tree set of regions

    TreeSet<Region> regions = new TreeSet<>();

    public RegionVector() {
    }

    public RegionVector generateIntrons () {
        long intronStart = 0;
        long intronEnd = 0;
        boolean first = true;
        RegionVector inversed = new RegionVector();

        //for every entry in regions
        for (Region r : regions) {
            //first run only get intron start
            if (!first){
                //next exon end as intron start
                intronEnd = r.getStart();
                //make new_region
                Region new_region = new Region(intronStart, intronEnd);
                inversed.addRegion(new_region);
                //remember exon start as intron end
                intronStart = r.getEnd() + 1;
            }
            else {
                intronStart = r.getEnd() + 1;
            }
            first = false;
        }
        return inversed;
    }

    public RegionVector get_SV() {
        //check intron start - WT
        //check intron end - WT
        //both - SV
        RegionVector SV = new RegionVector();

        //THIS
        //loop through to find sv_candidate
        //System.out.println(regions.size());
        for (Region sv_candidate : regions) {
            //System.out.println(sv_candidate.getStart());

            long start_sv_candidate = sv_candidate.getStart();
            long end_sv_candidate = sv_candidate.getEnd();

            //compare against all other introns
            boolean first_wt_intron = false;
            boolean sec_wt_intron = false;

            for (Region intron : regions) {
                //inner introns
                long start_intron = intron.getStart();
                long end_intron = intron.getEnd();
                //both
                if (start_intron == start_sv_candidate && end_intron == end_sv_candidate) {
                    //System.out.println("1");
                    continue;
                }
                //check if start is the same and end is higher
                if (start_intron == start_sv_candidate && end_intron < end_sv_candidate) {
                    first_wt_intron = true;
                    //System.out.println("2");
                }
                //check if end is same and start is higher
                if (end_intron == end_sv_candidate && start_intron > start_sv_candidate) {
                    //System.out.println("3");
                    sec_wt_intron = true;
                }

                //all conditions met means this is an SV intron
                if (first_wt_intron && sec_wt_intron) {
                    SV.addRegion(sv_candidate);
                    //System.out.println(sv_candidate.getStart());
                }
            }
        }

        //System.out.println("SV size" + SV.getRegions().size());
        return SV;
    }

    public boolean has_wt_introns(Region sv_intron) {
        //vars
        boolean start_wt_intron = false;
        boolean end_wt_intron = false;

        for (Region wt_intron : regions) {
            //check if start is the same and end is higher
            //if so change bool to true
            if (wt_intron.getStart() == sv_intron.getStart() && wt_intron.getEnd() < sv_intron.getEnd()) {
                start_wt_intron = true;
            }
            if (wt_intron.getEnd() == sv_intron.getEnd() && wt_intron.getStart() > sv_intron.getStart()) {
                end_wt_intron = true;
            }
        }
        //found SV intron
        if (start_wt_intron && end_wt_intron) {
            return true;
        }
        //didnt find SV intron
        return false;
    }

    //testing
    @Override
    public String toString() {

        String new_string = null;

        boolean first = true;
        //for every region make output (mind the first)
        for (Region reg : regions) {
            if (!first) {
                new_string = new_string + "|" + reg.getStart() + ":" + reg.getEnd();
            } else {
                new_string = reg.getStart() + ":" + reg.getEnd();
            }
            first = false;
        }
        return new_string;
    }


    public void addRegion (Region r){
        regions.add(r);
    }

    public TreeSet<Region> getRegions() {
        return regions;
    }


}
