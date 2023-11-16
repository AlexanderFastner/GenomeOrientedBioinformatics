import java.util.ArrayList;
import java.util.Comparator;
import java.util.Objects;
import java.util.TreeSet;

import augmentedTree.Interval;
import augmentedTree.IntervalTree;
import net.sf.samtools.*;

public class RegionVector {

    //A region vector has a tree set of regions

    TreeSet<Region> regions = new TreeSet<>();

    public RegionVector() {
    }

    public double getlength(){
        double d = 0;
        for(Region r: regions){
            d += r.getLength();
        }
        return d;
    }

    //check if a regionVector has any overlap with a given region
    public ArrayList<Region> checkWhatOverlap(Region region){
        ArrayList<Region> affected = new ArrayList<>();
        for(Region r: regions){
            //if there is an overlap add that region to list
            if(region.getStart() >= r.getStart() && region.getStart() <= r.getStop()){
                affected.add(r);
                affected.add(region);
            }
            if(region.getStop() <= r.getStop() && region.getStop() >= r.getStart()){
                affected.add(r);
                affected.add(region);
            }
            if(region.getStart() == r.getStop() + 1 || region.getStop() == r.getStart() - 1){
                affected.add(r);
                affected.add(region);
            }
        }
        return affected;
    }

    public RegionVector remove_overlap() {
        RegionVector no_overlap = new RegionVector();
        Region previous_region = null;

        for (Region r : regions) {
            if (previous_region == null) {
                previous_region = r;
                no_overlap.addRegion(r);
                continue;
            }
            //same start
            if (previous_region.getStart() == r.getStart()){
                previous_region.setEnd(Math.max(r.getStop(), previous_region.getStop()));
                continue;
            }

            //regions directly after each other
            if (previous_region.getStop() + 1 == r.getStart()){
                previous_region.setEnd(r.getStop());
            }

            if (r.getStart() <= previous_region.getStop()) {   //start before end of previous region
                previous_region.setEnd(Math.max(r.getStop(), previous_region.getStop()));
                continue;
            }

            no_overlap.addRegion(r);
            previous_region = r;
        }
        return no_overlap;
    }



    public boolean subRV(RegionVector readrv) {

        RegionVector cuttranscript = cut(readrv.getStart(), readrv.getEnd());


        //compare cut transcript to read_rv
        if (cuttranscript.equals(readrv)) {
            return true;
        }
        return false;
    }

    //cut from x1_cutStart to x2_cutEnd
    public RegionVector cut(int x1_cutStart, int x2_cutEnd) {
        RegionVector cut_rv = new RegionVector();

        for (Region r : regions) {
            if (r.getStart() <= x1_cutStart && r.getStop() >= x1_cutStart) {
                //cut this region == first Region in cut_rv
                if (r.getStop() <= x2_cutEnd) {
                    Region first_in_cutRV = new Region(x1_cutStart, r.getStop());
                    cut_rv.addRegion(first_in_cutRV);
                }
                else {
                    Region foo = new Region(x1_cutStart, x2_cutEnd);
                    cut_rv.addRegion(foo);
                }
            }

            //add middle (non-cut) regions
            if (r.getStart() > x1_cutStart && r.getStart() <= x2_cutEnd) {
                if (r.getStop() <= x2_cutEnd) {
                    cut_rv.addRegion(r);
                }
                // get cutend
                else {
                    Region last_in_cutRV = new Region(r.getStart(), x2_cutEnd);
                    cut_rv.addRegion(last_in_cutRV);
                    break;
                }
            }
        }
        return cut_rv;
    }


    //for merged transcriptomic
    //"all regions of read_rv are within regions of tr_rv"
    //check if all regions from read_rv are contained completely in transcript_rv
    public boolean contained(RegionVector read_rv) {
        boolean is_within; //true if read is within transcript_region
        //iterate over all read regions
        for (Region read_region : read_rv.regions) {
            is_within = false;

            // iterate over all transcript regions, and check if one of them contains read_region
            for (Region transcript_region : regions) {
                if (transcript_region.getStart() <= read_region.getStart() && transcript_region.getStop() >= read_region.getStop()) {
                    is_within = true;
                }
            }

            if (!is_within) {
                return false;
            }
        }
        return true;
    }




    //end inklusive
    public RegionVector generateIntrons () {
        long intronStart = 0;
        boolean first = true;
        RegionVector inversed = new RegionVector();

        //for every entry in regions
        for (Region r : regions) {
            //first run only get intron start
            if (!first){
                //next exon end as intron start
                int startexon = r.getStart();
                int intronEnd = startexon - 1;
                //make new_region
                Region new_region = new Region(intronStart, intronEnd);
                inversed.addRegion(new_region);
                //remember exon start as intron end
                long endexon = r.getEnd();
                intronStart = endexon + 1;
            }
            else {
                long endexon = r.getEnd();
                intronStart = endexon + 1;
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

    //output
    @Override
    public String toString() {

        String new_string = "";

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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        RegionVector that = (RegionVector) o;
        return Objects.equals(regions, that.regions);
    }

    @Override
    public int hashCode() {
        return Objects.hash(regions);
    }

    public void addRegion (Region r){
        regions.add(r);
    }

    public TreeSet<Region> getRegions() {
        return regions;
    }

    public int getStart() {
        return regions.first().getStart();
    }

    public int getEnd() {
        return regions.last().getStop();
    }

}
