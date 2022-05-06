import java.util.Comparator;
public class Region implements Comparable <Region>{

    private long start;
    private long end;
    private long length;

    public Region (long Start, long End){
        this.start = Start;
        this.end = End;
        this.length = End-Start+1;
    }

    public long getStart() {
        return start;
    }

    public void setStart(long start) {
        this.start = start;
    }

    public long getEnd() {
        return end;
    }

    public void setEnd(long end) {
        this.end = end;
    }

    public long getLength() {
        return length;
    }

    public void setLength(int length) {
        this.length = length;
    }
    @Override
    public String toString() {
        return start + ":" + end;
    }

    @Override
    public int compareTo(Region o) {
        if (this.start == o.start && this.end == o.end){
            return 0;
        }
        else if (this.start < o.start){
            return -1;
        }
        return 1;
    }
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        Region region = (Region) o;
        return (start == region.start && end == region.end);
    }

}