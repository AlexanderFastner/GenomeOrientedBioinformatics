import java.util.Comparator;
public class Region implements Comparable <Region>{

    private int start;
    private int end;
    private int length;

    public Region (int Start, int End){
        this.start = Start;
        this.end = End;
        this.length = End-Start;
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

    public int getLength() {
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