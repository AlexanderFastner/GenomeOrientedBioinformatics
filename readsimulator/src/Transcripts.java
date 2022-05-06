import java.util.HashMap;
import java.util.TreeSet;

public class Transcripts {

    private String transcriptID;
    //treeset of protIDs
    TreeSet<String> protIds = new TreeSet<>();

    public Transcripts (String TranscriptID) {
        this.transcriptID = TranscriptID;

    }

    public String getTranscriptID() {
        return transcriptID;
    }

    public void setTranscriptID(String transcriptID) {
        this.transcriptID = transcriptID;
    }

    public TreeSet<String> getProtIds() {
        return protIds;
    }

    public void setProtIds(TreeSet<String> protIds) {
        this.protIds = protIds;
    }
}
