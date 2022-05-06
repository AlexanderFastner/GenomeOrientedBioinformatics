import java.util.HashMap;

public class Gene {

    //has attributes
    private int nProt, nTrans;
    private String id, strand, type, name, chr;
    //has a hashmap of many RegionVectors
    //String (ProtID) RregionVector sorted by their ProtID
    private HashMap<String, RegionVector> RegionVectors = new HashMap();

    //ProtID RV
    HashMap<String, RegionVector> CDS = new HashMap<>();

    //hashmap of RegionVectors which have Treesets (all introns in this gene)
    private HashMap<String, RegionVector> allIntrons = new HashMap<>();
    private HashMap<String, Transcripts> allTranscripts = new HashMap<>();


    public Gene(String id, String strand, String type, String name){
        this.id = id;
        this.strand = strand;
        this.type = type;
        this.name = name;
    }


    public void addProtID_to_transcript (String transID, String protID){
        //temp transcript
        Transcripts currentTrans;
        try {
            currentTrans = allTranscripts.get(transID);
            currentTrans.getProtIds().add(protID);

        } catch (Exception e) {
            currentTrans = new Transcripts(transID);
            allTranscripts.put(transID, currentTrans);
            currentTrans.getProtIds().add(protID);
        }
    }

    public void addCDS(String transID, int start, int end) {
        //temp var
        Region r = new Region(start, end);
        RegionVector current_RV;

        try {
            current_RV = CDS.get(transID);
            current_RV.addRegion(r);

        } catch (Exception e) {
            current_RV = new RegionVector();
            current_RV.addRegion(r);
            CDS.put(transID, current_RV);
        }
    }

    public HashMap<String, RegionVector> getCDS() {
        return CDS;
    }

    public void setCDS(HashMap<String, RegionVector> CDS) {
        this.CDS = CDS;
    }

    public void addTranscript(String t, Transcripts transcripts){
        allTranscripts.put(t, transcripts);
    }

    public HashMap<String, Transcripts> getAllTranscripts() {
        return allTranscripts;
    }

    public void setAllTranscripts(HashMap<String, Transcripts> allTranscripts) {
        this.allTranscripts = allTranscripts;
    }

    public void addRegionVector (String protID, RegionVector r){
        RegionVectors.put(protID, r);
    }

    public HashMap<String, RegionVector> getRegionVectors() {
        return RegionVectors;
    }

    public void setRegionVectors(HashMap<String, RegionVector> regionVectors) {
        RegionVectors = regionVectors;
    }

    public void addAllIntrons(String protIDINTRON, RegionVector r){
        allIntrons.put(protIDINTRON, r);
    }

    public HashMap<String, RegionVector> getAllIntrons() {
        return allIntrons;
    }

    public void setAllIntrons(HashMap<String, RegionVector> allIntrons) {
        this.allIntrons = allIntrons;
    }

    public int getnProt() {
        return nProt;
    }

    public void setnProt(int nProt) {
        this.nProt = nProt;
    }

    public int getnTrans() {
        return nTrans;
    }

    public void setnTrans(int nTrans) {
        this.nTrans = nTrans;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }
}
