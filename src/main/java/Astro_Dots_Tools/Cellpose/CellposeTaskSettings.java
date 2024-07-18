package Astro_Dots_Tools.Cellpose;


public class CellposeTaskSettings {
    
    // Values defined from https://cellpose.readthedocs.io/en/latest/api.html
    String datasetDir;
    String model;
    int ch1;
    int ch2 = -1;

    int diameter;
    double flow_threshold = 0.4;
    double cellprob_threshold = 0.0;
    double stitch_threshold = -1;
   
    String cellposeEnvDirectory;
    String envType = "conda";
    String version = "2.0";
    boolean use3D = false;
    boolean useGpu = false;
    boolean useFastMode = false;
    boolean useResample = false;
    boolean omni = false;
    boolean cluster = false;
    boolean verbose = false;
    
    
    public CellposeTaskSettings(String model, int ch1, int diameter, String cellposeEnvDirectory) {
        this.model = model;
        this.ch1 = ch1;
        this.diameter = diameter;
        this.cellposeEnvDirectory = cellposeEnvDirectory;
    }
    
    public CellposeTaskSettings setDatasetDir(String datasetDir) {
        this.datasetDir = datasetDir;
        return this;
    }

    public CellposeTaskSettings setModel(String model) {
        this.model = model;
        return this;
    }

    public CellposeTaskSettings setChannel1(int ch1) {
        this.ch1 = ch1;
        return this;
    }

    public CellposeTaskSettings setChannel2(int ch2) {
        this.ch2 = ch2;
        return this;
    }

    public CellposeTaskSettings setDiameter(int diameter) {
        this.diameter = diameter;
        return this;
    }

    public CellposeTaskSettings setFlowTh(double flow_threshold) {
        this.flow_threshold = flow_threshold;
        return this;
    }

    public CellposeTaskSettings setCellProbTh(double cellprob_threshold) {
        this.cellprob_threshold = cellprob_threshold;
        return this;
    }
    
    public CellposeTaskSettings setStitchThreshold(double stitch_threshold) {
        this.stitch_threshold = stitch_threshold;
        return this;
    }
        
    public CellposeTaskSettings setCellposeEnvDirectory(String cellposeEnvDirectory) {
        this.cellposeEnvDirectory = cellposeEnvDirectory;
        return this;
    }
    
    public CellposeTaskSettings setVersion(String version) {
        this.version = version;
        return this;
    }
    
    public CellposeTaskSettings setEnvType(String envType) {
        this.envType = envType;
        return this;
    }

    public CellposeTaskSettings use3D(boolean use3D) {
        this.use3D = use3D;
        return this;
    }

    public CellposeTaskSettings useGpu(boolean useGpu) {
        this.useGpu = useGpu;
        return this;
    }
    
    public CellposeTaskSettings useFastMode(boolean useFastMode) {
        this.useFastMode = useFastMode;
        return this;
    }
    
    public CellposeTaskSettings useResample(boolean useResample) {
        this.useResample = useResample;
        return this;
    }

    public CellposeTaskSettings setOmni(boolean omni) {
        this.omni = omni;
        return this;
    }

    public CellposeTaskSettings setCluster(boolean cluster) {
        this.cluster = cluster;
        return this;
    }
    
    public CellposeTaskSettings setVerbose(boolean verbose) {
        this.verbose = verbose;
        return this;
    }
       
}
