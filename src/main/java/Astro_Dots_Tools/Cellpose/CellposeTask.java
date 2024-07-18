package Astro_Dots_Tools.Cellpose;

import Astro_Dots_Tools.Cellpose.CellposeTaskSettings;
import java.util.ArrayList;

public class CellposeTask {
    
    protected CellposeTaskSettings settings;

    public void setSettings(CellposeTaskSettings settings) {
        this.settings = settings;
    }

    public void run() throws Exception {
        ArrayList<String> options = new ArrayList<>();

        options.add("--dir");
        options.add("" + settings.datasetDir);

        options.add("--pretrained_model");
        options.add("" + settings.model);

        options.add("--chan");
        options.add("" + settings.ch1);

        if (settings.ch2 > -1) {
            options.add("--chan2");
            options.add("" + settings.ch2);
        }

        options.add("--diameter");
        options.add("" + settings.diameter);

        options.add("--flow_threshold");
        options.add("" + settings.flow_threshold);
       

        System.out.println("Cellpose version is set to:" + settings.version);
        if (settings.version.equals("0.6") || settings.version.equals("2.0")) {
            options.add("--cellprob_threshold");
        } else if (settings.version.equals("0.7") || settings.version.equals("1.0")) {
            options.add("--mask_threshold"); // supposed to be new flag name for 0.7 and 1.0 but not anymore in 2.
        }
        options.add("" + settings.cellprob_threshold);

        if (!settings.version.equals("0.6")) {
            if (settings.stitch_threshold > -1) {
                options.add("--stitch_threshold");
                options.add("" + settings.stitch_threshold);
                settings.use3D(false); // has to be 2D!
            }
    
            if (settings.omni)
                options.add("--omni");

            if (settings.cluster)
                options.add("--cluster");
        }
        
        if (settings.use3D) 
            options.add("--do_3D");
        
        if (settings.useGpu) 
            options.add("--use_gpu");
                
        if (settings.useFastMode) 
            options.add("--fast_mode");
        
        if (settings.useResample && !settings.version.equals("1.0") && !settings.version.equals("2.0"))
            options.add("--resample");
        
        if (settings.verbose)
            options.add("--verbose");

        options.add("--save_tif");

        options.add("--no_npy");
        
        
        Cellpose.execute(options, settings, null);
    }
}
