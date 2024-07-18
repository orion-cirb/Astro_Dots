package Astro_Dots_Tools.Cellpose;

import Astro_Dots_Tools.Cellpose.CellposeTaskSettings;
import Astro_Dots_Tools.Cellpose.CellposeTask;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.process.ImageConverter;
import net.imagej.ImageJ;
import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class CellposeSegmentImgPlusAdvanced {
    
    CellposeTaskSettings settings;
    ImagePlus imp;
    ImagePlus cellpose_imp;
    
    public CellposeSegmentImgPlusAdvanced(CellposeTaskSettings settings, ImagePlus imp) {
        this.settings = settings;
        this.imp = imp;
    }
    
    public ImagePlus run() {
        
        Calibration cal = imp.getCalibration();

        // Create temp folder to temporarily save the current time-point of the imp
        String tempDir = IJ.getDirectory("Temp");
        File cellposeTempDir = new File(tempDir, "cellposeTemp");
        cellposeTempDir.mkdir();

        // When plugin crashes, image files can pile up in the temp folder, so we make sure to clear everything
        File[] contents = cellposeTempDir.listFiles();
        if (contents != null) {
            for (File f : contents) {
                f.delete();
            }
        }

        // Add temp folder to the settings
        settings.setDatasetDir(cellposeTempDir.toString());


        if (settings.use3D == true) {
            if (imp.getNSlices() == 1) 
                System.out.println("WARNING: Can't use 3D mode on 2D image. 2D mode will be used");
                settings.use3D(false);
        }

        // Settings are done, we can now process the imp with Cellpose
        CellposeTask cellposeTask = new CellposeTask();
        cellposeTask.setSettings(settings);
        
        try {
            // Can't process time-lapse directly, so we'll save one time-point after another
            int impFrames = imp.getNFrames();

            // We'll use lists to store paths of saved input, output masks and outlines
            List<File> t_imp_paths = new ArrayList<>();
            List<File> cellpose_masks_paths = new ArrayList<>();
            List<File> cellpose_outlines_paths = new ArrayList<>();

            for (int t_idx = 1; t_idx <= impFrames; t_idx++) {
                // Duplicate all channels and all z-slices for a defined time-point
                ImagePlus t_imp = new Duplicator().run(imp, 1, imp.getNChannels(), 1, imp.getNSlices(), t_idx, t_idx);
                // Save the current t_imp into the cellposeTempDir
                File t_imp_path = new File(cellposeTempDir, imp.getShortTitle() + "-t" + t_idx + ".tif");
                FileSaver fs = new FileSaver(t_imp);
                fs.saveAsTiff(t_imp_path.toString());
                System.out.println(t_imp_path.toString());
                // Add to list of paths to delete at the end of operations
                t_imp_paths.add(t_imp_path);

                // Prepare path of the cellpose mask output
                File cellpose_imp_path = new File(cellposeTempDir, imp.getShortTitle() + "-t" + t_idx + "_cp_masks" + ".tif");
                cellpose_masks_paths.add(cellpose_imp_path);
                // Cellpose also creates a txt file (probably to be used with a script to import ROI in imagej), so we'll delete it too
                // (to generate ROIs from the label image we can use https://github.com/BIOP/ijp-larome)
                File cellpose_outlines_path = new File(cellposeTempDir, imp.getShortTitle() + "-t" + t_idx + "_cp_outlines" + ".txt");
                cellpose_outlines_paths.add(cellpose_outlines_path);
            }

            // Run CellPose
            cellposeTask.run();

            // Open all the cellpose_mask and store each imp within an ArrayList
            ArrayList<ImagePlus> imps = new ArrayList<>(impFrames);
            for (int t_idx = 1; t_idx <= impFrames; t_idx++) {
                ImagePlus cellpose_t_imp = IJ.openImage(cellpose_masks_paths.get(t_idx - 1).toString());
                // make sure to make a 16-bit imp
                // (issue with time-lapse, first frame have less than 254 objects and latest have more)
                if (cellpose_t_imp.getBitDepth() != 16) {
                    if (cellpose_t_imp.getNSlices() > 1) {
                        new ImageConverter(cellpose_t_imp).convertToGray16();
                    } else {
                        cellpose_t_imp.setProcessor(cellpose_t_imp.getProcessor().convertToShort(false));
                    }
                }
                imps.add(cellpose_t_imp.duplicate());
            }
            // Convert the ArrayList to an imp
            // https://stackoverflow.com/questions/9572795/convert-list-to-array-in-java
            ImagePlus[] impsArray = imps.toArray(new ImagePlus[0]);
            cellpose_imp = Concatenator.run(impsArray);
            cellpose_imp.setCalibration(cal);
            cellpose_imp.setTitle(imp.getShortTitle() + "-cellpose");

            // Delete the created files and folder
            for (int t_idx = 1; t_idx <= impFrames; t_idx++) {
                t_imp_paths.get(t_idx - 1).delete();
                cellpose_masks_paths.get(t_idx - 1).delete();
                cellpose_outlines_paths.get(t_idx - 1).delete();
            }
            cellposeTempDir.delete();
       
        } catch (Exception e) {
            e.printStackTrace();
        }
        return cellpose_imp;
    }
}