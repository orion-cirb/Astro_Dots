import Astro_Dots_Tools.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.MetadataTools;
import loci.formats.meta.IMetadata;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;


/**
 *  Detect astrocyte nucleus and processes
 *  Detect RNA dots
 *  Compute RNA dots volume and intensity in nucleus, processes and outside astrocyte 
 *  @author ORION-CIRB
 */
public class Astro_Dots implements PlugIn {
    
    private Tools tools = new Tools();
    
    public void run(String arg) {
        try {
            if (!tools.checkInstalledModules()) {
                return;
            }
            
            String imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            
            // Find extension of first image in input folder
            String fileExt = tools.findImageType(new File(imageDir));
            // Find all images with corresponding extension in folder
            ArrayList<String> imageFiles = tools.findImages(imageDir, fileExt);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("Error", "No images found with " + fileExt + " extension");
                return;
            }
            
            // Instantiate metadata and reader
            IMetadata meta = MetadataTools.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findCalibration(meta);

            // Find channels name
            String[] channels = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Generate dialog box
            String[] channelsOrdered = tools.dialog(channels);
            if (channelsOrdered == null) {
                return;
            }
           
            // Create output folder
            String outputDir = imageDir + File.separator + "Results_" + new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            File outDir = new File(outputDir);
            if (!Files.exists(Paths.get(outputDir))) {
                outDir.mkdir();
            }
            
            // Write headers results files
            // Global results
            FileWriter fileWriter = new FileWriter(outputDir + "results.csv", false);
            BufferedWriter results = new BufferedWriter(fileWriter);
            results.write("Image name\tROI name\tROI volume (µm3)\tSoma volume (µm3)\tProcesses volume (µm3)\tDots bg intensity"
                    + "\tDots total volume in soma (µm3)\tDots bg-corr total intensity in soma"
                    + "\tDots total volume in processes (µm3)\tDots bg-corr total intensity in processes"
                    + "\tDots total volume outside astrocyte (µm3)\tDots bg-corr total intensity outside astrocyte\n");
            results.flush();
            
            for (String f : imageFiles) {
                reader.setId(f);
                String imgName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + imgName + " ---");
                
                // Load ROIs
                tools.print("- Loading ROIs -");
                List<Roi> rois = tools.loadRois(imageDir, imgName);
                if (rois == null) continue;
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setSplitChannels(true);
                options.setQuiet(true);
                
                // Analyze nuclei channel
                tools.print("- Analyzing nuclei channel -");
                int indexCh = ArrayUtils.indexOf(channels, channelsOrdered[0]);
                ImagePlus imgNuc = BF.openImagePlus(options)[indexCh];
                ImagePlus maskNuc = tools.cellposeDetection(imgNuc);
                
                // Analyze astrocytes channel
                tools.print("- Analyzing astrocytes channel -");
                indexCh = ArrayUtils.indexOf(channels, channelsOrdered[1]);
                ImagePlus imgAstro = BF.openImagePlus(options)[indexCh];
                ImagePlus maskAstro = tools.astroSegmentation(imgAstro);
                        
                // Analyze RNA dots channel
                tools.print("- Analyzing RNA dots channel -");
                indexCh = ArrayUtils.indexOf(channels, channelsOrdered[2]);
                ImagePlus imgDots = BF.openImagePlus(options)[indexCh];
                ImagePlus maskDots = tools.dotsSegmentation(imgDots);
                double bgDots = tools.computeBackgroundNoise(imgDots, maskDots);
                
                for (Roi roi: rois) {
                    tools.print("- Analyzing ROI " + roi.getName() + " -");
                    int zStart = (Integer.parseInt(roi.getProperty("zStart")) < 1)? 1 : Integer.parseInt(roi.getProperty("zStart"));
                    int zStop = (Integer.parseInt(roi.getProperty("zStop")) > imgNuc.getNSlices())? imgNuc.getNSlices(): Integer.parseInt(roi.getProperty("zStop"));
                    
                    ImagePlus imgNucCrop = tools.cropImage(imgNuc, roi, zStart, zStop);
                    ImagePlus imgAstroCrop = tools.cropImage(imgAstro, roi, zStart, zStop);
                    ImagePlus imgDotsCrop = tools.cropImage(imgDots, roi, zStart, zStop); 
                        
                    Objects3DIntPopulation popNuc = tools.getPopulationInRoi(maskNuc, roi, zStart, zStop, tools.minNucVol, tools.maxNucVol);
                    Object3DInt objNuc = tools.nucleusSelection(popNuc, imgNucCrop, imgAstroCrop, imgDotsCrop);
                    if(objNuc == null) {
                        tools.print("WARNING: No nucleus found in ROI " + roi.getName() + ", ROI not analyzed");
                        continue;
                    }
                    Object3DInt objSoma = tools.dilateObject(objNuc, imgNucCrop, tools.nucDilation);
                    
                    ImagePlus maskAstroCrop = tools.getImageInRoi(maskAstro, roi, zStart, zStop);
                    objSoma.drawObject(ImageHandler.wrap(maskAstroCrop), 0);
                    Object3DInt objAstro = new Object3DInt(ImageHandler.wrap(maskAstroCrop));
                    ImagePlus locThickAstro = tools.localThickness3D(maskAstroCrop);
                    
                    ImagePlus maskDotsCrop = tools.getImageInRoi(maskDots, roi, zStart, zStop);
                    List<Object3DInt> objsDots = tools.dotsClassification(maskDotsCrop, objSoma, objAstro);
                    
                    tools.writeResults(roi, objSoma, objAstro, objsDots, bgDots, imgDotsCrop, results, imgName);
                    tools.drawResults(objSoma, maskAstroCrop, objsDots, imgNucCrop, imgAstroCrop, imgDotsCrop, locThickAstro, outputDir, imgName, roi.getName());
                    
                    tools.closeImage(imgNucCrop);
                    tools.closeImage(imgAstroCrop);
                    tools.closeImage(imgDotsCrop);
                    tools.closeImage(maskAstroCrop);
                    tools.closeImage(locThickAstro);
                    tools.closeImage(maskDotsCrop);
                }

                tools.closeImage(imgNuc);
                tools.closeImage(maskNuc);
                tools.closeImage(imgAstro);
                tools.closeImage(maskAstro);
                tools.closeImage(imgDots);
                tools.closeImage(maskDots);
            }
            
            results.close();
            tools.print("--- All done! ---");
        } catch (DependencyException | ServiceException | IOException | FormatException ex) {
            Logger.getLogger(Astro_Dots.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
