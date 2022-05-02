package Astro_Dot;




import static Tools.Astro_Dot_Tools.*;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */



public class Astro_Dot implements PlugIn {
    
    // Image directory
    public static String imageDir;
    private final boolean canceled = false;
    // output directory
    public static String outDirResults = "";
    public static Calibration cal = new Calibration();
    // Nucleus min and max filter volume
    public static double minNucSize = 50;
    public static double maxNucSize = 800;
    // dots min and max filter volume
    public final double minDotsSize = 0.03;
    public final double maxDotsSize = 20;
    // max volume of dot cluster
    public static double maxMacroDotsSize = 2000;
    // mean volume of all dots without clusters
    public static double meanSEMDotsSize = 0;
    // result buffermaxDotsSize
    private BufferedWriter outPutResults;
    
    
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            String[] imageFile = inDir.list();
            if (imageFile == null) {
                return;
            }
            // create output folder
            outDirResults = inDir + File.separator+ "Out"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write headers for global results file
            FileWriter fileResults = new FileWriter(outDirResults + "AstroDot_" + inDir.getName() + ".xls", false);
            outPutResults = new BufferedWriter(fileResults);
            outPutResults.write("ImageName\tRoi Name\tbg Intensity\tAstrocyte Volume\tDensity dots in Astro"
                    + "\tPercentage of dots not in astro\tPercentage of dots in soma\tPercentage of dots in fine processes\tPercentage of dots in large processes"
                    + "\tDots mean intensity in Astro\tMean process diameter\tPercentage of dot in GFAP neg\n");
            outPutResults.flush();
            
            // Write headers for dots classification results file
            FileWriter fileDotsResults = new FileWriter(outDirResults + "AstroDotClassification_" + inDir.getName() + ".xls", false);
            BufferedWriter outPutDotsResults = new BufferedWriter(fileDotsResults);
            outPutDotsResults.write("rootName\tRoi Name\t% dots in soma\t% dots in large processes\t% dots in fine processes\t% dots out of astrocyte\n");
            outPutDotsResults.flush();
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            Arrays.sort(imageFile);
            int imageNum = 0;
            ArrayList<String> ch = new ArrayList();
            for (int i = 0; i < imageFile.length; i++) {
                // For all ics or nd files
                String fileExt = FilenameUtils.getExtension(imageFile[i]);
                if ( "ics".equals(fileExt) || "nd".equals(fileExt)) {
                    String imageName = imageDir + imageFile[i];
                    String rootName = FilenameUtils.getBaseName(imageFile[i]);
                    // Find ROI file
                    String roi_file = imageDir+rootName+".zip";
                    if (!new File(roi_file).exists()) {
                        IJ.showMessage("No ROI file found !\nRoi file should be named as <nd or ics filename>.zip");
                        return;
                       }
                    else {
                        reader.setId(imageName);
                        int series = 0;
                        reader.setSeries(series);
                        int sizeC = reader.getSizeC();
                        imageNum++;
                        boolean showCal = false;
                        String[] channels = new String[sizeC+1];
                        if (imageNum == 1) {
                            // Find channels name
                            String channelNames = "None/";
                            if ("ics".equals(fileExt)) {
                                channels[0] = "None";
                                for (int n = 0; n < sizeC; n++) {
                                    channelNames += meta.getChannelExcitationWavelength(series, n).value().toString().replace(".0", "/");
                                }
                            }
                            else {
                                channelNames += meta.getImageName(0).replace("CSU_", "");
                            }
                            channels = channelNames.split("/");
                            // Check calibration
                            cal.pixelWidth = meta.getPixelsPhysicalSizeX(series).value().doubleValue();
                            cal.pixelHeight = cal.pixelWidth;
                            // problem to read calibration with nd files
                            if ((meta.getPixelsPhysicalSizeZ(series) == null) || (cal.pixelWidth == 1))
                                showCal = true;
                            else
                                cal.pixelDepth = meta.getPixelsPhysicalSizeZ(series).value().doubleValue();
                            
                            // return the index for channels DAPI, Astro, Dots and ask for calibration if needed 
                            ch = dialog(channels, showCal);
                            
                            if (ch == null) {
                                IJ.showStatus("Plugin cancelled !!!");
                                return;
                            }
                            cal.setUnit("microns");
                            System.out.println("x/y cal = " +cal.pixelWidth+", z cal = " + cal.pixelDepth);
                        }
                        
                        // find rois
                        RoiManager rm = new RoiManager(false);
                        rm.runCommand("Open", roi_file);
                        int index = 0;
                        
                        ImporterOptions options = new ImporterOptions();
                        options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                        options.setId(imageName);
                        options.setSplitChannels(true);
                        reader.setSeries(0); 
                        options.setQuiet(true);
                        
                        /**
                         * Open channels
                         */
                        
                        // Nucleus channel 
                        int channelIndex;
                        ImagePlus imgNuc, imgAstro, imgDots;
                        // Huygens ICS files
                        if ("ics".equals(fileExt)) {
                            channelIndex = ArrayUtils.indexOf(channels, ch.get(0)) -1 ;
                            System.out.println("Opening Nucleus channel : "+ rootName + "("+ch.get(0)+")");
                            imgNuc = BF.openImagePlus(options)[channelIndex];
                            
                            // Astrocyte channel
                            channelIndex = ArrayUtils.indexOf(channels, ch.get(1)) - 1;
                            System.out.println("Opening Astrocyte channel : " + rootName + "("+ch.get(1)+")");
                            imgAstro = BF.openImagePlus(options)[channelIndex];

                            // Dots channel
                            channelIndex = ArrayUtils.indexOf(channels, ch.get(2)) - 1;
                            System.out.println("Opening Dots channel : " + rootName + "("+ch.get(2)+")");
                            imgDots = BF.openImagePlus(options)[channelIndex];
                        }
                        // Huygens cmle tif files
                        else {
                            // Nucleus channel
                            String chName = "_w" + ArrayUtils.indexOf(channels, ch.get(0)) + "CSU-" + ch.get(0) ;
                            String dapiFile = inDir + File.separator + rootName + chName + ".TIF";
                            System.out.println("Opening Nucleus channel : "+ dapiFile);
                            imgNuc = IJ.openImage(dapiFile);
                            imgNuc.setCalibration(cal);
                        
                            // Astrocyte channel
                            chName = "_w" + ArrayUtils.indexOf(channels, ch.get(1)) + "CSU-" + ch.get(1);
                            String astroFile = inDir + File.separator + rootName + chName + "_cmle.tif";
                            System.out.println("Opening Astrocyte channel : " + astroFile);
                            imgAstro = IJ.openImage(astroFile);
                            imgAstro.setCalibration(cal);

                            // Dots channel
                            chName = "_w" + ArrayUtils.indexOf(channels, ch.get(2)) + "CSU-" + ch.get(2);
                            String dotsFile = inDir + File.separator + rootName + chName + "_cmle.tif";
                            System.out.println("Opening Dots channel : " + dotsFile);
                            imgDots = IJ.openImage(dotsFile);
                            imgDots.setCalibration(cal);
                        }
                        
                        // find background
                        find_background(imgAstro);
                        
                        // for each roi open image and crop
                        for (int r = 0; r < rm.getCount(); r++) {
                            index++;                            
                            // Find in roi name the desired top and bottom stack 
                            // roi name should be roi_number-ztop-zbottom
                            String[] regExp = rm.getName(r).split("-");
                            if (Integer.parseInt(regExp[1]) < 1)
                                regExp[1] = "1";
                            if (Integer.parseInt(regExp[2]) > imgNuc.getNSlices())
                                regExp[2] = Integer.toString(imgNuc.getNSlices());
                            int zStart = Integer.parseInt(regExp[1]);
                            int zStop = Integer.parseInt(regExp[2]);
                            
                            rm.select(imgNuc,r);
                            imgNuc.updateAndDraw();
                            Roi roiAstro = imgNuc.getRoi();
                            
                            /** 
                             * Nucleus
                             */
                            // make substack
                            ImagePlus imgNucZCrop = new Duplicator().run(imgNuc, zStart, zStop);
                            // bug in roi manager recenter roi
                            imgNucZCrop.setRoi(roiAstro);
                            roiAstro.setLocation(0, 0);
                            imgNucZCrop.updateAndDraw();
                            roiAstro = imgNucZCrop.getRoi();
                            
                            // find nucleus population
                            Objects3DPopulation nucPop = segmentNucleus2(imgNucZCrop, roiAstro);
                            System.out.println("Roi = "+index+" of " + rm.getCount());
                            System.out.println("Nucleus number = "+nucPop.getNbObjects());
                            objectsSizeFilter(minNucSize, maxNucSize, nucPop, imgNucZCrop, true);
                            System.out.println("After size filter Nucleus number = "+nucPop.getNbObjects());
                            flush_close(imgNucZCrop);
                            
                            if (nucPop.getNbObjects() != 0) {      
                                // astro channel  
                                rm.select(imgAstro, r);
                                imgAstro.updateAndDraw();
                                // make substack
                                ImagePlus imgAstroZCrop = new Duplicator().run(imgAstro, zStart, zStop); 
                                median_filter(imgAstroZCrop, 0.5);
                                imgAstroZCrop.setTitle(rootName+"_Astro");
                                // dots channel
                                rm.select(imgDots, r);
                                imgDots.updateAndDraw();
                                ImagePlus imgDotsZCrop = new Duplicator().run(imgDots, zStart, zStop);
                                
                                Object3D nucAstro = null;
                                // Find nucleus with intensity in astrocyte channel
                                // and ask to choose
                                    nucAstro = nucleusSelect(imgAstroZCrop, imgDotsZCrop, nucPop);                                  
                                    
                                if (nucAstro != null) {
                                    // compute distance map image
                                    ImagePlus imgAstroZCropMap = localThickness3D(imgAstroZCrop);
                                    ImagePlus imgAstroZMap = imgAstroZCropMap.duplicate();
                                    IJ.run(imgAstroZMap, "Enhance Contrast...", "saturated=0.3");
                                    IJ.run(imgAstroZMap, "Fire", "");
                                    IJ.resetMinAndMax(imgAstroZMap);
                                    IJ.run(imgAstroZMap, "Calibration Bar...", "location=[Upper Left] fill=White label=Black number=5 decimal=2 font=12 zoom=1 overlay show");
                                    FileSaver astroMapFile = new FileSaver(imgAstroZMap);
                                    astroMapFile.saveAsTiff(outDirResults + rootName + "_Roi"+ (r+1) + "_Map.tif");
                                    imgAstroZMap.close();
                                    
                                    // Find dots population
                                    Objects3DPopulation dotsPop = find_dots(imgDotsZCrop, roiAstro, thMethod);
                                    System.out.println("Dots number = "+dotsPop.getNbObjects());
                                    // Duplicate dots Population to filter dots without "macro-dots" to calculate the mean dot size
                                    Objects3DPopulation dotsPopNotMacro = dotsPop;
                                    objectsSizeFilter(minDotsSize, maxDotsSize, dotsPopNotMacro, imgDotsZCrop, false);
                                    meanSEMDotsSize = find_mean_dots_volume(dotsPopNotMacro);
                                    System.out.println("Dots mean size volume + sem = " + meanSEMDotsSize);
                                    // min size filter including macro-dots
                                    objectsSizeFilter(minDotsSize, Integer.MAX_VALUE, dotsPop, imgDotsZCrop, false);
                                    System.out.println("After min size filter dots number = "+dotsPop.getNbObjects());


                                    // calculate parameters
                                    double noAstroDot = classify_dots(nucAstro, dotsPop, imgAstro, imgAstroZCrop, imgAstroZCropMap);
                                    nucPop.updateNamesAndValues();
                                    // draw objects
                                    tagsObjects(nucAstro, dotsPop, imgAstroZCrop, outDirResults, rootName, r);                               

                                    // write a global parameters table by image
                                    compute_Image_parameters(roiAstro, r, rm.getCount(), imgAstroZCrop, imgAstroZCropMap, nucAstro, dotsPop, outPutResults,
                                            rootName, noAstroDot);
                                    // write dots parameters
                                    computeDotsParams(dotsPop, nucAstro, roiAstro.getName(), outPutDotsResults, rootName);
                                    flush_close(imgAstroZCrop);
                                    flush_close(imgAstroZCropMap);
                                }
                                else 
                                   System.out.println("No nucleus found !"); 
                                
                            }
                            else 
                                System.out.println("No nucleus found !");
                        }
                        
                        flush_close(imgNuc);
                        flush_close(imgAstro);
                        flush_close(imgDots);
                    }
                }
            }
            outPutResults.close();
            outPutDotsResults.close();
        }
        catch (IOException | FormatException | DependencyException | ServiceException ex) {
            Logger.getLogger(Astro_Dot.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.showStatus("Process done ...");
    }
}
