package Astrocytes;




import static Tools.Astro_Dots_Tools.*;
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
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
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



public class Astro_Dots implements PlugIn {
    
    // Image directory
    public static String imageDir;
    private final boolean canceled = false;
    // output directory
    public static String outDirResults = "";
    public static Calibration cal = new Calibration();
    // Nucleus min and max filter volume
    private final double minNucSize = 20;
    private final double maxNucSize = 500;
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
            FileWriter fileResults = new FileWriter(outDirResults + "AstroDots_" + inDir.getName() + ".xls", false);
            outPutResults = new BufferedWriter(fileResults);
            outPutResults.write("ImageName\tRoi Name\tbg Intensity\tAstrocyte Volume\tDensity dots in Astro"
                    + "\tPercentage of dots not in astro\tPercentage of dots in soma\tPercentage of dots in fine processes\tPercentage of dots in large processes"
                    + "\tDots mean intensity in Astro\tMean astro dot diameter\tAstro Mean Diameter\tAstro Std diameter\tAstro median diameter\tPercentage of dot in GFAP neg\n");
            outPutResults.flush();
            
            // Write headers for dots parameters results file
            FileWriter fileDotsResults = new FileWriter(outDirResults + "AstroDotsPop_" + inDir.getName() + ".xls", false);
            BufferedWriter outPutDotsResults = new BufferedWriter(fileDotsResults);
            outPutDotsResults.write("rootName\tRoi Name\t#dots\tDot type\tdot volume\tDot Mean Int.\tAstro diameter\tDistance to nucleus\n");
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
                // Find nd files
                if (imageFile[i].endsWith(".nd")) {
                    String imageName = inDir+ File.separator+imageFile[i];
                    String rootName = imageFile[i].replace(".nd", "");
                    // Find ROI file
                    String roi_file = imageDir+rootName+".zip";
                    if (!new File(roi_file).exists()) {
                        IJ.showStatus("No ROI file found !");
                        return;
                       }
                    else {
                        reader.setId(imageName);
                        int series = 0;
                        reader.setSeries(0);
                        imageNum++;
                        boolean showCal = false;
                        String channelNames = "None/"+meta.getImageName(0);
                        String[] channels = channelNames.replace("_", "-").split("/");
                        // Check calibration
                        if (imageNum == 1) {
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
                        
                        // Find in calibration template astromeanInt and astrobgInt
                        //astroMeanIntTh = findCalibration(templateFile, "AstroMeanInt");
                        //astroMeanBg = findCalibration(templateFile, "BgMeanInt");
                        
                        // find rois
                        RoiManager rm = new RoiManager(false);
                        rm.runCommand("Open", roi_file);
                        int index = 0;
                        
                        /**
                         * Open channels
                         */
                        
                        // Nucleus channel 
                        String chName = "_w" + ArrayUtils.indexOf(channels, ch.get(0)) + ch.get(0);
                        String dapiFile = inDir + File.separator + rootName + chName + ".TIF";
                        System.out.println("Opening Nucleus channel : "+ dapiFile);
                        ImagePlus imgNuc = IJ.openImage(dapiFile);
                        imgNuc.setCalibration(cal);
                        
                        // Astrocyte channel
                        chName = "_w" + ArrayUtils.indexOf(channels, ch.get(1)) + ch.get(1);
                        String astroFile = inDir + File.separator + rootName + chName + "_cmle.tif";
                        System.out.println("Opening Astrocyte channel : " + astroFile);
                        ImagePlus imgAstro = IJ.openImage(astroFile);
                        imgAstro.setCalibration(cal);
                        
                        // find background
                        find_background(imgAstro);
                        
                        // Dots channel
                        chName = "_w" + ArrayUtils.indexOf(channels, ch.get(2)) + ch.get(2);
                        String dotsFile = inDir + File.separator + rootName + chName + "_cmle.tif";
                        System.out.println("Opening Dots channel : " + dotsFile);
                        ImagePlus imgDots = IJ.openImage(dotsFile);
                        imgDots.setCalibration(cal);
                        
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
                            ImagePlus imgNucSeg = segmentNucleus(imgNucZCrop, roiAstro);
                            
                            // WaterShed slipt
                            ImagePlus imgNucSplit = watershedSplit(imgNucSeg, 10);
                            imgNucSplit.setCalibration(cal);
                            IJ.run(imgNucSplit, "Fill Holes", "stack");
                            
                            // find nucleus population
                            Objects3DPopulation nucPop = getPopFromImage(imgNucSplit);
                            System.out.println("Roi = "+index+" of " + rm.getCount());
                            System.out.println("Nucleus number = "+nucPop.getNbObjects());
                            objectsSizeFilter(minNucSize, maxNucSize, nucPop, imgNucSplit, true);
                            System.out.println("After size filter Nucleus number = "+nucPop.getNbObjects());
                            flush_close(imgNucSeg);
                            flush_close(imgNucSplit);
                            
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
                                // if more than one nucleus find nucleus with intensity in astrocyte channel
                                // and ask to choose
                                //if (nucPop.getNbObjects() > 1) 
                                    nucAstro = nucleusSelect(imgAstroZCrop, imgDotsZCrop, nucPop);                                  
                                //else
                                //    nucAstro = nucPop.getObject(0);
                                if (nucAstro != null) {
                                    // compute distance map image
                                    ImagePlus imgAstroZCropMap = localThickness3D(imgAstroZCrop);
                                    imgAstroZCropMap.setCalibration(cal);
                                    IJ.run(imgAstroZCropMap, "Fire", "");
                                    FileSaver astroMapFile = new FileSaver(imgAstroZCropMap);
                                    astroMapFile.saveAsTiff(outDirResults + rootName + "_Roi"+ (r+1) + "_Map.tif");

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
                                    // draw objects
                                    tagsObjects(nucAstro, dotsPop, imgAstroZCrop, outDirResults, rootName, r);                               

                                    // write a global parameters table by image
                                    compute_Image_parameters(roiAstro, r, rm.getCount(), imgAstroZCrop, imgAstroZCropMap, nucAstro, dotsPop, outPutResults,
                                            rootName, noAstroDot);
                                    // write dots parameters
                                    computeDotsParams(imgAstroZCropMap, imgAstroZCrop, dotsPop, nucAstro, roiAstro.getName(), outPutDotsResults, rootName);
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
            Logger.getLogger(Astro_Dots.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.showStatus("Process done ...");
    }
}
