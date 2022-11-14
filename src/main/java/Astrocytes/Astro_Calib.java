package Astrocytes;



import static Tools.Astro_Dots_Tools.bg;
import static Tools.Astro_Dots_Tools.find_background;
import static Tools.Astro_Dots_Tools.find_dots;
import static Tools.Astro_Dots_Tools.flush_close;
import static Tools.Astro_Dots_Tools.median_filter;
import static Tools.Astro_Dots_Tools.objectsSizeFilter;
import static Tools.Astro_Dots_Tools.threshold;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.Thresholder;
import java.awt.Color;
import java.awt.Font;
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
import mcib3d.image3d.ImageHandler;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
        
        
public class Astro_Calib implements PlugIn {
    
    private boolean canceled;
    private String imageDir;
    private static String outDirResults;
    private static Calibration cal = new Calibration();
    // min and max dots filter volume
    private final double minDotsSize = 0.03;
    private final double maxDotsSize = 20;
    // Threshold method for dots
    private static String thMethod = "";
    
    
    /**
     * Dialog ask for channels order
     * @param channels
     * @param showCal
     * @param cal
     * @return ch
     */
    public static ArrayList dialog(String[] channels, boolean showCal, Calibration cal) {
        ArrayList ch = new ArrayList();
        String[] thMethods = new Thresholder().methods;
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("Astrocyte : ", channels, channels[2]);
        gd.addChoice("Dots : ", channels, channels[3]);
        gd.addMessage("Other cellular type", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("No Astrocyte : ", channels, channels[0]);
        gd.addChoice("Dots Threshold Method", thMethods, thMethods[15]);
        if (showCal) {
            gd.addNumericField("XY pixel size : ", cal.pixelWidth, 3);
            gd.addNumericField("Z pixel size : ", 1, 3);
        }
        gd.showDialog();
        ch.add(0, gd.getNextChoice());
        ch.add(1, gd.getNextChoice());
        ch.add(2, gd.getNextChoice());
        thMethod = thMethods[gd.getNextChoiceIndex()];
        if (showCal) {
            cal.pixelWidth = gd.getNextNumber();
            cal.pixelDepth = gd.getNextNumber();
        }
        
        if(gd.wasCanceled())
            ch = null;
        return(ch);
    }
    
    /**
    * draw population
    * dots no astro in red
    * all dotsin green
    * astro image in grey
    * @param dotsPop
    * @param dotsPopNoAstro
    * @param img
    */
    public static void tagsObjects(Objects3DPopulation dotsPop, Objects3DPopulation dotsPopNoAstro, ImagePlus img, String imageName) {  
        ImagePlus imAstro = img.duplicate();
        ImageHandler imgObjDots = ImageHandler.wrap(imAstro).createSameDimensions();
        dotsPop.draw(imgObjDots);
        ImageHandler imgObjDotsNotAstro = ImageHandler.wrap(imAstro).createSameDimensions();
        dotsPopNoAstro.draw(imgObjDotsNotAstro);
        // save image for objects population
        ImagePlus[] imgColors = {imgObjDotsNotAstro.getImagePlus(), imgObjDots.getImagePlus(), null,
            imAstro};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
        imgObjects.updateAndDraw();
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(outDirResults + imageName + "-Objects.tif"); 
        flush_close(imgObjects);
    }
    
    
    
    
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
            
            // Write headers for calibration parameters file
            FileWriter fileCalibResults = new FileWriter(outDirResults + "Calibration_results_" + inDir.getName() + ".xls", false);
            BufferedWriter outPutDotsResults = new BufferedWriter(fileCalibResults);
            outPutDotsResults.write("rootName\tAstroMeanInt\tAstroStdInt\tBgMeanInt\n");
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

                        // return the index for channels DAPI, Astro, Dots and calibration (No Astro)
                        ch = dialog(channels, showCal, cal);

                        if ("None".equals(ch.get(2))) {
                            IJ.showStatus("Channel for no astro is empty !!");
                            return;
                        }
                        cal.setUnit("microns");
                        System.out.println("x/y cal = " +cal.pixelWidth+", z cal = " + cal.pixelDepth);
                    }
                        
                        /** 
                         * open channels
                         */
                        // Astrocyte channel
                        String chName = "_w" + ArrayUtils.indexOf(channels, ch.get(0)) + ch.get(0);
                        String astroFile = inDir + File.separator + rootName + chName + "_cmle.tif";
                        System.out.println("Opening Astrocyte channel : " + astroFile);
                        ImagePlus imgAstroOrg = IJ.openImage(astroFile);
                        ImagePlus imgAstro = new Duplicator().run(imgAstroOrg, 3, imgAstroOrg.getNSlices());
                        median_filter(imgAstro, 0.5);
                        flush_close(imgAstroOrg);
                        
                        // Dots channel
                        chName = "_w" + ArrayUtils.indexOf(channels, ch.get(1)) + ch.get(1);
                        String dotsFile = inDir + File.separator + rootName + chName + "_cmle.tif";
                        System.out.println("Opening Dots channel : " + dotsFile);
                        ImagePlus imgDotsOrg = IJ.openImage(dotsFile);
                        ImagePlus imgDots = new Duplicator().run(imgDotsOrg, 3, imgDotsOrg.getNSlices());
                        flush_close(imgDotsOrg);
                        
                        // No astro channel
                        chName = "_w" + ArrayUtils.indexOf(channels, ch.get(2)) + ch.get(2);
                        String NoAstroFile = inDir + File.separator + rootName + chName + "_cmle.tif";
                        System.out.println("Opening No Astro channel : " + NoAstroFile);
                        ImagePlus imgNoAstroOrg = IJ.openImage(NoAstroFile);
                        ImagePlus imgNoAstro = new Duplicator().run(imgNoAstroOrg, 3, imgNoAstroOrg.getNSlices());
                        flush_close(imgNoAstroOrg);

                        /**
                         * Compute mask
                         * 
                         **/
                        
                        // astro image
                        // make substack
                        ImagePlus imgAstroMask = imgAstro.duplicate();
                        median_filter(imgAstroMask, 2);
                        threshold(imgAstroMask, "Li", false, false);
                        IJ.run(imgAstroMask, "Options...", "iterations=1 count=1 do=Open stack");
                        
                        // No astro image
                        median_filter(imgNoAstro, 2);
                        threshold(imgNoAstro, "Li", false, false);
                        IJ.run(imgNoAstro, "Options...", "iterations=1 count=1 do=Open stack");
                        
                        /**
                         * Save masks
                         */
                        FileSaver imgAstroMaskFile = new FileSaver(imgAstroMask);
                        imgAstroMaskFile.saveAsTiff(outDirResults + rootName + "_AstroMask.tif");
                        FileSaver imgNoAstroMaskFile = new FileSaver(imgNoAstro);
                        imgNoAstroMaskFile.saveAsTiff(outDirResults + rootName + "_NoAstroMask.tif");
                        
                        
                        /**
                         * Find dots population
                        **/
                        Objects3DPopulation dotsPop = find_dots(imgDots, null, thMethod);
                        System.out.println("Dots number = "+dotsPop.getNbObjects());
                        objectsSizeFilter(minDotsSize, maxDotsSize, dotsPop, imgDots, false);
                        System.out.println("After size filter dots number = "+dotsPop.getNbObjects());
                        
                        /**
                         * Find dots population colocalized with no astro image and not with astro image
                         * Measure mean intensity for all colocalized dots
                         */
                        Objects3DPopulation dotsNoAstroPop = new Objects3DPopulation();
                        double astroMeanInt = 0;
                        double noAstroMaskMeanInt = 0;
                        double astroMaskMeanInt = 0;
                        DescriptiveStatistics stats = new DescriptiveStatistics(); 
                        for (int dot = 0; dot < dotsPop.getNbObjects(); dot++) {
                            Object3D dotObj = dotsPop.getObject(dot);
                            astroMaskMeanInt = dotObj.getPixMeanValue(ImageHandler.wrap(imgAstroMask));
                            noAstroMaskMeanInt = dotObj.getPixMeanValue(ImageHandler.wrap(imgNoAstro));
                            astroMeanInt = dotObj.getPixMeanValue(ImageHandler.wrap(imgAstro));
                            if (noAstroMaskMeanInt != 0  && astroMaskMeanInt == 0) {
                                dotsNoAstroPop.addObject(dotObj);
                                stats.addValue(astroMeanInt);
                            }
                        } 
                        System.out.println("No astro dots number = "+dotsNoAstroPop.getNbObjects());
                        
                        find_background(imgAstro);        
                        
                        /**
                         * Save results
                         */    
                        outPutDotsResults.write(rootName+"\t"+stats.getMean()+"\t"+stats.getStandardDeviation()+"\t"+
                                bg+"\n");
                        outPutDotsResults.flush();
                        
                        // save dots image
                        tagsObjects(dotsPop, dotsNoAstroPop, imgAstro, rootName);
                        flush_close(imgAstro);
                        flush_close(imgAstroMask);
                        flush_close(imgDots);
                        flush_close(imgNoAstro);
                }
            } 
            outPutDotsResults.close();
            IJ.showStatus("Calibration done....");
        } catch (DependencyException | ServiceException | IOException | FormatException ex) {
            Logger.getLogger(Astro_Calib.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}

