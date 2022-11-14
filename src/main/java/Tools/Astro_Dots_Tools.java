package Tools;

import static Astrocytes.Astro_Dots.cal;
import static Astrocytes.Astro_Dots.meanSEMDotsSize;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.*;
import ij.plugin.Duplicator;
import ij.plugin.GaussianBlur3D;
import ij.plugin.ImageCalculator;
import ij.plugin.RGBStackMerge;
import ij.plugin.Thresholder;
import ij.plugin.ZProjector;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.RankFilters;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Vector;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.processing.FastFilters3D;
import mcib3d.image3d.regionGrowing.Watershed3D;
import mcib3d.utils.ThreadUtil;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import sc.fiji.localThickness.Clean_Up_Local_Thickness;
import sc.fiji.localThickness.EDT_S1D;
import sc.fiji.localThickness.Local_Thickness_Parallel;


 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author phm
 */



public class Astro_Dots_Tools {

// Threshold methods for dots    
public static String thMethod = "";
// Astrocyte channel background
public static double bg = 0;
// if pureArn exclude no dots
public static boolean pureArn = false;
// ratio astro int / dots int for nucleus selection
public static boolean ratioInt = false; 
// value to subtract to dots image to remove shadow dots
private static final double bgDots = 500;
// min diameter of thin astrocyte process
private static double thinDiam;

    
    /**
     * Dialog ask for channels order and template file for calibration
     * @param channels
     * @param showCal
     * @return channel names
     */
    public static ArrayList dialog(String[] channels, boolean showCal) {
        ArrayList ch = new ArrayList();
        String[] thMethods = new Thresholder().methods;
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("DAPI : ", channels, channels[1]);
        gd.addChoice("IF : ", channels, channels[2]);
        gd.addChoice("Dots : ", channels, channels[3]);
        //gd.addMessage("Calibration template", Font.getFont("Monospace"), Color.blue);
        //gd.addFileField("Template file : ", templateFile);
        gd.addChoice("Dots Threshold Method", thMethods, thMethods[15]);
        gd.addCheckbox(" Specific mRNA", pureArn);
        //gd.addCheckbox(" mRNA in nucleus", ratioInt);
        //gd.addNumericField("Macro dots max size :  : ", maxMacroDotsSize, 0);
        if (showCal) {
            gd.addNumericField("XY pixel size : ", cal.pixelWidth, 3);
            gd.addNumericField("Z pixel size : ", 1, 3);
        }
        gd.showDialog();
        ch.add(0, gd.getNextChoice());
        ch.add(1, gd.getNextChoice());
        ch.add(2, gd.getNextChoice());
        //templateFile = gd.getNextString();
        thMethod = thMethods[gd.getNextChoiceIndex()];
        pureArn = gd.getNextBoolean();
        //ratioInt = gd.getNextBoolean(); 
        //maxMacroDotsSize = gd.getNextNumber();
        if (showCal) {
            cal.pixelWidth = gd.getNextNumber();
            cal.pixelDepth = gd.getNextNumber();
        }
        if(gd.wasCanceled())
            ch = null;
        return(ch);
    }
    
    
    /**Gaussian filter 
     * 
     * @param img
     * @param size
     */ 
    public static void gs_filter(ImagePlus img, double size) {
        GaussianBlur gaussian = new GaussianBlur();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            gaussian.blurGaussian(img.getProcessor(), size, size, 0.02);
            img.updateAndDraw();
        }
    }
    
    /*Median filter 
     * 
     * @param img
     * @param size
     */ 
    public static void median_filter(ImagePlus img, double size) {
        RankFilters median = new RankFilters();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            median.rank(img.getProcessor(), size, RankFilters.MEDIAN);
            img.updateAndDraw();
        }
    }
    
    
    /**
     * Find Calibration intensity in template file
     * @param templateFile
     * @param type
     * @return mean
     * @throws java.io.IOException
     */
    
    public static double findCalibration(String templateFile, String type) throws IOException {
        ResultsTable rt = ResultsTable.open(templateFile);
        double[] meanInt = new double[rt.getCounter()];
        for (int n = 0; n < rt.getCounter(); n++) {
            meanInt[n] = rt.getValue(type, n);
        }
        DescriptiveStatistics stats = new DescriptiveStatistics(meanInt);
        return(stats.getMean());
    }
    
    
    /**
     * Display detected astrocyte nucleus in green, others in red
     * User can select nucleus with checkboxes
     * @param imgAstro
     * @param imgDots
     * @param nucPop
     * @return selected nucleus
     */
    
    public static Object3D nucleusSelect(ImagePlus imgAstro, ImagePlus imgDots, Objects3DPopulation nucPop) {
        int nucSelected = find_astroNuc(nucPop, imgAstro, imgDots);
        ImagePlus img = imgAstro.duplicate();
        ImageHandler imgObjNuc = ImageHandler.wrap(img).createSameDimensions();
        ImageHandler imgObjNucSelected = ImageHandler.wrap(img).createSameDimensions();
        String[] nucIndex = new String[nucPop.getNbObjects()];
        boolean[] values = new boolean[nucPop.getNbObjects()];
        for (int i = 0; i < nucPop.getNbObjects(); i++) {
            Object3D obj = nucPop.getObject(i);
            nucIndex[i] = Integer.toString(i);
            values[i] = false;
            if (i == nucSelected) {
                obj.draw(imgObjNucSelected);
                labelsObject(obj, imgObjNucSelected.getImagePlus(), i, 255, 20);
                values[i] = true;
            }
            else {
                obj.draw(imgObjNuc);
                labelsObject(obj, imgObjNuc.getImagePlus(), i, 255, 20);
            }
        }
        
        IJ.run(imgObjNuc.getImagePlus(), "Enhance Contrast", "saturated=0.35");
        IJ.run(imgObjNucSelected.getImagePlus(), "Enhance Contrast", "saturated=0.35");
        img.setSlice(img.getNSlices()/2);
        IJ.run(img, "Enhance Contrast", "saturated=0.35");
        img.updateAndDraw();
        ImagePlus[] imgColors = {imgObjNuc.getImagePlus(), imgObjNucSelected.getImagePlus(), null, img, null};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        int nucZCenter = (int)nucPop.getObject(nucSelected).getCenterZ();
        imgObjects.setZ(nucZCenter+1);
        imgObjects.updateAndDraw();
        imgObjects.show("Nucleus");
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Astrocyte nucleus");
        gd.addMessage("Choose astrocyte nucleus");
        gd.addCheckboxGroup(nucPop.getNbObjects(), 1, nucIndex, values);
        gd.showDialog();            
        Vector<?> checkboxes = gd.getCheckboxes();
        imgObjects.hide();
        flush_close(imgObjects);
        Object3DVoxels objVoxel = new Object3DVoxels(imgObjNuc);
        // if one object
        if (nucPop.getNbObjects() == 1) {
            objVoxel = (Object3DVoxels)nucPop.getObject(0);
        }
        else {
            ArrayList<Object3D> objSelected = new ArrayList<>();
            for (int n = 0; n < checkboxes.size(); n++) {
                Checkbox nucCheck = (Checkbox) checkboxes.get(n);
                if (nucCheck.getState()) {
                    objSelected.add(nucPop.getObject(n));
                    //System.out.println(n);
                }
            }
            // if more than one object selected merge objects.
            objVoxel.addVoxelsUnion(objSelected);
            objVoxel.setValue(1);
        }
        if(gd.wasCanceled())
            objVoxel = null;
        flush_close(img);
        imgObjNuc.closeImagePlus();
        imgObjNucSelected.closeImagePlus();
        return(objVoxel);
    }
    
    
    
    
    /**
     * Threshold images, fill holes
     * @param img
     * @param thMed
     * @param fill 
     * @param calculate 
     */
    public static void threshold(ImagePlus img, String thMed, boolean fill, boolean calculate) {
        //  Threshold and binarize
        img.setZ(find_max(img));
        img.updateAndDraw();
        IJ.setAutoThreshold(img, thMed + " dark");
        Prefs.blackBackground = false;
        String method;
        if (calculate)
            method = "method="+thMed+" background=Dark calculate";
        else
            method = "method="+thMed+" background=Dark";
        IJ.run(img, "Convert to Mask", method);
        if (fill) {
            IJ.run(img,"Fill Holes", "stack");
        }
    }
    
    /**
     * Get objects 3D population from image
     * @param img
     * @return 
     */
    public static Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
    
    /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public static void clearOutSide(ImagePlus img, Roi roi) {
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(roi);
            ip.setBackgroundValue(0);
            ip.setColor(0);
            ip.fillOutside(roi);
        }
        img.updateAndDraw();
    }
    
    
    
    /**
     * Find mean dot volume
     * return meanVol+SEM
     * @param dotsPop
     * @return mean Vol + SEM
     */
    public static double find_mean_dots_volume(Objects3DPopulation dotsPop) {
        DescriptiveStatistics dotVolStats = new DescriptiveStatistics();
        for (int n = 0; n < dotsPop.getNbObjects(); n++) {
            Object3D dotObj = dotsPop.getObject(n);
            dotVolStats.addValue(dotObj.getVolumeUnit());
        }
        double meanSemVol = dotVolStats.getMean() + dotVolStats.getStandardDeviation() / Math.sqrt(dotVolStats.getN());
        return(meanSemVol);
    }
                            
    /**
     * Size filter objects
     * remove touching border
     * 
     * @param min
     * @param max
     * @param objects
     * @param img
     * @param border
    */
    public static void objectsSizeFilter(double min, double max, Objects3DPopulation objects, ImagePlus img, boolean border) {
        ImageHandler imh = ImageInt.wrap(img).createSameDimensions();
        Object3D obj = null;
        boolean remove = false;
        if (objects.getNbObjects() > 0) {
            for (int n = 0; n < objects.getNbObjects(); n++) {
                remove = false;
                obj = objects.getObject(n);
                double vol = obj.getVolumeUnit();
                //System.out.println(vol);
                // remove if touching border
                if (border) {
                    if (obj.touchBorders(imh, false)) {
                        remove = true;
                    }
                }
                // Size filter remove small
                if ((vol < min) || (vol > max)) {
                    remove = true;
                }
                if (remove) {
                    objects.removeObject(n);
                    n--;
                }
            }
        }
    }
    
    /**
     * Do Z projection
     * @param img
     * @param projection parameter
     * @return 
     */
    private static ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    /**
    * Find background image intensity
    * Z project min intensity
    * read mean intensity
    * @param img 
    */
    public static void find_backgroundProj(ImagePlus img) {
      img.deleteRoi();
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      bg = imp.getStatistics().mean;
      System.out.println("Background =  " + bg);
      flush_close(imgProj);
    }
    
    
    /**
    * Find background image intensity
    * find mean intensity outside astrocyte mask in astrocyte image
    * don't take in account pixels with 0 value
    * return mean intensity
     * @param imgAstro
    */
    public static void find_background(ImagePlus imgAstro) {
        ImagePlus img = new Duplicator().run(imgAstro, 3,imgAstro.getNSlices());
        median_filter(img, 0.5);
        ImagePlus imgMask = img.duplicate();
        threshold(imgMask, "Li", false, true);
        IJ.run(imgMask, "Invert", "stack");
        ImageCalculator imgCal = new ImageCalculator();
        ImagePlus img1 = imgCal.run("Multiply create stack", img, imgMask);
        IJ.run(img1, "Divide...", "value=255 stack");
        double sectionInt = 0;
        for (int n = 1; n <= img1.getNSlices(); n++) {
            double pixelInt = 0;
            int voxelZero = 0;
            img1.setZ(n);
            ImageProcessor ip = img1.getProcessor();
            for (int x = 0; x < img1.getHeight(); x++) {
                for (int y = 0; y < img1.getWidth(); y++) {
                    double voxelInt = ip.getPixelValue(x, y);
                    if ( voxelInt > 0) {
                        pixelInt += voxelInt;
                    }
                    else {
                        voxelZero++;
                    }
                }
            }
            sectionInt += pixelInt/((img1.getHeight()*img1.getWidth()) - voxelZero);
        } 
        bg = sectionInt/img1.getNSlices();
        System.out.println("Background = "+bg);
        flush_close(img);
        flush_close(imgMask);
        flush_close(img1);
    }
    
   /**
    * Find statistics diameter in distance map image
    */
    
    public static DescriptiveStatistics findDiameterStats(ImagePlus imgMap) {
       DescriptiveStatistics diamStats = new DescriptiveStatistics();
       for (int n = 1; n <= imgMap.getNSlices(); n++) {
            imgMap.setZ(n);
            ImageProcessor ip = imgMap.getProcessor();
            for (int x = 0; x < imgMap.getHeight(); x++) {
                for (int y = 0; y < imgMap.getWidth(); y++) {
                    double voxelInt = ip.getPixelValue(x, y);
                    if ( voxelInt > 0)
                        diamStats.addValue(voxelInt*cal.pixelWidth);
                }
            }
        } 
        return(diamStats);
    }
    
    /**
    * Find if a dot is on a large process 
    * return true if the max intensity voxel of the dot in DM image have a 8 neighour diameter > thinDiam
    * @param imgMap
     * @return 
    */
    public static boolean largeProcess(ImagePlus imgMap, Object3D dot) {
        ImageStack imgMapStack = imgMap.getImageStack();
        boolean diamOk = false;
        // thin process diameter in pixel
        // image map is in pixel
        thinDiam = cal.pixelDepth / cal.pixelWidth;
        // Voxel with max intensity in image map
        double maxValue = dot.getPixMaxValue(ImageHandler.wrap(imgMap));
        Voxel3D VolxelMax = dot.getPixelMax(ImageHandler.wrap(imgMap));
        // look if one of all 8 neighours > thinDiam
        if (maxValue <= thinDiam) {
            diamOk = neighboursOK((int)VolxelMax.x, (int)VolxelMax.y, (int)VolxelMax.z, maxValue, imgMapStack);
        }
        else
            diamOk = true;
        return(diamOk);
    }
    
    /***
     * Find 8 neighours
     * @param x
     * @param y
     * @param z
     * @param value
     * @param imgStackMap
     * @return true if on neighbourg have an int > value
     */
    private static boolean neighboursOK(int x, int y, int z, double value, ImageStack imgStackMap) {
        boolean neighboursOK = false;
        for (int slice = z - 1; slice <= (z + 1); slice++) {
            if ((!neighboursOK) && (slice < imgStackMap.getSize()) && (slice >= 0)) {
                for (int col = y - 1 ; col <= (y + 1) ; col++) {
                    if ((!neighboursOK) && (col < imgStackMap.getHeight()) && (col >= 0)){
                        for (int row = x - 1 ; row <= (x + 1) ; row++) {
                            if(!((col == y) && (row == x)) && (row < imgStackMap.getWidth()) && (row >= 0)) {
                                //System.out.println(row+","+col+","+slice);
                                double intMap = imgStackMap.getVoxel(row, col, slice);
                                // check if intensity > minDiam
                                if (intMap > value) {
                                   neighboursOK = true;
                                   break;
                                }
                            }
                        }
                    }
                }
            }
        }
        return(neighboursOK);
    }
    
    /**
     * Find roi volume
     * @param roi
     * @param img
     * @return volume
     */
    private static double roi_volume(Roi roi, ImagePlus imgAstro) {
        ImagePlus img = new Duplicator().run(imgAstro, 1, 1);
        ImageProcessor ip = img.getProcessor();
        ip.setRoi(roi);
        ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.AREA, cal);
        double volume = stats.area * cal.pixelDepth * imgAstro.getNSlices();
        flush_close(img);
        return(volume);
    }
    
    
    /**
     * Nucleus segmentation
     * @param imgNucZCrop
     * @param roiAstro
     * @return image mask
     */
    public static ImagePlus segmentNucleus(ImagePlus imgNucZCrop, Roi roiAstro) {
        imgNucZCrop.deleteRoi();
        median_filter(imgNucZCrop, 4);
        IJ.run(imgNucZCrop, "Difference of Gaussians", "  sigma1=30 sigma2=28 stack");
        threshold(imgNucZCrop, "Otsu", true, false);
        IJ.run(imgNucZCrop, "Options...", "iterations=5 count=1 do=Open stack");
        imgNucZCrop.setRoi(roiAstro);
        imgNucZCrop.updateAndDraw();
        IJ.run("Colors...", "foreground=black background=white selection=yellow");
        IJ.run(imgNucZCrop,"Analyze Particles...", "  show=Masks exclude clear stack");
        ImagePlus imgNucMask = WindowManager.getImage("Mask of "+imgNucZCrop.getTitle());
        imgNucMask.hide();
        flush_close(imgNucZCrop);
        return(imgNucMask);    
    }
    
    
    /**
     * Find astrocyte nucleus keep nucleus 
     * that have the higher/lower ratio between intensity in astrocyte channel
     * and intensity in dots channel
     * @param nucPop
     * @param imgAstro
     * @param imgDots
     * @return nucAstro
     */
    
    public static int find_astroNuc( Objects3DPopulation nucPop, ImagePlus imgAstro, ImagePlus imgDots) {
        double maxRatioInt = 0;
        int index = 0;      
        if (nucPop.getNbObjects() > 1) {
            for (int o = 0; o < nucPop.getNbObjects(); o++) {
                Object3D obj = nucPop.getObject(o);
                double astroInt = obj.getPixMeanValue(ImageHandler.wrap(imgAstro));
                double dotsInt = obj.getPixMeanValue(ImageHandler.wrap(imgDots));
                double ratio;
                if (ratioInt)
                    ratio = dotsInt/astroInt;
                else
                    ratio = astroInt/dotsInt;
                if (ratio > maxRatioInt) {
                    maxRatioInt = ratio;
                    index = o;
                } 
            }   
        }
        return(index);
    }
    
    /**
     * compute local thickness
     * @param imgAstro
     * @return astroMap
    **/
    public static ImagePlus localThickness3D (ImagePlus imgAstro) {
        ImagePlus img = imgAstro.duplicate();
        img.setCalibration(cal);
        threshold(img, "Li", false, false);
        EDT_S1D edt = new EDT_S1D();
        edt.runSilent = true;
        edt.thresh = 1;
        edt.inverse = false;
        edt.showOptions = false;
        edt.setup("", img);
        edt.run(img.getProcessor());
        ImagePlus imgEDT = edt.getResultImage();
        imgEDT.setCalibration(cal);
        Local_Thickness_Parallel locThk = new Local_Thickness_Parallel();
        locThk.runSilent = true;
        locThk.setup("", imgEDT);
        locThk.run(imgEDT.getProcessor());
        ImagePlus imgLocThk = locThk.getResultImage();
        imgLocThk.setCalibration(cal);
        Clean_Up_Local_Thickness cleanUp = new Clean_Up_Local_Thickness();
        cleanUp.runSilent = true;
        cleanUp.setup("", imgLocThk);
        cleanUp.run(imgLocThk.getProcessor());
        ImagePlus astroMap = cleanUp.getResultImage();
        astroMap.setCalibration(cal);
        flush_close(img);
        flush_close(imgEDT);
        flush_close(imgLocThk);
        return(astroMap);
    }
    
    
    /**
     * Find Z with max intensity in stack
     * @param img
     * @return z
     */
    
    private static int find_max(ImagePlus img) {
        double max = 0;
        int zmax = 0;
        for (int z = 1; z <= img.getNSlices(); z++) {
            ImageProcessor ip = img.getStack().getProcessor(z);
            ImageStatistics statistics = new ImageStatistics().getStatistics(ip, ImageStatistics.MEAN, img.getCalibration());
            double meanInt = statistics.mean;
            if (meanInt > max) {
                max = meanInt;
                zmax = z;
            }
        }
        return(zmax);
    }
    
    
   /**
    * 3D watershed
    * @param imgBinmask
    * @param rad
    * @return distance map image
    */
    public static ImagePlus watershedSplit(ImagePlus imgBinmask, int rad) {
        int nbCpus = ThreadUtil.getNbCpus();
        float resXY = 1;
        float resZ = 1;
        float radXY = rad;
        float radZ = rad;
        if (cal != null) {
            resXY = (float) cal.pixelWidth;
            resZ = (float) cal.pixelDepth;
            radZ = radXY * (resXY / resZ);
        }
        IJ.showStatus("Computing EDT");
        ImageInt imgMask = ImageInt.wrap(imgBinmask);
        ImageFloat edt = EDT.run(imgMask, 1, resXY, resZ, false, nbCpus);
        ImageHandler edt16 = edt.convertToShort(true);
        ImagePlus edt16Plus = edt16.getImagePlus();
        GaussianBlur3D.blur(edt16Plus, 6.0, 6.0, 6.0);
        edt16 = ImageInt.wrap(edt16Plus);
        edt16.intersectMask(imgMask);
        // seeds
        ImageHandler seedsImg;
        seedsImg = FastFilters3D.filterImage(edt16, FastFilters3D.MAXLOCAL, radXY, radXY, radZ, 0, false);
        IJ.showStatus("Computing watershed");
        Watershed3D water = new Watershed3D(edt16, seedsImg, 10, 1);
        ImagePlus imp = water.getWatershedImage3D().getImagePlus();
        WindowManager.getWindow("Log").dispose();
        IJ.setThreshold(imp, 1, 65535);
        Prefs.blackBackground = false;
        IJ.run(imp, "Convert to Mask", "method=Default background=Dark");
        flush_close(imgMask.getImagePlus());
        flush_close(edt.getImagePlus());
        flush_close(edt16.getImagePlus());
        flush_close(edt16Plus);
        flush_close(seedsImg.getImagePlus());
        return(imp);
    }
    
    /**
     * Find dots population
     * @param imgDots
     * @param roi
     * @return dotsPop
     */
    
    public static Objects3DPopulation find_dots(ImagePlus imgDots, Roi roi, String thMethod) {
        ImagePlus imp = imgDots.duplicate();
        // substract bgDots 
        IJ.run(imp, "Subtract...", "value="+bgDots+" stack");
        IJ.run(imp, "Difference of Gaussians", "  sigma1=3 sigma2=1 stack");
        threshold(imp, thMethod, false, true);   
        if (roi != null) {
            imp.setRoi(roi);
            IJ.run(imp, "Clear Outside", "stack");
        }
        Objects3DPopulation dotsPop = getPopFromImage(imp);
        return(dotsPop);
    }
    
    
    /**
     * write object labels
     * @param obj
     * @param img
     * @param number
     * @param color
     * @param size
     */
    
    public static void labelsObject (Object3D obj, ImagePlus img, int number, int color, int size) {
        Font tagFont = new Font("SansSerif", Font.PLAIN, size);
        int[] box = obj.getBoundingBox();
        int z = (int)obj.getCenterZ();
        int x = box[0] - 2;
        int y = box[2] - 2;
        img.setSlice(z+1);
        ImageProcessor ip = img.getProcessor();
        ip.setFont(tagFont);
        ip.setColor(color);
        ip.drawString(Integer.toString(number), x, y);
        img.updateAndDraw();    
    }
    
   /**
    * draw population
    * nucleus in blue
    * large processes = 2 in green
    * not in astro = 0 in red
    * thin processes = 1 in yellow
    * @param nucObj
    * @param dotsPop
    * @param img
     * @param imgDir
     * @param imgName
     * @param roi
    */
    public static void tagsObjects(Object3D nucObj, Objects3DPopulation dotsPop, ImagePlus img, String imgDir, String imgName, int roi) {  
        ImagePlus imAstro = img.duplicate();
        ImageHandler imgObjNuc = ImageHandler.wrap(imAstro).createSameDimensions();
        nucObj.draw(imgObjNuc, 255);                        
        ImageHandler imgObjDotsNotAstro = ImageHandler.wrap(imAstro).createSameDimensions();
        ImageHandler imgObjDotsFineProcess = ImageHandler.wrap(imAstro).createSameDimensions();
        ImageHandler imgObjDotsLargeProcess = ImageHandler.wrap(imAstro).createSameDimensions();
        for (int n = 0; n < dotsPop.getNbObjects(); n++) {
            Object3D obj = dotsPop.getObject(n);
            switch (obj.getValue()) {
                case 0:
                    obj.draw(imgObjDotsNotAstro, 255);
                    break;
                case 1:
                    obj.draw(imgObjDotsFineProcess, 255);
                    break;
                default:
                    obj.draw(imgObjDotsLargeProcess, 255);
                    break;
            }
        }
        // save image for objects population
        ImagePlus[] imgColors = {imgObjDotsNotAstro.getImagePlus(), imgObjDotsLargeProcess.getImagePlus(), imgObjNuc.getImagePlus(),
            imAstro, null, null, imgObjDotsFineProcess.getImagePlus()};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
        imgObjects.updateAndDraw();
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(imgDir + imgName + "_Astro-"+(roi+1)+"-Objects.tif"); 
        flush_close(imgObjects);
    }
    
    
   
    /**
     * Flush and close images
     * @param img 
     */
    public static void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    } 
    
         
    /**
    * classify dots
    * @param nucObj nucleus
    * @param dotsPop dots population
     * @param imgAstro
    * @param imgAstroZcrop read dots intensity
    * @param astroMapZcrop
     * @param noAstroDot

    **/
    public static double classify_dots(Object3D nucObj, Objects3DPopulation dotsPop, ImagePlus imgAstro, ImagePlus imgAstroZcrop, ImagePlus astroMapZcrop) {
        IJ.showStatus("Classify dots ....");
        double bgThresholdInt = bg;
        double noAstroDot = 0;
        // if purearn exclude no dot
        // BgThresholdInt = -1
        if (pureArn)
            bgThresholdInt= -bg;
        // measure dots mean intensity and volume
        for (int n = 0; n < dotsPop.getNbObjects(); n++) {
            Object3D dotObj = dotsPop.getObject(n);
            double meanIntDotAstroimg = dotObj.getPixMeanValue(ImageHandler.wrap(imgAstroZcrop));
            
            /* classify dots
            * 
            * dots 0 (red) not in astro : 
            *   mean intensity dots <= bgThresholdInt 
            *   dist to nucleus > 2
            * dots 1 (yellow) Thin processes :
            *   mean intensity dots > bgThresholdInt, 
            *   dist to nucleus > 2, 
            *   largeProcess = false
            * dots 2 (green) Large processes and in soma :
            *   dist to nucleus < 2
            *   largeProcess = true
            */
            double distNuc = 0;
            if (nucObj.includes(dotObj))
                distNuc = 0;
            else
                distNuc = dotObj.distBorderUnit(nucObj);
            boolean largeProcess = largeProcess(astroMapZcrop, dotObj);
            // dots 0
            if (meanIntDotAstroimg <= bgThresholdInt && distNuc > 2) {
                dotObj.setValue(0);
            }

            // dots 1
            else if (meanIntDotAstroimg > bgThresholdInt && distNuc > 2 && !largeProcess) {
                dotObj.setValue(1);
            }
            // dots 2
            else if (distNuc <= 2 || largeProcess) {
                dotObj.setValue(2);
            }
            if (pureArn && meanIntDotAstroimg <= -bgThresholdInt && distNuc > 2) 
                noAstroDot++;
        }
        return(noAstroDot);
    }
    
    /**
     * Compute and write global image parameters
     * @param roi
     * @param roiIndex
     * @param totalRoi
     * @param imgAstro
     * @param imgAstroMap
     * @param nucObj
     * @param dotsPop
     * @param results
     * @param imageName
     * @param noAstroDot
     * @throws java.io.IOException
     */
    public static void compute_Image_parameters(Roi roi, int roiIndex, int totalRoi, ImagePlus imgAstro, ImagePlus imgAstroMap, Object3D nucObj, Objects3DPopulation dotsPop,
            BufferedWriter results, String imageName, double noAstroDot) throws IOException {
        /* write headers
        * largeBrDots (dots in large branches) = dots 2 - dotsinSoma
        * fineBrDots (dots in fine branches) = dots 1
        * dotsNoinAstro (dots outiside astrocyte) = dots 0
        * totalDots (all dots) = dots0 + dots1 + dots2
        * dotsinAstro (dots in astrocyte) = largeBrDots + fineBrDots + dotsinSoma
        * dotsinSoma (dots inside soma) = dots with dist to nucleus < 2
        * dotsinNuc (dots in nucleus) = dots with dist = 0
        * dotsdensinAstro (dots density in astrocyte) = dotsinAstro / astroVolume
        * percDotsNotinAstro (% dots not in astrocyte = dotsNotinAstro / totalDots
        * percDotsinSoma (% dots in soma) = dotsinSoma/totalDotsAstro
        * perDostFineBr (% dots in fine branches) = fineBrDots / totalDotsAstro
        * perDostLargeBr (% dots in large branches) = largeBrDots / totalDotsAstro
        */
        
        // measure roi volume
        double astroVol = roi_volume(roi, imgAstro);
        // Compute nucleus volume
        // add to astro parameters
        double nucVol = nucObj.getVolumeUnit(); 
        float dotsNoinAstro = 0;
        float dotsinFineBr = 0;
        float dotsinLargeBr = 0;
        float dotsinSoma= 0;
        DescriptiveStatistics dotVolStats = new DescriptiveStatistics();
        DescriptiveStatistics dotDiameterStats = new DescriptiveStatistics();
        DescriptiveStatistics dotMeanIntStats = new DescriptiveStatistics();
        String roiName = roi.getName();
        for (int n = 0; n < dotsPop.getNbObjects(); n++) {
            Object3D dotObj = dotsPop.getObject(n);
            double dotVol = dotObj.getVolumeUnit();
            dotVolStats.addValue(dotVol);
            double distNuc = 0;
            if (nucObj.includes(dotObj))
                distNuc = 0;
            else {
                distNuc = dotObj.distCenterBorderUnit(nucObj);
            }
            int Dottype = dotObj.getValue();
            double astroDiameter = dotObj.getPixMaxValue(ImageHandler.wrap(imgAstroMap)) * cal.pixelWidth;
            dotDiameterStats.addValue(astroDiameter);
            switch (Dottype) {
                case 0 :
                    if (dotVol >= meanSEMDotsSize) {
                        double ratioVol = Math.round(dotVol/meanSEMDotsSize);
                        dotsNoinAstro += ratioVol;
                    }
                    else
                        dotsNoinAstro++;
                    break;
                case 1 :
                    if (dotVol >= meanSEMDotsSize) {
                        double ratioVol = Math.round(dotVol/meanSEMDotsSize);
                        dotsinFineBr += ratioVol;
                    }
                    else
                        dotsinFineBr++;
                    break;
                case 2 : 
                    if (dotVol >= meanSEMDotsSize) {
                        double ratioVol = Math.round(dotVol/meanSEMDotsSize);
                        dotsinLargeBr += ratioVol;
                    }
                    else {
                        dotsinLargeBr++;
                        break;
                    }
            }
            if (distNuc < 2) 
                if (dotVol >= meanSEMDotsSize) {
                    double ratioVol = Math.round(dotVol/meanSEMDotsSize);
                    dotsinSoma+=ratioVol;  
                }   
                else 
                    dotsinSoma++;
            
            double IntdotsinAstro = dotObj.getPixMeanValue(ImageHandler.wrap(imgAstro));
            dotMeanIntStats.addValue(IntdotsinAstro);
        }
        // exclude dot in soma for large dots
        dotsinLargeBr = dotsinLargeBr - dotsinSoma;
        float totalDots = dotsinLargeBr + dotsinFineBr + dotsNoinAstro + dotsinSoma;
        float dotsinAstro = dotsinLargeBr + dotsinFineBr + dotsinSoma;
        double dotsdensinAstro = dotsinAstro / astroVol;
        double percDotsNotinAstro = 100*(dotsNoinAstro / totalDots);
        double percDotsinSoma = 100*(dotsinSoma / dotsinAstro);
        double percDotsFineBr = 100*(dotsinFineBr / dotsinAstro);
        double percDotsLargeBr = 100*((dotsinLargeBr) / dotsinAstro);
        double meanIntDotsinAstro = dotMeanIntStats.getMean();
        double meanAstroDotDiameter = dotDiameterStats.getMean();
        double percNoAstroDot = 100*(noAstroDot/totalDots);
        
        // find diameter statistics for image map
        DescriptiveStatistics diamStats = findDiameterStats(imgAstroMap);
        double astroMeanDiameter = diamStats.getMean();
        double astroStdDiameter = diamStats.getStandardDeviation();
        double astroMedianDiameter = diamStats.getPercentile(50);
        results.write(imageName + "\t" + roiName + "("+(roiIndex+1)+"/"+totalRoi+")" + "\t" + bg + 
                "\t" + astroVol + "\t" + dotsdensinAstro +
                "\t" + percDotsNotinAstro + "\t" + percDotsinSoma + "\t" + percDotsFineBr + "\t" + percDotsLargeBr +
                "\t" + meanIntDotsinAstro + "\t"+meanAstroDotDiameter+"\t"+astroMeanDiameter+"\t"+astroStdDiameter+"\t"+astroMedianDiameter+"\t"+percNoAstroDot+"\n");
        results.flush();
    }
    
    /**
    * For each image and for each roi write in results file 
    * dots type
    * dot volume
    * astrocyte diameter
     * @param imgAstroMap
     * @param imgAstro
     * @param dotsPop
     * @param astroNuc
     * @param roiName
     * @param results
     * @param imageName
     * @throws java.io.IOException
    */
    public static void computeDotsParams(ImagePlus imgAstroMap, ImagePlus imgAstro, Objects3DPopulation dotsPop, Object3D astroNuc, String roiName, 
            BufferedWriter results, String imageName) throws IOException {
        double astroDiam = 0;
        double dotInt = 0;
        double dotVol = 0; 
        double distToNuc = 0;
        int dotType = 0;
        for (int n = 0; n < dotsPop.getNbObjects(); n++) {
            Object3D obj = dotsPop.getObject(n);
            dotVol = obj.getVolumeUnit();
            astroDiam = obj.getPixMaxValue(ImageHandler.wrap(imgAstroMap)) * cal.pixelWidth;
            dotInt = obj.getPixMeanValue(ImageHandler.wrap(imgAstro));
            distToNuc = obj.distCenterBorderUnit(astroNuc);
            dotType = obj.getValue();
            results.write(imageName + "\t" + roiName + "\t" + n + "\t" + dotType + "\t" + dotVol + "\t" + dotInt + "\t" + astroDiam + "\t" + distToNuc + "\n");
            results.flush();
        }
    }
    
}
