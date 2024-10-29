package Astro_Dots_Tools;

import Astro_Dots_Tools.Cellpose.CellposeTaskSettings;
import Astro_Dots_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.*;
import ij.plugin.Duplicator;
import ij.plugin.ImageCalculator;
import ij.plugin.RGBStackMerge;
import ij.plugin.filter.Analyzer;
import ij.plugin.frame.RoiManager;
import ij.process.AutoThresholder;
import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Object3DPlane;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.VoxelInt;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;
import sc.fiji.localThickness.LocalThicknessWrapper;


/**
 * @author ORION-CIRB
 */
public class Tools {
    
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    private final String helpUrl = "https://github.com/orion-cirb/Astro_Dots/tree/version2";
    
    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
    private Calibration cal;
    public String[] channelsName = new String[]{"Nuclei", "Astrocytes", "RNA dots"};
    
    // Nuclei detection
    public String cellposeEnvDir = IJ.isWindows()? System.getProperty("user.home")+File.separator+"miniconda3"+File.separator+"envs"+File.separator+"CellPose" : "/opt/miniconda3/envs/cellpose";
    public String cellposeModel = "cyto2";
    public int cellposeDiam = 80;
    public double cellposeStitchTh = 0.75;
    public double minNucVol = 50;// µm3
    public double maxNucVol = 1200; // µm3
    public double nucDilation = 2; // µm
    
    // Astrocytes segmentation
    public String astroThMethod = "Huang";
    
    // RNA dots segmentation and classification
    public String dotsThMethod = "Triangle";
    private double dotEstimatedVol = 0.003; // µm3
   

    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom2.Object3DInt");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Flush and close an image
     */
    public void closeImage(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Get extension of the first image found in the folder
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd" :
                   ext = fileExt;
                   break;
                case "nd2" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }
    
    
    /**
     * Find images in folder
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    
    /**
     * Find image calibration
     */
    public void findCalibration(IMetadata meta) {
        cal = new Calibration();
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
    }
    
    
    /**
     * Find channels name and add None at the end of channels list
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels(String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);     
    }

    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 60, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        for (int n = 0; n < channelsName.length; n++) {
            gd.addChoice(channelsName[n], channels, channels[n]);
        }
        
        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min volume (µm3): ", minNucVol, 2);
        gd.addNumericField("Max volume (µm3): ", maxNucVol, 2);
        gd.addNumericField("Dilation (soma) (µm):", nucDilation, 2);

        gd.addMessage("Astrocytes segmentation", Font.getFont("Monospace"), Color.blue);
        String[] thMethods = AutoThresholder.getMethods();  
        gd.addChoice("Thresholding method: ", thMethods, astroThMethod);       
        
        gd.addMessage("RNA dots segmentation", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("Thresholding method:", thMethods, dotsThMethod);
        gd.addNumericField("Dot estimated volume (µm3)", dotEstimatedVol, 4);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY pixel size (µm):", cal.pixelWidth, 4);
        gd.addNumericField("Z pixel size (µm):", cal.pixelDepth, 4);
        gd.addHelp(helpUrl);
        gd.showDialog();
        
        String[] channelsOrdered = new String[channelsName.length];
        for (int n = 0; n < channelsOrdered.length; n++) {
            channelsOrdered[n] =  gd.getNextChoice();
        }
        
        minNucVol = gd.getNextNumber();
        maxNucVol = gd.getNextNumber();
        nucDilation = gd.getNextNumber();
        
        astroThMethod = gd.getNextChoice();
        
        dotsThMethod = gd.getNextChoice();
        dotEstimatedVol = gd.getNextNumber();
        
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        
        if (gd.wasCanceled())
            return(null);
        
        return(channelsOrdered);
    }
    
    
    /**
     * Load ROIs
     */
    public List<Roi> loadRois(String imageDir, String rootName) {
        String roiName = imageDir+rootName;
        roiName = new File(roiName+".zip").exists() ? roiName+".zip" : roiName+".roi";
        List<Roi> rois = new ArrayList<>();
        
        if (new File(roiName).exists()) {
            RoiManager rm = new RoiManager(false);
            rm.runCommand("Open", roiName);
            List<Roi> roisTemp = Arrays.asList(rm.getRoisAsArray());
            
            for(Roi roi: roisTemp) {
                String[] name = roi.getName().split("-");
                if(name.length == 3) {
                    // Find in ROI name the desired top and bottom stack
                    roi.setProperty("zStart", name[1]);
                    roi.setProperty("zStop", name[2]);
                    rois.add(roi);
                } else {
                    print("WARNING: ROI " + roiName + " should be named roiNumber-zTop-zBottom, ROI not analyzed");
                }
            }
            return(rois);
        } else {
            print("WARNING: No ROI file found for image " + rootName + ", image not analyzed");
            return(null);
        }
    }

    
    /*
     * Detect 3D cells in a Z-stack: 
     * - apply CellPose in 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitching threshold parameter
     */
    public ImagePlus cellposeDetection(ImagePlus imgIn) throws IOException{
       ImagePlus img = imgIn.duplicate();
       
       // Define CellPose settings
       CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModel, 1, cellposeDiam, cellposeEnvDir);
       settings.setStitchThreshold(cellposeStitchTh);
       settings.useGpu(true);
       
        // Run CellPose
       CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, img);
       ImagePlus imgOut = cellpose.run();
       imgOut.setCalibration(cal);
       
       closeImage(img);
       return(imgOut);
    }
    
    
    /**
     * Segment astrocytes with median filtering + thresholding + median filtering
     */
    public ImagePlus astroSegmentation(ImagePlus img) {
        ImagePlus imgMed = medianFilter(img, 2, 1);
        ImagePlus imgBin = threshold(imgMed, astroThMethod);
        ImagePlus imgOut = medianFilter(imgBin, 1, 1);
        imgOut.setCalibration(cal);
        
        closeImage(imgMed);
        closeImage(imgBin);
        return(imgOut);
    }
    
    
    /**
     * Segment RNA dots with DoG filtering + thresholding  + median filtering
     */
    public ImagePlus dotsSegmentation(ImagePlus img) {
        ImagePlus imgDog = DOG(img, 1, 2);
        ImagePlus imgBin = threshold(imgDog, dotsThMethod);
        ImagePlus imgOut = medianFilter(imgBin, 1, 1);
        imgOut.setCalibration(cal);
        
        closeImage(imgDog);
        closeImage(imgBin);
        return(imgOut);
    }
    
    
    /**
     * 3D median filtering using CLIJ2
     */ 
    public ImagePlus medianFilter(ImagePlus img, double sizeXY, double sizeZ) {
        ClearCLBuffer imgCL = clij2.push(img); 
        ClearCLBuffer imgCLMed = clij2.create(imgCL);
        clij2.median3DSphere(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
        ImagePlus imgMed = clij2.pull(imgCLMed);
        clij2.release(imgCL);
        clij2.release(imgCLMed);
        return(imgMed);
    }

    
    /**
     * Difference of Gaussians using CLIJ
     */ 
    public ImagePlus DOG(ImagePlus img, double sizeXY1, double sizeXY2) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian2D(imgCL, imgCLDOG, sizeXY1, sizeXY1, sizeXY2, sizeXY2);
        clij2.release(imgCL);
        ImagePlus imgDOG = clij2.pull(imgCLDOG); 
        clij2.release(imgCLDOG);
        return(imgDOG);
    }
    
    
    /**
     * Automatic thresholding using CLIJ2
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCL);
        clij2.release(imgCLBin);
        
        IJ.run(imgBin, "8-bit", "");
        return(imgBin);
    }
    
    
    /**
     * Compute image background noise
     */
    public double computeBackgroundNoise(ImagePlus img, ImagePlus mask) {
        Objects3DIntPopulation popDots = new Objects3DIntPopulation(ImageHandler.wrap(mask));
        ImageHandler imhBg = ImageHandler.wrap(img.duplicate());
        for(Object3DInt dot: popDots.getObjects3DInt())
            dot.drawObject(imhBg, 0);

        ImagePlus imgBg = imhBg.getImagePlus();
        double totalInt = 0;
        int nbZeroPixels = 0;
        for(int z = 0; z < imgBg.getNSlices(); ++z) {
            for(int x = 0; x < imgBg.getWidth(); ++x) {
                for(int y = 0; y < imgBg.getHeight(); ++y) {
                    double pixelInt = imhBg.getPixel(x, y, z);
                    if(pixelInt != 0)
                        totalInt += pixelInt;
                    else
                        nbZeroPixels++;
                }
            }
        }
        double bg = totalInt / (imgBg.getHeight() * imgBg.getWidth() * imgBg.getNSlices() - nbZeroPixels);
        System.out.println("Background noise = " + bg);
        
        closeImage(imgBg);
        return(bg);
    }
    
    
    public ImagePlus cropImage(ImagePlus imgIn, Roi roi, int zStart, int zStop) {
        imgIn.setRoi(roi);
        ImagePlus imgOut = new Duplicator().run(imgIn, zStart, zStop);
        imgIn.deleteRoi();
        return(imgOut);
    }
    
    
    /**
     * Return 3D object in ROI
     */
    public ImagePlus getImageInRoi(ImagePlus mask, Roi roi, int zStart, int zStop) throws IOException{
        ImagePlus maskCrop = cropImage(mask, roi, zStart, zStop);
        IJ.setBackgroundColor(0, 0, 0);
        Rectangle rect = roi.getBounds();
        roi.setLocation(0, 0);
        maskCrop.setRoi(roi);
        IJ.run(maskCrop, "Clear Outside", "stack");
        roi.setLocation(rect.x, rect.y);
        maskCrop.setCalibration(cal);
        return(maskCrop);
    }
    
    
    /**
     * Return 3D cells population in ROI
     */
    public Objects3DIntPopulation getPopulationInRoi(ImagePlus mask, Roi roi, int zStart, int zStop, double minVol, double maxVol) throws IOException{
        ImagePlus maskCrop = getImageInRoi(mask, roi, zStart, zStop);

        // Filter detections
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageInt.wrap(maskCrop));
        System.out.println(pop.getNbObjects() + " objects detected in ROI");
        zFilterPop(pop);
        sizeFilterPop(pop, minVol, maxVol);
        System.out.println(pop.getNbObjects() + " objects remaining after size filtering");
        pop.resetLabels();

        closeImage(maskCrop);
        return(pop);
    }
    

    /**
     * Remove objects that appear in only one z-slice
     */
    public void zFilterPop(Objects3DIntPopulation pop) {
        pop.getObjects3DInt().removeIf(p -> (p.getObject3DPlanes().size() == 1));
    }
    
    
    /**
     * Remove objects from population with size outside of given range
     */
    public void sizeFilterPop(Objects3DIntPopulation pop, double minVol, double maxVol) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < minVol) || (new MeasureVolume(p).getVolumeUnit() > maxVol));
    }
    

    /**
     * Display detected astrocyte nucleus in green, others in blue
     * User can select the proper nucleus in a dialog box
     */
    public Object3DInt nucleusSelection(Objects3DIntPopulation nucPop, ImagePlus imgNuc, ImagePlus imgAstro, ImagePlus imgDots) {
        Object3DInt preSelectedNuc = detectAstrocyteNucleus(nucPop);
        if(preSelectedNuc == null)
            return(null);
        
        ImageHandler imhNuc = ImageHandler.wrap(imgAstro).createSameDimensions();
        nucPop.drawInImage(imhNuc);
        IJ.run(imhNuc.getImagePlus(), "Enhance Contrast", "saturated=0.35");
        
        ImageHandler imhNucPreSelected = imhNuc.createSameDimensions();
        preSelectedNuc.drawObject(imhNucPreSelected);
        IJ.run(imhNucPreSelected.getImagePlus(), "Enhance Contrast", "saturated=0.35");

        ImagePlus[] imgColors = {null, imhNuc.getImagePlus(), imhNucPreSelected.getImagePlus(), imgAstro, imgNuc};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgObjects.setZ(preSelectedNuc.getObject3DPlanes().get(0).getZPlane()+1);
        imgObjects.show("");
        
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("");
        gd.addMessage("Select proper astrocyte nucleus:");
        String[] radioItems = new String[nucPop.getNbObjects()];
        for(int i = 0; i < nucPop.getNbObjects(); i++) {
            radioItems[i] = Integer.toString(i+1);
        }
        gd.addRadioButtonGroup("", radioItems, nucPop.getNbObjects(), 1, Integer.toString((int)preSelectedNuc.getLabel()));
        gd.showDialog();
        
        imgObjects.hide();
        Object3DInt selectedNuc = nucPop.getObjectByLabel(Float.parseFloat(gd.getNextRadioButton()));
        
        closeImage(imhNuc.getImagePlus());
        closeImage(imhNucPreSelected.getImagePlus());
        closeImage(imgObjects);
        
        if(gd.wasCanceled())
            return(null);
        
        return(selectedNuc);
    }
    
    
    /**
     * Return dilated object restricted to image borders
     */
    public Object3DInt dilateObject(Object3DInt obj, ImagePlus img, double dilSize) {
        Object3DInt objDil = new Object3DComputation(obj).getObjectDilated((float)(dilSize/cal.pixelWidth), (float)(dilSize/cal.pixelHeight),(float)(dilSize/cal.pixelDepth));
        
        // Check if object goes over image borders
        BoundingBox bbox = objDil.getBoundingBox();
        BoundingBox imgBbox = new BoundingBox(ImageHandler.wrap(img));
        int[] box = {imgBbox.xmin, imgBbox.xmax, imgBbox.ymin, imgBbox.ymax, imgBbox.zmin, imgBbox.zmax};
        if (bbox.xmin < 0 || bbox.xmax > imgBbox.xmax || bbox.ymin < 0 || bbox.ymax > imgBbox.ymax || bbox.zmin < 0 || bbox.zmax > imgBbox.zmax) {
            Object3DInt objDilInImg = new Object3DInt();
            for (Object3DPlane p: objDil.getObject3DPlanes()) {
                for (VoxelInt v: p.getVoxels()) {
                    if (v.isInsideBoundingBox(box))
                        objDilInImg.addVoxel(v);
                }
            }
            objDilInImg.setVoxelSizeXY​(cal.pixelWidth);
            objDilInImg.setVoxelSizeZ(cal.pixelDepth);
            return(objDilInImg);
        } else {
            objDil.setVoxelSizeXY​(cal.pixelWidth);
            objDil.setVoxelSizeZ(cal.pixelDepth);
            return(objDil);
        }
    }
    
    
    /**
     * Detect astrocyte nucleus as the biggest one
     */
    public Object3DInt detectAstrocyteNucleus(Objects3DIntPopulation nucPop) {
        float label = 0;
        double maxVol = 0;
        
        for (Object3DInt nuc: nucPop.getObjects3DInt()) {
            double vol = new MeasureVolume(nuc).getVolumePix();
            if (vol > maxVol) {
                maxVol = vol;
                label = nuc.getLabel();
            }
        }
        
        if(label == 0)
            return(null);
        
        return(nucPop.getObjectByLabel(label));
    }
    
    
    /**
     * Compute the local thickness of a 3D image stack
     */
    public ImagePlus localThickness3D(ImagePlus img) {
        System.out.println("Computing 3D local thickness...");
        
        IJ.run(img, "8-bit", "");
        img.setCalibration(cal);
        
        LocalThicknessWrapper locThick = new LocalThicknessWrapper();
        locThick.setup("", img);
        locThick.run(img.getProcessor());
        ImagePlus imgLocThick = locThick.getResultImage();
        imgLocThick.hide();
        
        return(imgLocThick);
    }
    
    
    /*
     * Classify dots in 3 categories: inside soma / inside processes / outside astrocyte
     */
    public List<Object3DInt> dotsClassification(ImagePlus maskDots, Object3DInt objSoma, Object3DInt objAstro) {
        System.out.println("Classifying dots...");
        
        ImageHandler imhOut = ImageHandler.wrap(maskDots.duplicate());
        objSoma.drawObject(imhOut, 0);
        objAstro.drawObject(imhOut, 0);
        Object3DInt objDotsOut = new Object3DInt(imhOut);
        
        ImagePlus imgIn = new ImageCalculator().run("subtract stack create", maskDots, imhOut.getImagePlus());
        
        ImageHandler imhSoma = ImageHandler.wrap(imgIn.duplicate());
        objAstro.drawObject(imhSoma, 0);
        Object3DInt objDotsSoma = new Object3DInt(imhSoma);
        
        ImageHandler imhProcess = ImageHandler.wrap(imgIn.duplicate());
        objSoma.drawObject(imhProcess, 0);
        Object3DInt objDotsProcess = new Object3DInt(imhProcess);

        closeImage(imhOut.getImagePlus());
        closeImage(imhSoma.getImagePlus());
        closeImage(imhProcess.getImagePlus());
        closeImage(imgIn);
        
        return(Arrays.asList(objDotsSoma, objDotsProcess, objDotsOut));  
        
    }
    
  
    /**
     * Compute and write results
     */
    public void writeResults(Roi roi, Object3DInt objSoma, Object3DInt objAstro, List<Object3DInt> objsDots, double bgDots, 
                             ImagePlus imgDots, BufferedWriter results, String imgName) throws IOException {
        
        double roiVol = getRoiVolume(roi, imgDots);
        double somaVol = new MeasureVolume(objSoma).getVolumeUnit();
        double processVol = new MeasureVolume(objAstro).getVolumeUnit();
        ImageHandler imhDots = ImageHandler.wrap(imgDots);
        double dotsSomaVol = new MeasureVolume(objsDots.get(0)).getVolumeUnit();
        double dotsSomaInt = new MeasureIntensity(objsDots.get(0), imhDots).getValueMeasurement(MeasureIntensity.INTENSITY_SUM) - bgDots * new MeasureVolume(objsDots.get(0)).getVolumePix();
        double dotsProcessVol = new MeasureVolume(objsDots.get(1)).getVolumeUnit();
        double dotsProcessInt = new MeasureIntensity(objsDots.get(1), imhDots).getValueMeasurement(MeasureIntensity.INTENSITY_SUM) - bgDots * new MeasureVolume(objsDots.get(1)).getVolumePix();
        double dotsOutVol = new MeasureVolume(objsDots.get(2)).getVolumeUnit();
        double dotsOutInt = new MeasureIntensity(objsDots.get(2), imhDots).getValueMeasurement(MeasureIntensity.INTENSITY_SUM) - bgDots * new MeasureVolume(objsDots.get(2)).getVolumePix();
        
        results.write(imgName+"\t"+roi.getName()+"\t"+roiVol+"\t"+somaVol+"\t"+processVol+"\t"+bgDots+
                      "\t"+dotsSomaVol+"\t"+(int)Math.round(dotsSomaVol/dotEstimatedVol)+"\t"+dotsSomaInt+
                      "\t"+dotsProcessVol+"\t"+(int)Math.round(dotsProcessVol/dotEstimatedVol)+"\t"+dotsProcessInt+
                      "\t"+dotsOutVol+"\t"+(int)Math.round(dotsOutVol/dotEstimatedVol)+"\t"+dotsOutInt+"\n");
        results.flush();
    }
    
    
    /**
     * Compute ROI volume
     */
    public double getRoiVolume(Roi roi, ImagePlus img) {
        roi.setLocation(0, 0);
        img.setRoi(roi);
        ResultsTable rt = new ResultsTable();
        Analyzer analyzer = new Analyzer(img, Analyzer.AREA, rt);
        analyzer.measure();
        return(rt.getValue("Area", 0) * img.getNSlices() * cal.pixelDepth);
    }
    

    /*
     * Draw results 
     */
    public void drawResults(Object3DInt objSoma, ImagePlus maskAstro, List<Object3DInt> objsDots, ImagePlus imgNuc,
            ImagePlus imgAstro, ImagePlus imgDots, ImagePlus locThickAstro, String imgDir, String imgName, String roiName) {
        
        IJ.run(locThickAstro, "Calibrate...", "function=None unit=µm");
        IJ.run(locThickAstro, "Calibration Bar...", "location=[Upper Left] fill=White label=Black number=5 decimal=2 font=13 zoom=0.8 overlay show");
        new FileSaver(locThickAstro).saveAsTiff(imgDir + imgName + "_" + roiName + "_diam.tif");
        
        ImageHandler imhSoma = ImageHandler.wrap(maskAstro).createSameDimensions();
        objSoma.drawObject(imhSoma, 255);
        ImageHandler imhDotsSoma = imhSoma.createSameDimensions();
        objsDots.get(0).drawObject(imhDotsSoma, 255);
        ImageHandler imhDotsProcess = imhSoma.createSameDimensions();
        objsDots.get(1).drawObject(imhDotsProcess, 255);
        ImageHandler imhDotsOut = imhSoma.createSameDimensions();
        objsDots.get(2).drawObject(imhDotsOut, 255);
        
        ImagePlus[] resultsAstro = {null, maskAstro, imhSoma.getImagePlus(), imgAstro, imgNuc};
        ImagePlus imgResultsAstro = new RGBStackMerge().mergeHyperstacks(resultsAstro, true);
        imgResultsAstro.setCalibration(cal);
        imgResultsAstro.setC(3);
        IJ.run(imgResultsAstro, "Enhance Contrast", "saturated=0.35");
        imgResultsAstro.setC(4);
        IJ.run(imgResultsAstro, "Enhance Contrast", "saturated=0.35");
        new FileSaver(imgResultsAstro).saveAsTiff(imgDir + imgName + "_" + roiName + "_astro.tif");
        
        ImagePlus[] resultsDots = {imhDotsSoma.getImagePlus(), imhDotsProcess.getImagePlus(), imhDotsOut.getImagePlus(), imgAstro, imgDots};
        ImagePlus imgResultsDots = new RGBStackMerge().mergeHyperstacks(resultsDots, true);
        imgResultsDots.setCalibration(cal);
        imgResultsDots.setC(4);
        IJ.run(imgResultsDots, "Enhance Contrast", "saturated=0.35");
        imgResultsDots.setC(5);
        IJ.run(imgResultsDots, "Enhance Contrast", "saturated=0.35");
        new FileSaver(imgResultsDots).saveAsTiff(imgDir + imgName + "_" + roiName + "_dots.tif");
        
        closeImage(imhSoma.getImagePlus());
        closeImage(imhDotsSoma.getImagePlus());
        closeImage(imhDotsProcess.getImagePlus());
        closeImage(imhDotsOut.getImagePlus());
        closeImage(imgResultsAstro);
        closeImage(imgResultsDots);
    }
    
}
