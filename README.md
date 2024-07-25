# Astro_Dots

* **Developed for:** Katia
* **Team:** Cohen-Salmon
* **Date:** July 2024
* **Software:** Fiji

### Images description

3D images taken with a 63x objective.

3 channels:
  1. *440:* Nuclei
  2. *525:* Astrocytes
  3. *570:* RNA dots

With each image should be provided a *.roi* or *.zip* file containing one or multiple ROI(s).
A ROI, named *..._zStart_zStop*, should surround each astrocyte that needs to be analyzed.

### Plugin description

* Detect nuclei with Cellpose
* Detect astrocytes with median filtering + thresholding + median filtering 
* Detect RNA dots with DoG filtering + thresholding + median filtering
* Compute RNA dots channel background noise
* For each ROI provided:
  * Open dialog box allowing the user to select the proper astrocyte nucleus
  * Dilate the corresponding nucleus to get the astrocyte soma
  * Fill the obtained soma in black in astrocyte mask to get the mask of astrocyte processes only
  * Compute RNA dots volume and background-corrected intensity in each compartment: astrocyte nucleus / astrocyte processes / outside astrocyte

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ2** Fiji plugin
* **Cellpose** conda environment + *cyto2* model

### Version history

Version 2 released on July 25, 2024.
