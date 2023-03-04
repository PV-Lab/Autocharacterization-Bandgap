# Author: Alexander (Aleks) E. Siemenn <asiemenn@mit.edu>
# Date:   04 March 2023

# import
import cv2
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import ndimage
from spectral import *

def segmentation(bil, hdr, rotate_crop_params, lower_pixel_thresh=50, upper_pixel_thresh=3500):
    '''
    =======================================
    == VISION SEGMENTATION OF HYPERCUBE  ==
    =======================================

    Inputs:
        bil:                  Path to the .bil hyperspectral datacube file
        hdr:                  Path to the .hdr hyperspectral datacube file
        rotate_crop_params:   User-defined dictionary of values used to rotate and crop the image. E.g., {'theta': -0.5, 'x1': 45, 'x2': 830, 'y1': 120, 'y2': 550}
        lower_pixel_thresh:   Elements containing fewer than the specified number of pixels will be removed
        upper_pixel_thresh:   Elements containing more than the specified number of pixels will be removed

    Outputs:
        data:                 An array containing the reflectance data extracted from the vision-segmented hyperspectral datacube
    '''

    # load spectral hypercube
    hypercube = envi.open(hdr, bil).load()  # load datacude
    singleband = hypercube[:, :, 20]  # grab a single band to use for segmentation
    singleband = (singleband / singleband.max() * 255).astype("uint8")
    # plot raw image
    imshow(hypercube, bands=(154, 78, 27))  # R: 154, G: 78, B: 27
    plt.title('Raw Image')
    plt.show()
    # crop/rotate single band
    singleband = ndimage.rotate(singleband, rotate_crop_params['theta'])  # reads image and rotates
    singleband = singleband[rotate_crop_params['y1']:rotate_crop_params['y2'],
                 rotate_crop_params['x1']:rotate_crop_params['x2']]  # crops image
    # crop/rotate hypercube
    hypercube = ndimage.rotate(hypercube, rotate_crop_params['theta'])  # reads image and rotates
    hypercube = hypercube[rotate_crop_params['y1']:rotate_crop_params['y2'],
                rotate_crop_params['x1']:rotate_crop_params['x2']]  # crops image
    # plot cropped image
    imshow(hypercube, bands=(154, 78, 27))  # R: 154, G: 78, B: 27
    plt.title('Cropped Image')
    plt.show()
    # transpose to segment along print direction
    singlebandT = ndimage.rotate(singleband, -90)
    hypercubeT = ndimage.rotate(hypercube, -90)
    # segment image with watershed
    watershed, edges = watershed_segment(singlebandT, small_elements_pixels=lower_pixel_thresh,
                                         large_elements_pixels=upper_pixel_thresh)
    # plot segmented image
    fig, ax = plt.subplots()
    ax.imshow(watershed.T, cmap='viridis')
    ax.invert_yaxis()
    plt.title('Segmented Image')
    plt.show()
    # extract spectra from segmentation and format data
    idx_min = np.unique(watershed)[1:].min()
    idx_max = np.unique(watershed)[1:].max()
    idx_norms = []
    spectra_all = []
    for drop in np.unique(watershed)[1:]:  # iterate through all segmented droplets, excluding background=0
        cube = hypercubeT[watershed == drop, :]
        spectra = np.median(cube, axis=0)  # get median to avoid outliers
        spectra_all.append(spectra)
        idx_norm = (drop - idx_min) / (idx_max - idx_min)
        idx_norms.append(idx_norm)
    data = pd.DataFrame(np.concatenate(
        (np.array(pd.read_csv('wavelength.txt', sep='\t', header=None)[0]).reshape(1, -1), np.array(spectra_all))).T)
    data = data.rename(columns={0: 'wavelength'})
    # plot extracted reflectance data from segmented hypercube
    fig, ax = plt.subplots(figsize=(8, 4))
    colormap = mpl.cm.get_cmap('viridis')
    for n in range(len(spectra_all) - 1):  # plot lower peaks in front (FAPbI purple, MAPbI yellow)
        m = len(spectra_all) - 1 - n
        rgba = colormap(idx_norms[m])
        spectra = spectra_all[m]
        plt.plot(spectra / spectra.max(), c=rgba)  # plot normalization for visualization
    ax.minorticks_on()
    ax.grid(which='minor', color='gray', linestyle='--', alpha=0.35)
    ax.grid(which='major', color='gray', linestyle='--', alpha=0.35)
    ax.set_yticklabels([])
    plt.xlabel(r'Wavelength, $\lambda$ (nm)')
    plt.ylabel(r'Reflectance, $R$ (a.u.)')
    fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=100), cmap=mpl.cm.viridis), label='%MA')
    plt.title('Extracted Reflectance Spectra')
    plt.show()

    return data


# Sub-functions for segmentation()

def segment_on_dt(a, img, threshold):
    '''
    Implements watershed segmentation.

    Inputs:
        a:         Image input
        img:       Threshold binned image
        threshold: RGB threshold value

    Outputs:
        lbl:       Borders of segmented droplets
        wat:       Segmented droplets via watershed
        lab:       Indices of each segmented droplet
    '''
    # estimate the borders of droplets based on known and unknown background + foreground (computed using dilated and erode)
    kernel = np.ones((5, 5), np.uint8)
    border = cv2.dilate(img, None, iterations=1)
#     border = cv2.erode(border, kernel)
    border = border - cv2.erode(border, kernel)
    # segment droplets via distance mapping and thresholding
    dt = cv2.distanceTransform(img, 2, 3)
    dt = ((dt - dt.min()) / (dt.max() - dt.min()) * 255).astype(np.uint8)
    _, dt = cv2.threshold(dt, threshold, 255, cv2.THRESH_BINARY)
    # obtain the map of segmented droplets with corresponding indices
    lbl, ncc = ndimage.label(dt)
    lbl = lbl * (255 / (ncc + 1))
    lab = lbl
    # Completing the markers now.
    lbl[border == 255] = 255
    lbl = lbl.astype(np.int32)
    a = cv2.cvtColor(a,cv2.COLOR_GRAY2BGR)  # we must convert grayscale to BGR because watershed only accepts 3-channel inputs
    wat = cv2.watershed(a, lbl)
    lbl[lbl == -1] = 0
    lbl = lbl.astype(np.uint8)
    return 255 - lbl, wat, lab  # return lab, the segmented and indexed droplets

def watershed_segment(image, small_elements_pixels=0, large_elements_pixels=999999):
    '''
    Applies watershed image segmentation to separate droplet pixels from background pixels.

    Inputs:
        image:                   Input droplet image to segment
        small_elements_pixels:   Removes small elements that contain fewer than the specified number of pixels.
        large_elements_pixels:   Removes large elements that contain more than the specified number of pixels.

    Outputs:
        water:                   Watershed segmented droplets
        labs:                    Edges to each watershed segmented droplet
    '''
    RGB_threshold = 0
    img = image.copy()
    img = 255 - img
    _, img_bin = cv2.threshold(img, 0, 255,
                               # threshold image using Otsu's binarization # https://docs.opencv.org/4.x/d7/d4d/tutorial_py_thresholding.html
                               cv2.THRESH_OTSU)
    img_bin = cv2.morphologyEx(img_bin, cv2.MORPH_OPEN,
                               np.ones((12, 12), dtype=int))
    img_bin = cv2.erode(img_bin, np.ones((3, 3), np.uint8))
    img_bin = cv2.medianBlur(img_bin, 7)
    result, water, labs = segment_on_dt(a=img, img=img_bin,
                                        threshold=RGB_threshold)  # segment droplets from background and return indexed droplets
    water = cv2.dilate(water.astype('uint8'), np.ones((5, 5), np.uint8))
    water = cv2.medianBlur(water,5)
    # remove small/large elements
    uniq_full, uniq_counts = np.unique(water,
                                       return_counts=True)  # get all unique watershed indices with pixel counts
    large_elements = uniq_full[uniq_counts > large_elements_pixels]  # mask large elements based on number of pixels
    small_elements = uniq_full[uniq_counts < small_elements_pixels] # mask small elements based on number of pixels
    for n in range(len(large_elements)):
        water[water == large_elements[n]] = 0  # remove all large elements
    for n in range(len(small_elements)):
        water[water == small_elements[n]] = 0  # remove all small elements
    return water, labs
