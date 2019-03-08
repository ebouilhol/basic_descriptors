#!/usr/bin/python

import skimage
import numpy as np
from skimage import measure
from skimage import feature
from skimage import exposure
from skimage import filters
from skimage import morphology
from skimage.morphology import disk, square
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from matplotlib.patches import Circle
from skimage.filters.rank import entropy
import os
import skimage.io as io
from skimage.feature import blob_log
import csv

def get_max_projection(image_3d):
    num_slices = len(image_3d[:,0,0])
    max_proj = image_3d[0,:,:]
    for slice_idx in range(1, num_slices):
        slice = image_3d[slice_idx,:,:]
        max_proj = np.maximum(max_proj, slice)
    min = np.min(max_proj)
    max_proj[max_proj < min + np.sqrt(min)]=0
    return max_proj


def vignette(image):
    # print image.shape
    width, height = image.shape[:2]
    xs = np.arange(width)
    ys = np.arange(height)
    distance_squared = (xs - width / 3.0)[..., np.newaxis] ** 2 + (ys - height / 2.0) ** 2
    # print distance_squared
    sigma_squared = (width / 2.0) ** 2 + (height / 2.0) ** 2
    # print sigma_squared
    falloff = np.exp(-distance_squared / sigma_squared)
    result = image * falloff
    return result


def max_contour(contours):
    max = 0
    pos = 0
    for n, contour in enumerate(contours):
        if len(contour) > max:
            max = len(contour)
            pos = n
    return contours[pos]

def morpho_math(image, dilatation, closing, erosion):
    image = skimage.morphology.binary_dilation(image, disk(dilatation))
    image = skimage.morphology.closing(image, square(closing))
    image = ndi.binary_fill_holes(image)
    image = skimage.morphology.binary_erosion(image, disk(erosion))
    return image

def find_edges_canny(image, disk, gabor, canny, sigma):
    image = image.astype(int)
    image = skimage.filters.rank.median(image, disk(disk))
    image, toto = skimage.filters.gabor(image, gabor)
    image = skimage.feature.canny(image, sigma=sigma)
    return image

def find_edges_canny_dapi(image, sigma):
    image = skimage.feature.canny(image, sigma=sigma)
    return image


def find_edges_canny_dapi_prot(image, sigma, low, high):
    image = skimage.feature.canny(image, sigma=sigma,low_threshold=low, high_threshold=high, use_quantiles=True)
    return image


def stretch_hist(tubulin):
    tubulin = tubulin.astype(int)
    bins_percent=[1,10,20,30,40,50,60,70,80,90,100]
    percent = np.percentile(tubulin, bins_percent)
    #etirement histogramme
    tubulin = ((tubulin - percent[9])*255) / (percent[10]-percent[9])
    tubulin = skimage.exposure.rescale_intensity(tubulin, out_range=(0,255))
    return tubulin

def tubulin_treatment(tubulin_3d, case, val, entropy_size, dilatation, closing, erosion, sigma, low, high):
    # tubulin = get_max_projection(tubulin_3d)
    tubulin = tubulin_3d
    tubulin = vignette(tubulin)
    tubulin=stretch_hist(tubulin)
    tubulin = tubulin.astype(int)
    tubulin= entropy(tubulin, square(entropy_size))

    ###### non mandatory step
    if case == 1:
        bins_percent=[1,10,20,30,40,50,60,75,76,val,100]
        percent = np.percentile(tubulin, bins_percent)
        tubulin[tubulin <= percent[9]] = 0
        ######

    tubulin = find_edges_canny_dapi_prot(tubulin, sigma, low, high)
    tubulin = morpho_math(tubulin, dilatation, closing, erosion)
    return tubulin



#if case = 0 then val is not needed, put the value you want !
def dapi_treatment(dapi, case, val, dilatation, closing, erosion, sigma, low, high):
    if case == 0:
        dapi = skimage.exposure.rescale_intensity(dapi)
        dapi_threshold = skimage.filters.threshold_otsu(dapi)
        dapi[dapi <= dapi_threshold] = 0
        # dapi = find_edges_canny_dapi_prot(dapi, sigma, low, high)
        dapi = morpho_math(dapi, dilatation, closing, erosion)
    elif case == 1:
        dapi = dapi.astype(int)
        dapi= entropy(dapi, square(20))
        bins_percent=[1,10,20,30,40,50,60,75,76,val,100]
        percent = np.percentile(dapi, bins_percent)
        dapi[dapi <= percent[9]] = 0
        dapi = skimage.feature.canny(dapi, sigma=sigma, low_threshold=low, high_threshold=high, use_quantiles=True)
    dapi = morpho_math(dapi, dilatation, closing, erosion)
    return dapi

def save_mask(name, image):
    skimage.io.imsave(os.path.join(name+".png"), image)


def plot_nucleus_contours(raw_dapi, dapi_contour, x,y, dapi_out_file_name):
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize = (8, 8))
    # ax.imshow(tubulin_keep, interpolation='nearest', cmap=plt.cm.gray)
    ax.imshow(raw_dapi, interpolation='nearest', cmap=plt.cm.gray)
    ax.plot(dapi_contour[:, 1], dapi_contour[:, 0], linewidth=2, color = 'b')
    ax.axis('image')
    centroid = Circle((x, y), 3, color=(0, 0, 1)) # small circle at the center (256, 256)
    ax.add_patch(centroid)
    plt.xlim((0,512))
    plt.ylim((512, 0))
    # plt.show()
    fig.savefig(os.path.join(dapi_out_file_name+".png"))
    plt.close()

def plot_cell_contours(raw_tubulin, tubulin_contour, tubulin_out_file_name):
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize = (8, 8))
    ax.imshow(raw_tubulin, interpolation='nearest', cmap=plt.cm.gray)
    ax.plot(tubulin_contour[:, 1], tubulin_contour[:, 0], linewidth=2, color = 'r')
    ax.axis('image')
    plt.xlim((0,512))
    plt.ylim((512, 0))
    # plt.show()
    fig.savefig(tubulin_out_file_name)
    plt.close()


def get_nucleus_centroid(dapi_regions, dapi_mask):
    max_region = dapi_regions[0]
    max_region_area = dapi_regions[0].area;
    for region in dapi_regions:
        if region.area > max_region_area:
            max_region = region
            max_region_area = region.area

    cell_mask = dapi_mask
    for region in dapi_regions:
        if region.area < max_region_area:
            start_x, start_y, end_x, end_y = region.bbox
            for x in range(start_x, end_x):
                for y in range(start_y, end_y):
                    cell_mask[x,y] = 0

    # y, x = max_region.centroid
    return max_region.centroid

def get_spots(img, mask, h5_grp):
    tab_spot_verif = []
    sobel = filters.sobel(img)
    blurred = filters.gaussian(sobel, sigma=1.0)
    tophat = morphology.white_tophat(blurred, square(3))
    rescale = exposure.rescale_intensity(tophat, out_range=(0, 255))
    blobs_log = blob_log(rescale, min_sigma=4, max_sigma=5, num_sigma=1, threshold=3.5)
    for blob in blobs_log:
        x = int(blob[0])
        y = int(blob[1])
        r = int(blob[2])
        if (mask[y, x] == 1):
            tab_spot_verif.append([x, y, r])
    h5_grp.create_dataset('spots', data=np.array(tab_spot_verif))


def add_mtoc_from_csv(mtoc_file, mask, root, h5_grp):
    with open(mtoc_file, 'rb') as f:
        reader = csv.reader(f)
        mtoc_pos = list(reader)
    mtocs_names = [row[0] for row in mtoc_pos]
    mtocs_x = [int(row[1]) for row in mtoc_pos]
    mtocs_y = [int(row[2]) for row in mtoc_pos]
    folder_name = os.path.basename(os.path.normpath(root))
    pos = [i for i, val in enumerate(mtocs_names) if val == folder_name]
    tab_mtoc = []
    for i in pos:
        # if (mask[mtocs_y[i], mtocs_x[i]] == 0):
        tab_mtoc.append([mtocs_x[i], mtocs_y[i]])
    if len(tab_mtoc) != 0:
        h5_grp.create_dataset('mtoc_position', data=tab_mtoc)



def dapi_all_slices(img_3d, dapi_case, dapi_treshold, dapi_mask_output_name, dapi_out_file_name, dapi_dilatation, dapi_closing, dapi_erosion, sigma, low, high, h5_grp):

    img = get_max_projection(img_3d)
    dapi_mask = dapi_treatment(img, dapi_case, dapi_treshold, dapi_dilatation, dapi_closing, dapi_erosion, sigma, low, high)
    dset3 = h5_grp.create_dataset('nucleus_mask', data=np.array(dapi_mask).astype(np.int8))
    # io.imshow(dapi_mask)
    # io.show()
    save_mask(dapi_mask_output_name, dapi_mask)
    dapi_contours = measure.find_contours(dapi_mask, 0.8)
    dapi_contour = max_contour(dapi_contours)
    dapi_label = measure.label(dapi_mask)
    dapi_regions = measure.regionprops(dapi_label)
    y, x  = get_nucleus_centroid(dapi_regions, dapi_mask)
    h5_grp.create_dataset('nucleus_centroid', data=[x, y])

    plot_nucleus_contours(img, dapi_contour, x, y, dapi_out_file_name)


def tubulin_all_slices(root, img_3d, case, treshold, entropy_size, mask_output_name, out_file_name, dilatation, closing, erosion, sigma, low, high, h5_grp, path):
    img = get_max_projection(img_3d)
    cell_mask = tubulin_treatment(img, case, treshold, entropy_size, dilatation, closing, erosion, sigma, low, high)
    h5_grp.create_dataset("cell_mask", data=np.array(cell_mask).astype(np.int8))
    # io.imshow(cell_mask)
    # io.show()
    save_mask(mask_output_name, cell_mask)
    cell_contours = measure.find_contours(cell_mask, 0.8)
    cell_contour = max_contour(cell_contours)
    add_mtoc_from_csv(os.path.join('mtoc', path+'.csv'), cell_mask, root, h5_grp)
    get_spots(img, cell_mask, h5_grp)

    plot_cell_contours(img, cell_contour, out_file_name)