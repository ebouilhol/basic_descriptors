#!/usr/bin/python

import argparse
import os
# import skimage
from skimage import io, exposure, img_as_uint, img_as_float
import numpy as np
import glob

parser = argparse.ArgumentParser(description='Create 3d image file from slices')
parser.add_argument('-t','--type', help='IF or FISH', required=True)
parser.add_argument('-i','--in', help='Input path', required=True)
parser.add_argument('-o','--out', help='Output path', required=True)

args = vars(parser.parse_args())
imageType = args['type']
out = args['out']
path = args['in']

extensioncy3 = "CY3*.tif"
extensioncy5 = "CY5*.tif"
extensiondapi = "DAPI*.tif"
extensionYFP = "YFP*.tif"

def stack_function(path, ext, out, imageName):
    # for dirs in os.listdir(path):
        # print(dirs)
    path2 = os.path.join(path) + "/*" + ext
        # print(path2)
    num_slice = 0
    image_3d = []
    if glob.glob(path2) !=[]:
        # print("glob glob !")
        for filename in np.sort(glob.glob(path2)):
            # print(filename)
            slice_temp = io.imread(filename)
            image_3d = np.append(image_3d, slice_temp)
            num_slice += 1
        image_3d = image_3d.reshape(num_slice, 512, 512)
        image_3d = exposure.rescale_intensity(image_3d, out_range=(0, 255))
        image_3d = image_3d.astype(int)
        image_3d = img_as_uint(image_3d)

        if not os.path.exists(os.path.join(out)):
            os.makedirs(os.path.join(out))
        io.imsave(os.path.join(out,imageName+".tif"), image_3d)

print("Process FISH or IF...")
stack_function(path, extensioncy3, out, imageType)

print("Process Tubulin...")
stack_function(path, extensioncy5, out, "tubulin")

print("Process Tubulin...")
stack_function(path, extensionYFP, out, "micropattern")

print("Process DAPI...")
stack_function(path, extensiondapi, out, "dapi")