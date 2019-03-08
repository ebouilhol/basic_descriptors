#!/usr/bin/python
import skimage.io as io
import os
from utils import *
import h5py
import warnings

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()


if __name__ == '__main__':
    path = './29.4.2016_RACK1siRNA-BactinFISH'
    dirs = os.listdir(path)

    h5_file = h5py.File("basic_descriptors.h5", 'w')

    tubulin_case = 0
    tubulin_treshold = 78

    dapi_case = 0
    dapi_treshold = 80

    for root, subdirs, files in os.walk(path):
        if subdirs == []:
            print("Running : "+root)
            grp = h5_file.create_group(root)
            tubulinfile= os.path.join(root, "tubulin.tif")
            tubulin_out_file_name =os.path.join(root, "tubulinfile.png")
            tubulin_3d = io.imread(tubulinfile)
            cell_out_file_name = os.path.join(root, "cell_plot")
            cell_mask_output_name = os.path.join(root, "cell_mask")

            tubulin_all_slices(root, tubulin_3d, 0, 0, 30, cell_mask_output_name, cell_out_file_name, 15, 15, 23, 10, 0.1, 0.95, grp, path)


            dapifile = os.path.join(root, "dapi.tif")
            dapi_out_file_name =os.path.join(root, "dapi_plot")
            dapi_mask_output_name =os.path.join(root, "nucleus_mask")
            dapi_3d = io.imread(dapifile)
            dapi_all_slices(dapi_3d, dapi_case, dapi_treshold, dapi_mask_output_name, dapi_out_file_name, 8, 10, 8, 10, 0.1, 0.95, grp)

