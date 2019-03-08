# basic_descriptors
===


# Create a conda environment :

Install Anaconda following the website instruction: https://www.anaconda.com/distribution/

`conda create -n dypfish python=2.7`

(This creates a dypfish environment for python)
Activate : 

`conda activate dypfish`

# Install packages

`pip install matplotlib`
  
  `conda install numpy`
  
  `conda install scipy`
  
  `pip install scikit-image`
  
 `pip install h5py`
  
  `conda install os`
  
  `conda install csv`

# Download this repository :
`git clone https://github.com/ebouilhol/basic_descriptors.git`

`cd basic_descriptors`

Open the "run_create_3d_images.py" file, change the path value according to the localisation of your image folder. 

This script creates stack of images for tubulin, nucleus, FISH based on CY3/CY5/DAPI images, simply launch it : 

`python run_create_3d_images.py`

After that, we are going to generate the HDF5 file based on the previous stacked images. 
Open the "main.py" file, change the path value according to the localisation of your folder. 

This script generates cell mask, nucleus mask, nucleus position, spots position. 
It does not provide MTOC position. You need to provide it in a CSV file, and put it in the mtoc/ folder. To create the CSV file, take a look at the example provided. It should be named as the folder. 

Launch the main : 

`python main.py`


It should create a basic_descriptors.h5 file. You may get some warnings.

If you want to visualize the results, uncomment the function plot_nucleus_contours in dapi_all_slices (utils.py)
And plot_cell_contours in tubulin_all_slices (utils.py)

The scripts needs to be adjusted for each series of images.

You can adjust the values in the main.py file

    tubulin_case = 0
    tubulin_treshold = 78

    dapi_case = 0
    dapi_treshold = 80

Tubulin case if 0 : try an auto treshold. If the result is not good, switch the case to 1 and try different treshold values until the results are good. 
Same thing for dapi case and dapi treshold. 



