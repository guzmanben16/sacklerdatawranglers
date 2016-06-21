#!/usr/bin/env python

from skimage import io, filter, measure
import mahotas
from scipy import ndimage
import numpy as np

#this script will perform some simple image processing for Super Resolution Microscopy Images
# (1) read in the image (tif)
# (2) separate the colors - 'magenta' and 'green'
# (3) threshold the image - produces image masks
# (4) detect clusters
# (5) measure clusters
# (6) output measured clusters data into .csv file (can be opened as a table in excel)
# (7) a little bit more advanced - colocalization analysis

pixels_to_nm = 20 # (1 px = 20 nm for reconstructed images)
image_file = 'image_for_basic_image_analysis.tif'

# (1) read in the image file - use scikit-image io.imread function
image_rgb_array = io.imread(image_file)

# (2) separate the colors
print 'Image array shape: ' + str(image_rgb_array.shape)
print 'Image array dtype: ' + str(image_rgb_array.dtype)
red_image = image_rgb_array[:,:,0] # take the red color only
green_image = image_rgb_array[:,:,1] #take the green color only

#save the red and green images as individual greayscale images
io.imsave('image_red.tif', red_image)
io.imsave('image_green.tif', green_image)

# (3) threshold the image, global methods return the threshold value:
#Ridler-Calvard method, this is default ImageJ 'threshold': 
# http://fiji.sc/Auto_Threshold#Default
th = mahotas.rc(red_image)
red_mask = red_image > th

th  = mahotas.rc(green_image)
green_mask = green_image > th

print 'RC threshold: ' + str(th)

#save red and green mask - white out True values and cast the data type
io.imsave('image_red_mask.tif', 255*red_mask.astype('uint8')) #, plugin='freeimage')
io.imsave('image_green_mask.tif', 255*green_mask.astype('uint8')) #, plugin='freeimage')

# (4) detect the clusters - cluster regions are labeled by integers, 1 to num_labels
#use 8-connectivity, diagonal pixels will be included as part of a structure
#this is ImageJ default but we have to specify this for Python, or 4-connectivity will be used
s = [[1,1,1],
	[1,1,1],
	[1,1,1]]

red_labeled_clusters, num_labels = ndimage.label(red_mask, structure=s)
green_labeled_clusters, num_labels = ndimage.label(green_mask, structure=s)

# (5) measure cluster properties
propList = ['Area',
            'Image', #sliced binary region image same size as bounding box
            'BoundingBox',
            'MajorAxisLength',
            'MinorAxisLength',
            'Perimeter',
            'MinIntensity',
            'MeanIntensity',
            'MaxIntensity']

red_props = measure.regionprops(red_labeled_clusters, red_image) #send in original image for Intensity measurements
green_props = measure.regionprops(green_labeled_clusters, green_image) #send in original image for Intensity measurements

# create image with an outline of the labeled clusters (red and green), with overlap area in white
outline_image = np.zeros((len(red_image), len(red_image[0]), 3), dtype='uint8')

# (6) output data to file
#remove BBox and Image, we are not printing those to the csv file
propList.remove('Image')
propList.remove('BoundingBox')

# get data for the red clusters
for cluster_props in red_props[0:2]:
    print str(cluster_props['Label'])
    for i,prop in enumerate(propList):
        if(prop == 'Area'):
            to_print = cluster_props[prop]*pixels_to_nm**2
        elif(prop.find('Intensity') < 0):
            to_print = cluster_props[prop]*pixels_to_nm
        else:
            to_print = cluster_props[prop]
        print '\t' + prop + ':\t' + str(to_print) + '\n'



