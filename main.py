# imports for band gap extraction
import vision
import bandextractor
import compextractorb
import os
import requests
# imports for visualization
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

# Download data from OSF database
bandextractor.download(url='https://osf.io/download/aebuq',filename='FAMAPbI.hdr')
bandextractor.download(url='https://osf.io/download/27h3e',filename='FAMAPbI.bil')
bil = 'data/FAMAPbI.bil'
hdr = 'data/FAMAPbI.hdr'

# Read motor speed and gcode files used for composition extraction 
datapath = os.getcwd() 
with open(datapath+'/data/motor_speeds.txt', 'r') as file:
    motor_speeds = file.read()
gcode = pd.read_csv(os.getcwd()+'/data/gcode_XY.csv')

# Variables
rotate_crop_params = {
    'theta': 1, # Image rotation angle
    'x1': 55, # left x bound 
    'x2': 870, # Right x bound 
    'y1': 270, # Top y bound 
    'y2': 720 # Bottom y bound 
} # USER DEFINED ROTATE/CROP PARAMETERS

compextractor_params = {
    'rasterSpeed' : 38, # mm/s gocde raster speed
    'rasterScaleY' : 0.8, # scale gcode raster pattern in Y-dim
    'rasterScaleX' : 0.95, # scale gcode raster pattern in X-dim
    'rasterOffsetY' : 30, # offset gcode raster pattern in X-dim
    'rasterOffsetX' : -10 # offset gcode raster pattern in X-dim
}


# Segment Droplets
data, sgmnt = vision.segmentation(bil=bil, hdr=hdr, rotate_crop_params=rotate_crop_params, savepath=f'{datapath}/Bandgap/', return_segment = True) 

# Extract Compositions 
idx, comp, vals = compextractorb.get_compositions(data, segmentation=sgmnt, motor_speeds=motor_speeds, gcode=gcode, rasterSpeed=compextractor_params['rasterSpeed'], rasterScaleY=compextractor_params['rasterScaleY'], rasterScaleX=compextractor_params['rasterScaleX'], rasterOffsetX=compextractor_params['rasterOffsetX'], rasterOffsetY=compextractor_params['rasterOffsetY'], savepath=f'{datapath}/Bandgap/')

# run band gap extractor
egs = bandextractor.autoextract_direct(data=vals, savepath=f'{datapath}/Bandgap/', verbose=False)

# Create x and Y Vectors to Plot the Band Gaps 
x, y = bandextractor.mapping(comp,egs)

# plot the extracted band gaps
bandextractor.plot_ev(x,y,datapath)

# save sorted raw spectra as csv
bandextractor.save_csv(comp,vals,datapath)