#功能：将sentinel-2的L2a img数据转换为tif格式，并将12个波段合成一个
# Compare this snippet from 2023%E6%8F%90%E5%8F%96/%E9%95%BF%E6%B1%9F%E6%B0%B4%E7%94%9F%E6%80%81/1_img2tif.py:
import os
from . import config
from osgeo import gdal
import numpy as np
from rasterio.merge import merge
import rasterio
import glob

def pre_deal(preProcessing_floder=config.preProcessing_floder,output_file=config.origin_multi_tif):

    # Change this to the directory containing your input single-band TIFFs
    input_directory = preProcessing_floder

    output_directory = config.work_floder

    output_file = os.path.join(output_directory, 'totalTif.tif')

    tiff_files = sorted(glob.glob(os.path.join(input_directory, '*.img')) + 
                   glob.glob(os.path.join(input_directory, '*.tif')))


    # Open the first TIFF file to get the profile
    with rasterio.open(tiff_files[0]) as src:
        profile = src.profile.copy()

    # Update the profile for the output file
    profile.update(count=len(tiff_files))

    # Create the output multi-band TIFF
    with rasterio.open(output_file, 'w', **profile) as dst:
        for index, tiff_file in enumerate(tiff_files, start=1):
            with rasterio.open(tiff_file) as src:
                band_data = src.read(1)
                dst.write(band_data, index)

    print(f'Multi-band TIF saved as {output_file}')
