
import os
import config 
import rasterio
from rasterio.features import shapes
from osgeo import gdal, ogr, osr
import numpy as np
import glob
from skimage.morphology import dilation,disk
from skimage.filters import threshold_multiotsu
from skimage import exposure
import geopandas as gpd
from scipy.ndimage import binary_opening, binary_closing


def calculate_indices(input_tif=config.origin_multi_tif_cut, output_folder=config.index_total_floder):
    
    with rasterio.open(input_tif) as src:
        profile = src.profile
        profile.update(dtype=rasterio.float32, count=1, compress='lzw')
        B2 = src.read(2).astype(np.float32)
        B3 = src.read(3).astype(np.float32)
        B4 = src.read(4).astype(np.float32)
        B8 = src.read(8).astype(np.float32)
        B8a = src.read(9).astype(np.float32)
        B11 = src.read(11).astype(np.float32)
        B12 = src.read(12).astype(np.float32)
        
        #NDVI
        ndvi = (B8 - B4) / (B8 + B4)
        profile.update(nodata=None)
        with rasterio.open(f"{output_folder}/NDVI.tif", 'w', **profile) as dst:
            dst.write(ndvi.astype(rasterio.float32), 1) 
            
        # MNDWI
        mndwi = (B3 - B11) / (B3 + B11)
        with rasterio.open(f"{output_folder}/MNDWI.tif", 'w', **profile) as dst:
            dst.write(mndwi.astype(rasterio.float32), 1)
            
        #TC2 
        brightness = 0.3561 * B2 + 0.3972 * B3 + 0.3904 * B4 + 0.6966 * B8 + 0.2286 * B11 + 0.1596 * B12
        greenness = -0.3344 * B2 - 0.3544 * B3 - 0.4556 * B4 + 0.6966 * B8 - 0.0242 * B11 - 0.2630 * B12
        
        svsi = brightness - greenness
        with rasterio.open(f"{output_folder}/SVSI.tif", 'w', **profile) as dst:
            dst.write(svsi.astype(rasterio.float32), 1)

        # EVSI
        evsi = (B8a - B4) / (B8a + B4)
        with rasterio.open(f"{output_folder}/EVSI.tif", 'w', **profile) as dst:
            dst.write(evsi.astype(rasterio.float32), 1)

        # FAI 从177可改为40
        fai = B8 - (B4 + 177 * ((B11 - B4) / 945))
        with rasterio.open(f"{output_folder}/FAI.tif", 'w', **profile) as dst:
            dst.write(fai.astype(rasterio.float32), 1)

        profile.update(count=5)
        with rasterio.open(f"{output_folder}/totalIndex.tif", 'w', **profile) as dst:
            dst.write(np.stack([evsi,fai,mndwi,ndvi,svsi], axis=0).astype(rasterio.float32))

        print(f"result save as {output_folder}")
        

def clip_raster_with_shapefile_updated(input_raster_path=config.index_total_file, output_raster_path=config.index_total_file_cut, shapefile_path=config.border_shp, nodata=-9999):
    # Open the shapefile
    shapefile = ogr.Open(shapefile_path)
    layer = shapefile.GetLayer()
    
    # Use gdal warp to align the output raster with the shapefile boundary
    warp_options = gdal.WarpOptions(cutlineDSName=shapefile_path, 
                                    cropToCutline=True, 
                                    dstNodata=nodata, 
                                    xRes=10, yRes=10,  # assuming 10x10 resolution, modify as per your data
                                    resampleAlg=gdal.GRIORA_NearestNeighbour)  # Nearest neighbour resampling, modify if needed
    
    # Warp (clip and align) the input raster using the shapefile
    gdal.Warp(output_raster_path, input_raster_path, options=warp_options)

    print("clip finish.")


def plant_classify(input_tif=config.index_total_file_cut,input_shp=config.border_shp,output_directory=config.result_floder):

    # Read the shapefile
    gdf = gpd.read_file(input_shp)
    # Read the TIFF file
    with rasterio.open(input_tif) as src:
        data = src.read().astype(np.float32) 
        profile = src.profile
    data[data == -9999] = np.nan

    # Scale each band to 0-255 and convert to uint8
    scaled_data = np.empty_like(data, dtype=np.uint8)
    for i in range(data.shape[0]):
        scaled_data[i] = ((data[i] - np.nanmin(data[i])) / (np.nanmax(data[i]) - np.nanmin(data[i])) * 255).astype(np.uint8)

    nbins_values = [min(1024, len(np.unique(band))) for band in scaled_data]

    thresholds = [threshold_multiotsu(band, nbins=n)[-1] for band, n in zip(scaled_data, nbins_values)]

    thresholds[1] = threshold_multiotsu(scaled_data[1], nbins=nbins_values[1])[-2]

    output_txt = os.path.join(output_directory, 'thresholds.txt')
    with open(output_txt, 'w') as txt_file:
        txt_file.write(f'EVSI Threshold: {thresholds[0]}\n')
        txt_file.write(f'FAI Threshold: {thresholds[1]}\n')
        txt_file.write(f'MNDWI Threshold: {thresholds[2]}\n')
        txt_file.write(f'NDVI Threshold: {thresholds[3]}\n')
        txt_file.write(f'SVSI Threshold: {thresholds[4]}\n')

    # Apply decision tree classification
    #[mndwi,TC2,evsi,svsi];
    classification = np.full_like(scaled_data[0], 5, dtype=np.uint8)  # Initialize with 'Water'
    classification[np.isnan(data[0])] = 5  # Set to water

    classification[(scaled_data[0] > thresholds[0])] = 2  # Emergent vegetation
    classification[(scaled_data[0] <= thresholds[0]) & (scaled_data[2] > thresholds[2])] =5  # Water
    classification[(scaled_data[0] <= thresholds[0]) & (scaled_data[2] <= thresholds[2]) & (scaled_data[1] > thresholds[1])] =4#水华
    classification[(scaled_data[0] <= thresholds[0]) & (scaled_data[1] <= thresholds[1]) & (scaled_data[2] <= thresholds[2]) & (scaled_data[3] > thresholds[3])] = 1  # Floating vegetation
    classification[(scaled_data[0] <= thresholds[0]) & (scaled_data[1] <= thresholds[1]) &(scaled_data[2] <= thresholds[2]) & (scaled_data[3] <= thresholds[3]) & (scaled_data[4] < thresholds[4])] = 3 # Submerged vegetation


    classification[np.isnan(data[0])] = 255  # Set nodata values


    # Output the classified TIFF file
    output_tif = os.path.join(output_directory, 'result.tif')
    profile.update(dtype=rasterio.uint8, count=1, nodata=255)  # Update the data type to uint8, band count to 1, and nodata value to 255
    with rasterio.open(output_tif, 'w', **profile) as dst:
        dst.write(classification, 1)
        
    output_tif2 = os.path.join(output_directory, 'result_mask.tif')
    original_classification = classification.copy()

    for class_value in [1, 2, 3, 4, 5]:
        mask = classification == class_value
        opened = binary_opening(mask, structure=disk(3))
        closed = binary_closing(opened, structure=disk(3))
        classification[mask] = 0
        classification[closed] = class_value

 
    dilated_classification = dilation(classification)
    classification[classification == 0] = dilated_classification[classification == 0]


    with rasterio.open(output_tif2, 'w', **profile) as dst:
        dst.write(classification, 1)

    print("classification finish")
    
