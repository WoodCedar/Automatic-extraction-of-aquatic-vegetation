import config
import rasterio
import numpy as np
from skimage.filters import threshold_multiotsu
from osgeo import gdal, ogr
import os
import geopandas as gpd
from scipy.ndimage import binary_opening, binary_closing
from skimage.morphology import disk, dilation
from rasterio.features import shapes
from shapely.geometry import shape

def calculate_indices(input_tif, output_folder):
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
        # FAI 
        fai = B8 - (B4 + 177 * ((B11 - B4) / 945))
        with rasterio.open(f"{output_folder}/FAI.tif", 'w', **profile) as dst:
            dst.write(fai.astype(rasterio.float32), 1)
        profile.update(count=5)
        with rasterio.open(f"{output_folder}/totalIndex.tif", 'w', **profile) as dst:
            dst.write(np.stack([evsi,fai,mndwi,ndvi,svsi], axis=0).astype(rasterio.float32))

        
def clip_raster_with_shapefile(input_raster_path, output_raster_path, shapefile_path, nodata=-9999):
    shapefile = ogr.Open(shapefile_path)
    warp_options = gdal.WarpOptions(cutlineDSName=shapefile_path, cropToCutline=True, dstNodata=nodata, xRes=10, yRes=10, resampleAlg=gdal.GRIORA_NearestNeighbour)
    gdal.Warp(output_raster_path, input_raster_path, options=warp_options)
def raster_to_shapefile(raster_path, output_shp_path):
    """
    将栅格数据转换为 Shapefile 文件。
    :param raster_path: 输入的栅格数据（TIFF）文件路径。
    :param output_shp_path: 输出的 Shapefile 文件路径。
    """
    with rasterio.open(raster_path) as src:
        image = src.read(1)  # 假设分类结果在第一个波段
        results = ({'properties': {'raster_val': v}, 'geometry': s}
                    for i, (s, v) in enumerate(shapes(image, transform=src.transform)))

        gdf = gpd.GeoDataFrame.from_features(list(results))
        gdf.crs = src.crs
        gdf.to_file(output_shp_path)

def plant_classify(input_tif, input_shp, output_directory):
    gdf = gpd.read_file(input_shp)
    with rasterio.open(input_tif) as src:
        data = src.read().astype(np.float32) 
        profile = src.profile

    data[data == -9999] = np.nan
    scaled_data = np.empty_like(data, dtype=np.uint8)
    for i in range(data.shape[0]):
        scaled_data[i] = ((data[i] - np.nanmin(data[i])) / (np.nanmax(data[i]) - np.nanmin(data[i])) * 255).astype(np.uint8)
    
    nbins_values = [min(1024, len(np.unique(band))) for band in scaled_data]
    thresholds = [threshold_multiotsu(band, nbins=n)[-1] for band, n in zip(scaled_data, nbins_values)]
    thresholds[1] = threshold_multiotsu(scaled_data[1], nbins=nbins_values[1])[-2]

    output_txt = os.path.join(output_directory, 'thresholds.txt')
    with open(output_txt, 'w') as txt_file:
        txt_file.write(f'NDVI Thresholds\n')
        for i, threshold in enumerate(thresholds):
            txt_file.write(f'Band {i+1}: {threshold}\n')

    classification = np.full_like(scaled_data[0], 5, dtype=np.uint8)
    classification[np.isnan(data[0])] = 5
    classification[(scaled_data[0] > thresholds[0])] = 2  # Emergent vegetation
    classification[(scaled_data[0] <= thresholds[0]) & (scaled_data[2] > thresholds[2])] =5  # Water
    classification[(scaled_data[0] <= thresholds[0]) & (scaled_data[2] <= thresholds[2]) & (scaled_data[1] > thresholds[1])] =4#水华
    classification[(scaled_data[0] <= thresholds[0]) & (scaled_data[1] <= thresholds[1]) & (scaled_data[2] <= thresholds[2]) & (scaled_data[3] > thresholds[3])] = 1  # Floating vegetation
    classification[(scaled_data[0] <= thresholds[0]) & (scaled_data[1] <= thresholds[1]) &(scaled_data[2] <= thresholds[2]) & (scaled_data[3] <= thresholds[3]) & (scaled_data[4] < thresholds[4])] = 3 # Submerged vegetation


    classification[np.isnan(data[0])] = 255  # Set nodata values


    output_tif = os.path.join(output_directory, 'result.tif')
    profile.update(dtype=rasterio.uint8, count=1, nodata=255)
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


def ndvi_process():
    try:
        print("NDVI - 计算指数")
        calculate_indices(config.origin_multi_tif_cut, config.ndvi_index_total_folder)

        print("NDVI - 影像裁剪")
        clip_raster_with_shapefile(config.ndvi_index_total_file, config.ndvi_index_total_file_cut, config.border_shp)

        print("NDVI - 植被分类")
        plant_classify(config.ndvi_index_total_file_cut, config.border_shp, config.ndvi_result_folder)

        print("NDVI - 处理完成")
    except Exception as e:
        print(f"NDVI 处理过程中出错: {e}")
    
    try:
        raster_to_shapefile(config.ndvi_result_folder + '/result.tif', config.ndvi_result_folder + '/result.shp')
        raster_to_shapefile(config.ndvi_result_folder + '/result_mask.tif', config.ndvi_result_folder + '/result_mask.shp')
        print("Shapefile 保存成功")
    except Exception as e:
        print(f"Shapefile 转换过程中出错: {e}")

if __name__ == "__main__":
    ndvi_process()