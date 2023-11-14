import config
import rasterio
import numpy as np
from osgeo import gdal, ogr
import geopandas as gpd
from scipy.ndimage import binary_opening, binary_closing
from skimage.morphology import disk, dilation
from skimage.filters import threshold_multiotsu
import os
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
        
        # MNDWI
        mndwi = (B3 - B11) / (B3 + B11)
        with rasterio.open(f"{output_folder}/MNDWI.tif", 'w', **profile) as dst:
            dst.write(mndwi.astype(rasterio.float32), 1)
        #TC2 and svsi 
        brightness = 0.3561 * B2 + 0.3972 * B3 + 0.3904 * B4 + 0.6966 * B8 + 0.2286 * B11 + 0.1596 * B12
        greenness = -0.3344 * B2 - 0.3544 * B3 - 0.4556 * B4 + 0.6966 * B8 - 0.0242 * B11 - 0.2630 * B12
        svsi = brightness - greenness
        with rasterio.open(f"{output_folder}/SVSI.tif", 'w', **profile) as dst:
            dst.write(svsi.astype(rasterio.float32), 1)
        with rasterio.open(f"{output_folder}/TC2.tif", 'w', **profile) as dst:
            dst.write(greenness.astype(rasterio.float32), 1)
        # EVSI
        evsi = (B8a - B4) / (B8a + B4)
        with rasterio.open(f"{output_folder}/EVSI.tif", 'w', **profile) as dst:
            dst.write(evsi.astype(rasterio.float32), 1)
        profile.update(count=4)
        with rasterio.open(f"{output_folder}/totalIndex.tif", 'w', **profile) as dst:
            dst.write(np.stack([mndwi, greenness, evsi, svsi], axis=0).astype(rasterio.float32))

        


def clip_raster_with_shapefile(input_raster_path, output_raster_path, shapefile_path, nodata=-9999):
    # 栅格数据裁剪
    shapefile = ogr.Open(shapefile_path)
    warp_options = gdal.WarpOptions(cutlineDSName=shapefile_path, cropToCutline=True, dstNodata=nodata)
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

    # Apply decision tree classification
    #[mndwi,TC2,evsi,svsi];1为浮叶，2为挺水，3为沉水，4为水华，5为水体
    classification = np.full_like(scaled_data[0], 5, dtype=np.uint8)  # Initialize with 'Water'

    # 根据新的分类逻辑进行分类
    mask_mndwi_le_a = scaled_data[0] <= thresholds[0]
    mask_evsi_gt_c = scaled_data[2] > thresholds[2]
    classification[mask_mndwi_le_a & mask_evsi_gt_c] = 2  # 挺水植被

    mask_evsi_le_c = scaled_data[2] <= thresholds[2]
    classification[mask_mndwi_le_a & mask_evsi_le_c] = 1  # 浮叶植被

    mask_mndwi_gt_a = scaled_data[0] > thresholds[0]
    mask_tc2_le_b = scaled_data[1] <= thresholds[1]
    classification[mask_mndwi_gt_a & mask_tc2_le_b] = 5  # 水体

    mask_tc2_gt_b = scaled_data[1] > thresholds[1]
    mask_svi_lt_d = scaled_data[3] < thresholds[3]
    classification[mask_mndwi_gt_a & mask_tc2_gt_b & mask_svi_lt_d] = 3  # 沉水植被

    mask_svi_ge_d = scaled_data[3] >= thresholds[3]
    classification[mask_mndwi_gt_a & mask_tc2_gt_b & mask_svi_ge_d] = 4  # 水华

    classification[np.isnan(data[0])] = 255  # Set nodata values


    # Output the classified TIFF file
    output_tif = os.path.join(output_directory, 'result.tif')
    profile.update(dtype=rasterio.uint8, count=1, nodata=255)  # Update the data type to uint8, band count to 1, and nodata value to 255
    with rasterio.open(output_tif, 'w', **profile) as dst:
        dst.write(classification, 1)
        
    output_tif2 = os.path.join(output_directory, 'result_mask.tif')
    original_classification = classification.copy()

    for class_value in range(1, 6):  # 假设有 5 类植被
        mask = classification == class_value
        opened = binary_opening(mask, structure=disk(3))
        closed = binary_closing(opened, structure=disk(3))
        classification[mask] = 0
        classification[closed] = class_value

    # 使用膨胀操作填充未分类的像素
    dilated_classification = dilation(classification)
    classification[classification == 0] = dilated_classification[classification == 0]

    # 保存分类结果
    profile.update(dtype=rasterio.uint8, count=1, nodata=255)
    with rasterio.open(output_tif2, 'w', **profile) as dst:
        dst.write(classification, 1)

    print("分类完成，结果保存在:", output_tif2)
    

def tc2_process():
    try:
        print("TC2 - 计算指数")
        calculate_indices(config.origin_multi_tif_cut, config.tc2_index_total_folder)

        print("TC2 - 影像裁剪")
        clip_raster_with_shapefile(config.tc2_index_total_file, config.tc2_index_total_file_cut, config.border_shp)

        print("TC2 - 植被分类")
        plant_classify(config.tc2_index_total_file_cut, config.border_shp, config.tc2_result_folder)

        print("TC2 - 处理完成")
        

    except Exception as e:
        print(f"TC2 处理过程中出错: {e}")
    
    try:
        raster_to_shapefile(config.tc2_result_folder + '/result.tif', config.tc2_result_folder + '/result.shp')
        raster_to_shapefile(config.tc2_result_folder + '/result_mask.tif', config.tc2_result_folder + '/result_mask.shp')
        print("Shapefile 保存成功")
    except Exception as e:
        print(f"Shapefile 转换过程中出错: {e}")

if __name__ == "__main__":
    tc2_process()