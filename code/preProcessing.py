import os
import config
import rasterio
import glob
from osgeo import gdal


def pre_deal(preprocessing_folder=config.preProcessing_folder, output_file=config.origin_multi_tif):
    print("预处理文件夹路径:", preprocessing_folder)  # 打印文件夹路径

    tiff_files = sorted(glob.glob(os.path.join(preprocessing_folder, '*.img')) +
                        glob.glob(os.path.join(preprocessing_folder, '*.tif')))

    print("找到的文件列表:", tiff_files)  # 打印找到的文件列表



    if not tiff_files:
        print("没有找到任何 IMG 或 TIF 文件进行预处理。不会进行多波段合成")
 
    if not os.path.exists(output_file):
        # 创建多波段 TIFF 文件
        with rasterio.open(tiff_files[0]) as src:
            profile = src.profile.copy()
            profile.update(count=len(tiff_files))

        with rasterio.open(output_file, 'w', **profile) as dst:
            for index, tiff_file in enumerate(tiff_files, start=1):
                with rasterio.open(tiff_file) as src:
                    dst.write(src.read(1), index)

        print(f'多波段 TIFF 文件已保存为 {output_file}')

    # 裁剪 TIFF 文件
    clip_raster_with_shapefile(output_file, config.origin_multi_tif_cut, config.border_shp)

def clip_raster_with_shapefile(input_raster_path, output_raster_path, shapefile_path):
    print(f"开始裁剪操作，输入文件: {input_raster_path}, 输出文件: {output_raster_path}, 形状文件: {shapefile_path}")

    warp_options = gdal.WarpOptions(cutlineDSName=shapefile_path, cropToCutline=True, dstNodata=-9999)
    result = gdal.Warp(output_raster_path, input_raster_path, options=warp_options)

    if result:
        print(f'裁剪后的 TIFF 文件已保存为 {output_raster_path}')
    else:
        print("裁剪操作失败")

if __name__ == "__main__":
    pre_deal()
