import os

# 工作文件夹路径 - 可以根据需要修改
work_folder = r"E:\工作\2023提取工作\23年数据\草海230705完成\20230705"

# 边界文件路径
border_shp = os.path.join(work_folder, "边界.shp")

# 预处理文件夹路径
preProcessing_folder = os.path.join(work_folder, "preProcessing")
if not os.path.exists(preProcessing_folder):
    os.makedirs(preProcessing_folder)

# 原始多波段文件路径
origin_multi_tif = os.path.join(work_folder, "totalTif.tif")
origin_multi_tif_cut = os.path.join(work_folder, "totalTif_cut.tif")

# NDVI 指数文件和结果文件夹路径
ndvi_index_total_folder = os.path.join(work_folder, "ndviIndex")
if not os.path.exists(ndvi_index_total_folder):
    os.makedirs(ndvi_index_total_folder)
ndvi_index_total_file = os.path.join(ndvi_index_total_folder, "totalIndex.tif")
ndvi_index_total_file_cut = os.path.join(ndvi_index_total_folder, "totalIndexCut.tif")

# TC2 指数文件夹路径
tc2_index_total_folder = os.path.join(work_folder, "tc2Index")
if not os.path.exists(tc2_index_total_folder):
    os.makedirs(tc2_index_total_folder)
tc2_index_total_file = os.path.join(tc2_index_total_folder, "totalIndex.tif")
tc2_index_total_file_cut = os.path.join(tc2_index_total_folder, "totalIndexCut.tif")

# 结果文件夹路径
tc2_result_folder = os.path.join(work_folder, "tc2result")
if not os.path.exists(tc2_result_folder):
    os.makedirs(tc2_result_folder)
tc2_result_file = os.path.join(tc2_result_folder, "tc2result.tif")

ndvi_result_folder = os.path.join(work_folder, "ndviresult")
if not os.path.exists(ndvi_result_folder):
    os.makedirs(ndvi_result_folder)
ndvi_result_file = os.path.join(ndvi_result_folder, "ndviresult.tif")
