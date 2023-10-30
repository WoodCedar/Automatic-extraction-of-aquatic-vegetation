import os


#工作文件夹、湖泊边界
work_floder = r"F:\新建文件夹"
border_shp =work_floder+r"\边界.shp"


#预处理wenjianjia
preProcessing_floder = work_floder+r"\preProcessing"
if not os.path.exists(preProcessing_floder):
    os.makedirs(preProcessing_floder)

#原始多波段文件
origin_multi_tif =work_floder+r"\totalTif.tif"
origin_multi_tif_cut = work_floder+r"\totalTif_cut.tif"

#old index和结果
old_index_total_floder = work_floder+r"\oldIndex"
if not os.path.exists(old_index_total_floder):
    os.makedirs(old_index_total_floder)
old_index_total_file = old_index_total_floder+r"\totalIndex.tif"
old_index_total_file_cut = old_index_total_floder+r"\totalIndexCut.tif"




#index floder
index_total_floder = work_floder+r"\index"
if not os.path.exists(index_total_floder):
    os.makedirs(index_total_floder)
index_total_file = index_total_floder+r"\totalIndex.tif"
index_total_file_cut = index_total_floder+r"\totalIndexCut.tif"

#resutl floder
result_floder = work_floder+r"\result"
if not os.path.exists(result_floder):
    os.makedirs(result_floder)


old_result_floder = work_floder+r"\oldresult"
if not os.path.exists(old_result_floder):
    os.makedirs(old_result_floder)
old_result_file = old_result_floder+r"\result.tif"
