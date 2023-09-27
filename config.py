import os


#work floder and border
work_floder = r"*****"
border_shp =work_floder+r"\border.shp"


#preProcessing floder
preProcessing_floder = work_floder+r"\preProcessing"


#original image file
origin_multi_tif =work_floder+r"\totalTif.tif"
origin_multi_tif_cut = work_floder+r"\totalTif_cut.tif"



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

