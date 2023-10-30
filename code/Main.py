from . import config
from . import Processing
from . import preProcessing
from . import oldProcessing


import os



def tc2Processing():
    
    if os.listdir(config.preProcessing_floder) and not os.path.exists(config.origin_multi_tif):
        print("预处理")
        preProcessing.pre_deal()
        
    print("影像裁剪")   
    Processing.clip_raster_with_shapefile_updated(config.origin_multi_tif,config.origin_multi_tif_cut)
    print("指数计算")
    Processing.calculate_indices()
    print("指数裁剪")
    Processing.clip_raster_with_shapefile_updated()
    print("分类")
    Processing.plant_classify()

    print("All operations completed successfully.")
    
    
    
    
    


def ndviProcessing():
     
    if os.listdir(config.preProcessing_floder) and not os.path.exists(config.origin_multi_tif) and os.path.exists(config.preProcessing_floder):
        print("预处理")
        preProcessing.pre_deal()
        
    print("影像裁剪")   
    oldProcessing.clip_raster_with_shapefile_updated(config.origin_multi_tif,config.origin_multi_tif_cut)
    print("指数计算")
    oldProcessing.calculate_indices()
    print("指数裁剪")
    oldProcessing.clip_raster_with_shapefile_updated()
    print("分类")
    oldProcessing.plant_classify()
    



#tc2Processing() #shaoyifan tc2
#ndviOldProcessing() #zhanghaobin ndvi\evsi