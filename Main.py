import config
import Processing
import preProcessing
import os



def main():
    
    if os.listdir(config.preProcessing_floder) and not os.path.exists(config.origin_multi_tif):
        print("preProcessing")
        preProcessing.pre_deal()
        
    print("image cut")   
    Processing.clip_raster_with_shapefile_updated(config.origin_multi_tif,config.origin_multi_tif_cut)
    print("cal index")
    Processing.calculate_indices()
    print("index cut")
    Processing.clip_raster_with_shapefile_updated()
    print("classification")
    Processing.plant_classify()

    print("All operations completed successfully.")

if __name__ == '__main__':
    main()
