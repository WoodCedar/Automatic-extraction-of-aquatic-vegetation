import config
import ndviProcessing
import tc2Processing
import preProcessing
import os

def process_program(process_func, program_name):
    """
    通用程序处理函数
    :param process_func: 处理函数
    :param program_name: 程序名称
    if os.listdir(config.preProcessing_folder) and not os.path.exists(config.origin_multi_tif) and not os.path.exists(config.origin_multi_tif_cut):
        print(f"{program_name} - 预处理")
    """
    
    preProcessing.pre_deal()
    print(f"{program_name} - 执行程序")
    process_func()

def ndvi_program():
    # NDVI 程序逻辑
    ndviProcessing.ndvi_process()

def tc2_program():
    # TC2 程序逻辑
    tc2Processing.tc2_process()

# 执行 NDVI 程序
print("\n开始进行ndvi分类\n")
process_program(ndvi_program, "NDVI Program")
print("\n开始进行tc2分类\n")
# 执行 TC2 程序
process_program(tc2_program, "TC2 Program")
