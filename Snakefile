rule extract:
    output:
        "src/data/COMPAS_Output_wWeights.h5"
    script:
        "src/scripts/ExtractZenodoData.py"


rule CosmicIntegration:
    output:
        "src/data/RateData/"
    cache:
        True
    script:
        "src/scripts/CosmicIntegration/CallCosmicIntegration.py"

rule LinkHdf5Files:
    output:
        "src/data/RateData/small_Rate_info.h5"
    cache:
        True
    script:
        "src/scripts/CosmicIntegration/LinkHdf5Files.py"

        # "src/data/1_small_Rate_info.h5"
        # "src/data/2_small_Rate_info.h5"
        # "src/data/3_small_Rate_info.h5"
        # "src/data/4_small_Rate_info.h5"
        # "src/data/5_small_Rate_info.h5"
        # "src/data/6_small_Rate_info.h5"
        # "src/data/7_small_Rate_info.h5"
        # "src/data/8_small_Rate_info.h5"
        # "src/data/9_small_Rate_info.h5"
        # "src/data/10_small_Rate_info.h5"
        # "src/data/11_small_Rate_info.h5"
        # "src/data/12_small_Rate_info.h5"
        # "src/data/13_small_Rate_info.h5"
        # "src/data/small_Rate_info.h5"