rule extract:
    output:
        "src/data/COMPAS_Output_wWeights.h5"
    script:
        "src/scripts/ExtractZenodoData.py"


rule CosmicIntegration:
    output:
        "src/data/RateData/CI_job_IDs.txt"
    script:
        "src/scripts/CosmicIntegration/CallCosmicIntegration.py"

rule combineFiles:
    output:
        "src/data/RateData/Rate_info.h5"
    cache:
        True
    script:
        "src/scripts/CosmicIntegration/CheckCompletionAndCombine.py"


# rule LinkHdf5Files:
#     output:
        # "src/data/RateData/"
#         "src/data/RateData/Rate_info.h5"
#     cache:
#         True
#     script:
#         "src/scripts/CosmicIntegration/LinkHdf5Files.py"

