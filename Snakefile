# rule extract:
#     priority: 10
#     output:
#         "src/data/COMPAS_Output_wWeights.h5"
#     script:
#         "src/scripts/ExtractZenodoData.py"

rule CosmicIntegration:
    priority: 5
    output:
        "src/data/RateData/CI_job_IDs.txt"
    script:
        "src/scripts/CosmicIntegration/CallCosmicIntegration.py"

rule combineFiles:
    priority: 1
    output:
        "src/data/RateData/small_Rate_info.h5"
    cache:
        True
    script:
        "src/scripts/CosmicIntegration/CheckCompletionAndCombine.py"


