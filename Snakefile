rule extract:
    output:
        "src/data/COMPAS_Output_wWeights.h5"
    script:
        "src/scripts/ExtractZenodoData.py"
    priority: 3

rule CosmicIntegration:
    output:
        "src/data/RateData/CI_job_IDs.txt"
    script:
        "src/scripts/CosmicIntegration/CallCosmicIntegration.py"
    priority: 2

rule combineFiles:
    output:
        "src/data/RateData/Rate_info.h5"
    cache:
        True
    script:
        "src/scripts/CosmicIntegration/CheckCompletionAndCombine.py"
    priority: 1


