rule extract:
    output:
        "src/data/COMPAS_Output_wWeights.h5"
    script:
        "src/scripts/ExtractZenodoData.py"

rule CosmicIntegration:
    output:
        "src/data/RateData/CI_job_IDs.txt"
    cache:
        True
    script:
        "src/scripts/CosmicIntegration/CallCosmicIntegration.py"

rule combineFiles:
    input: 
        "src/data/RateData/CI_job_IDs.txt"
    output:
        "src/data/RateData/Rate_info.h5"
    cache:
        True
    script:
        "src/scripts/CosmicIntegration/CheckCompletionAndCombine.py"


