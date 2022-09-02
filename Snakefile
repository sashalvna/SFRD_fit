rule fit_data:
    priority: 10
    input: "src/data/SFRMetallicityFromGasTNG100.hdf5"
    output:
        "src/data/test_best_fit_parameters.txt"
    script:
        "src/scripts/Fit_model_to_sfrdzZ.py"

rule extract:
    priority: 10
    input:"src/data/COMPAS_Output_wWeights.tar.gz"
    output:
        "src/data/COMPAS_Output_wWeights.h5"
    script:
        "src/scripts/ExtractZenodoData.py"

rule extractFig5:
    input:
        "src/data/Figure5.tar.gz"
    output:
        directory("src/data/Figure5")
    script:
        "src/scripts/Extract_fig5data.py"


rule CosmicIntegration:
    priority: 5
    input: 
        "src/data/test_best_fit_parameters.txt"
    output:
        "src/data/RateData/CI_job_IDs.txt"
    script:
        "src/scripts/CosmicIntegration/CallCosmicIntegration.py"


rule combineFiles:
    priority: 1
    input: "src/data/RateData/CI_job_IDs.txt"
    output:
        "src/data/RateData/Rate_info.h5"
    cache:
        True
    script:
        "src/scripts/CosmicIntegration/CheckCompletionAndCombine.py"


