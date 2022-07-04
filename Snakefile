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
        "src/data/RateData/Rate_info.h5"
    cache:
        True
    script:
        "src/scripts/CosmicIntegration/LinkHdf5Files.py"

