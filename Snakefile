rule extract:
    output:
        "src/data/COMPAS_Output_wWeights.h5"
    script:
        "src/scripts/ExtractZenodoData.py"


rule CosmicIntegration:
    output:
        "src/data/Rate_info.hdf5"
    cache:
        True
    script:
        "src/scripts/CosmicIntegration/CallCosmicIntegration.py"