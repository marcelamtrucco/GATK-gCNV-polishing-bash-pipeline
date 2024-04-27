# GATK-gCNV-polishing-bash-pipeline

This script has been developed at the Institute of Biological Chemistry of the School of Exact and Natural Science (IQUIBICEN)
of University of Buenos Aires.

This script contains several steps for CNV priorization after gGATK-CNVs pipelines for cohort and case samples.

The first step involves: 
  * Identifying and remvoing Common CNVs in the VCF output of the pipeline (with cnv-filter.sh)
  *  After Common CNVs filtering, each CNV interval is extended by 50% length at both sides (objetive is merge close segments that were split by GATK's Viterby algorithm)
  *  Overlapping segments, after the extension, are merged using bedtools merge tool.
  *  Then, a dynamic QS is defined based on the Copy Number (CN) and Number of Exons (NP) for each CNV. If the estimated QS is lower than the defined threshold
    for each condition, the CNV is filtered out, just keeping high-quality CNVs for further analysis.
    
