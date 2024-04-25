# GATK-gCNV-polishing-bash-pipeline

This script has been developed at the Institute of Biological Chemistry of the School of Exact and Natural Science (IQUIBICEN)
of University of Buenos Aires.

This script contains several steps for CNV priorization after gGATK-CNVs pipelines for cohort and case samples.

The first step involves: 
  _ * Identifying and remvoing Common CNVs in the VCF output of the pipeline (with cnv-filter.sh)
  _ *  After Common CNVs filtering,each CNV interval is extended by 50% length at each side (objetive is merged close segments split by GATK's Viterby algorithm)
  _ *  The segments that overlapped, after the extension, are merged using bedtools merge tool.
  _ *  Then, a dynamic QS is defined based on Copy Number (CN) and Number of Exons (NP) for each CNV
  _ *  then if the estimated CNV QS is lower that defined threshold for each condition, the CNV is filtered out, just keeping high-quality CNVs for further analysis.
