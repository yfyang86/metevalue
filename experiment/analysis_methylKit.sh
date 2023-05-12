analysis_type=p #p or bonferroni or BH or e
tools_type=methylKit
evalue_dir=output_data/methylKit/methylKit_evalue
DMR_dir=simulation_data/bed/DMRs_unDMRs_signal.bed
output_dir=output_data/methylKit/methylKit_analysis_${analysis_type}

Rscript test.R ${analysis_type} ${tools_type} ${evalue_dir} ${DMR_dir} ${output_dir}



