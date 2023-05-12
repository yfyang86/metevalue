analysis_type=e #p or bonferroni or BH or e
tools_type=metilene
evalue_dir=output_data/metilene/metilene_evalue
DMR_dir=simulation_data/bed/DMRs_unDMRs_signal.bed
output_dir=output_data/metilene/metilene_analysis_${analysis_type}

Rscript test.R ${analysis_type} ${tools_type} ${evalue_dir} ${DMR_dir} ${output_dir}



