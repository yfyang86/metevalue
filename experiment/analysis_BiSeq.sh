analysis_type=BH #p or bonferroni or BH or e
tools_type=BiSeq
evalue_dir=output_data/BiSeq/BiSeq_evalue
DMR_dir=simulation_data/bed/DMRs_unDMRs_signal.bed
output_dir=output_data/BiSeq/BiSeq_analysis_${analysis_type}

Rscript test.R ${analysis_type} ${tools_type} ${evalue_dir} ${DMR_dir} ${output_dir}



