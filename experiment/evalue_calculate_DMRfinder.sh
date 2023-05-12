tools=DMRfinder
methy_rate_dir=simulation_data/bed/rate_combine_8
tools_output_dir=output_data/DMRfinder/DMRfinder_DMR
output_dir=output_data/DMRfinder/DMRfinder_evalue

Rscript evalue_calculate.R ${tools} ${methy_rate_dir} ${tools_output_dir} ${output_dir}

