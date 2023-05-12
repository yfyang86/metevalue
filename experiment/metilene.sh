data=simulation_data/bed/rate_combine_8
output_dir=output_data/metilene/metilene_DMR
difference=0.05
metilened=tools/metilene_v0.2-8

${metilened}/metilene -a g1 -b g2 -d ${difference} ${data} > ${output_dir}

