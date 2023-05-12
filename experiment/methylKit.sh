number=8
methylKitd=tools/methylKit
difference=0.05
data=simulation_data/bed
output_dir=output_data/methylKit/methylKit_DMR

Rscript ${methylKitd}/methylKit.R ${number} ${difference} ${data} ${output_dir}



