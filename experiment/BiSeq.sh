number=8
BiSeqd=tools/BiSeq
difference=0.05
data=simulation_data/bed
output_dir=output_data/BiSeq/BiSeq_DMR

Rscript ${BiSeqd}/BiSeq.R ${number} ${difference} ${data} ${output_dir}



