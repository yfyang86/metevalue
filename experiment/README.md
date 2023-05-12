This is about the pipeline of data process of 4 different tools.

First, in order to combine the bed data, you should run combine_bed.sh first.

Please use bash command to run toolsname.sh to get the DMRs calculated by different tools.
You can change the sample size, the at least difference of methylation rates, or the directory of data files by changing the parameters at the head of .sh files.

Then, you can run evalue_calculate_toolsname.sh to calculate the evalue of DMRs of different tools.

Finally, please run analysis_toolsname.sh to get the evaluation indicator of different methods and different tools.
