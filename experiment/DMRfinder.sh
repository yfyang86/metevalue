number=8
DMRfinderd=tools/DMRfinder
data=simulation_data/bed
output_dir=output_data/DMRfinder
filename=test${number}_DMRfinder.cov
difference=0.05

if [ ! -f "$filename" ]; then
    for((index=1;index<=${number};index++));
    do	
        cp ${data}/test${index}.bed ${data}/test${index}_DMRfinder.cov
        cp ${data}/control${index}.bed ${data}/control${index}_DMRfinder.cov
    done 
fi 

 
python3 ${DMRfinderd}/combine_CpG_sites.py -v  -o ${output_dir}/combined.${number}.csv `for((i=1;i<=${number};i++));do echo -ne "${data}/test${i}_DMRfinder.cov ";done``for((i=1;i<=${number};i++));do echo -ne "${data}/control${i}_DMRfinder.cov ";done`

Rscript ${DMRfinderd}/findDMRs.r -d ${difference} -i ${output_dir}/combined.${number}.csv -o ${output_dir}/DMRfinder_DMR -v -n Control,Exptl test1_DMRfinder,`for((i=2;i<=${number};i++));do echo -ne ",test${i}_DMRfinder";done`  control1_DMRfinder,`for((i=2;i<=${number};i++));do echo -ne ",control${i}_DMRfinder";done`


