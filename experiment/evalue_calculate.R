library(metevalue)
arg = commandArgs(T)
tools = arg[1]
methy_rate_dir = arg[2]
tools_output_dir = arg[3]
output_dir = arg[4]

if(tools=='DMRfinder'){
    data_evalue =  metevalue.DMRfinder(methy_rate_dir, tools_output_dir, bheader=TRUE)
    write.table(data_evalue, output_dir, row.names=F, col.names=T, quote=F, sep='\t')
}

if(tools=='methylKit'){
    data_evalue =  metevalue.methylKit(methy_rate_dir, tools_output_dir, bheader=TRUE)
    write.table(data_evalue, output_dir, row.names=F, col.names=T, quote=F, sep='\t')
}

if(tools=='BiSeq'){
    data_evalue =  metevalue.biseq(methy_rate_dir, tools_output_dir, bheader=TRUE)
    write.table(data_evalue, output_dir, row.names=F, col.names=T, quote=F, sep='\t')
}

if(tools=='metilene'){
    data_evalue =  metevalue.metilene(methy_rate_dir, tools_output_dir)
    write.table(data_evalue, output_dir, row.names=F, col.names=T, quote=F, sep='\t')
}