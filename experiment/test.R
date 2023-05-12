arg = commandArgs(T)
analysis_type = arg[1]
tools_type = arg[2]
evalue_dir = arg[3]
DMR_dir = arg[4]
output_dir = arg[5]
#analysis_type='e'
#tools_type='BiSeq'
#evalue_dir='output_data/BiSeq/BiSeq_evalue'
#DMR_dir='simulation_data/bed/DMRs_unDMRs_signal.bed'
library(regioneR)
library(dplyr)
library(pROC)
out1=read.csv(evalue_dir,header=T,sep='\t')
if(analysis_type=='BH' | analysis_type=='bonferroni'){
  if(tools_type=='BiSeq'){
    out1$median.p = p.adjust(out1$median.p, method=analysis_type)
  }
  else{
    out1$p = p.adjust(out1$p, method=analysis_type)
  }
}
out2=read.table(DMR_dir,header=F,sep='\t')
names(out2)=c("chr2","start","end","isDMR")
wendy=overlapRegions(out1,out2)
print("over")
wendy=cbind(wendy,0)
colnames(wendy)[7]="covlen"
for(i in 1:nrow(wendy))
{
if(wendy[i,"type"]=="ArightB")
wendy[i,"covlen"]=wendy[i,"endB"]-wendy[i,"startA"]
else if(wendy[i,"type"]=="AleftB")
wendy[i,"covlen"]=wendy[i,"endA"]-wendy[i,"startB"]
else if(wendy[i,"type"]=="AinB")
wendy[i,"covlen"]=wendy[i,"endA"]-wendy[i,"startA"]
else if(wendy[i,"type"]=="BinA")
wendy[i,"covlen"]=wendy[i,"endB"]-wendy[i,"startB"]
}
if(analysis_type=='e'){
  out1_1=cbind(out1$start,out1$end,out1$e_value)
  out1_1=data.frame(out1_1)
  names(out1_1)=c("start","end","e_mean")
  wendy_all=left_join(wendy,out1_1,c("startA"="start","endA"="end"))
  wendy_all=left_join(wendy_all,out2,c("startB"="start","endB"="end"))
  wendy_all=cbind(wendy_all,qxz=0)
  wendy_all[wendy_all$e_mean>20,"qxz"]=1
}
if(analysis_type=='p' | analysis_type=='bonferroni' | analysis_type=='BH'){
  if(tools_type=='BiSeq'){
    out1_1=cbind(out1$start,out1$end,out1$median.p)
  }
  else{
    out1_1=cbind(out1$start,out1$end,out1$p)
  }
  out1_1=data.frame(out1_1)
  names(out1_1)=c("start","end","e_mean")
  wendy_all=left_join(wendy,out1_1,c("startA"="start","endA"="end"))
  wendy_all=left_join(wendy_all,out2,c("startB"="start","endB"="end"))
  wendy_all=cbind(wendy_all,qxz=0)
  wendy_all[wendy_all$e_mean<0.05,"qxz"]=1  
}
wendy_group=group_by(wendy_all,startB)
wendy_group=summarise(wendy_group,count=n())
wendy_group=data.frame(wendy_group)
wendy_group=cbind(wendy_group,0,0)
names(wendy_group)=c("startB","count","isFind","covrate")
for(i in 1:nrow(wendy_group))
{
Bvalue=wendy_group[i,"startB"]
wendy_temp=wendy_all[wendy_all$startB==Bvalue,]
sum_temp=sum(wendy_temp$covlen)
Blen_temp=wendy_temp[1,"endB"]-wendy_temp[1,"startB"]
q_temp=0
for(j in 1:nrow(wendy_temp))
{
  q_temp=q_temp+wendy_temp[j,"qxz"]*(wendy_temp[j,"covlen"])
}
if(q_temp>Blen_temp/2)
  wendy_group[i,"isFind"]=1
wendy_group[i,"covrate"]=q_temp/Blen_temp
}
out2_all=left_join(out2,wendy_group,c("start"="startB"))
out2_all[is.na(out2_all$isFind),"isFind"]=0
out2_all[is.na(out2_all$covrate),"covrate"]=0
out2_group=group_by(out2_all,isDMR,isFind)
out2_group=summarise(out2_group,count=n())
out2_group=data.frame(out2_group)
beta=matrix(nrow=2,ncol=2)
out2_group_ex=data.frame(isDMR=c(0,1,0,1),isFind=c(0,0,1,1),count2=c(0,0,0,0))
out2_group=merge(out2_group_ex,out2_group,all.x=T)
out2_group[is.na(out2_group$count), "count"]=0
beta[1,1]=out2_group[out2_group$isDMR==1 & out2_group$isFind==1,"count"]
beta[2,1]=out2_group[out2_group$isDMR==1 & out2_group$isFind==0,"count"]
beta[1,2]=out2_group[out2_group$isDMR==0 & out2_group$isFind==1,"count"]
beta[2,2]=out2_group[out2_group$isDMR==0 & out2_group$isFind==0,"count"]
beta=data.frame(beta)
rownames(beta)=c("Y","N")
colnames(beta)=c("Y","N")
rownames(beta)=c("Y_m","N_m")

print("complete calculate")
############################################################################################################


#Analysis
data_ana = data.frame(1:1)
data_matrix = beta
data_auc = out2_all

ACC = (data_matrix[1,1]+data_matrix[2,2])/(data_matrix[1,1]+data_matrix[1,2]+data_matrix[2,1]+data_matrix[2,2])
FDR = (data_matrix[1,2])/(data_matrix[1,2]+data_matrix[1,1])
FPR = (data_matrix[1,2])/(data_matrix[1,2]+data_matrix[2,2])
FNR = (data_matrix[2,1])/(data_matrix[2,1]+data_matrix[1,1])
roc1 = roc(data_auc$isDMR,data_auc$covrate,levels=c(0,1))
AUC = as.numeric(auc(roc1))

data_ana = cbind(data_ana,ACC,FDR,FPR,1-FNR,AUC)
names(data_ana) = c('a','ACC','FDR','Type_I_error','Power','AUC')
data_ana = data_ana[,-1]

write.table(data_ana,output_dir,row.names=F,col.names=T,quote=F,sep='\t')
