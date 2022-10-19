library(stringr)
evalue_inn = function(x, c1, c2){
  miu1 = mean(x[c1])
  miu2 = mean(x[c2])
  sigma1 = sd(x[c1])
  sigma2 = sd(x[c2])
  miu = mean(x[c(c1,c2)])
  sigma = sd(x[c(c1,c2)])
  norm_value_1 = dnorm(x[c1],miu1,sigma1)/dnorm(x[c1],miu,sigma)
  norm_value_2 = dnorm(x[c2],miu2,sigma2)/dnorm(x[c2],miu,sigma)
  miu
}

evalue_norm = function(data){
  name = names(data)
  pos_1 = str_count(name, 'g1')
  pos_2 = str_count(name, 'g2')
  site_1 = name[pos_1==1]
  site_2 = name[pos_2==1]
  if(length(site_1)<2 | length(site_2)<2)
    warning("Format error!!! ")
  norm_evalue = apply(data, 1, evalue_inn, site_1, site_2)
  norm_evalue
}

