#a_b is the data like the input type of metilene
#use chr, start, end to determine the region should be calculated

library(stringr)
evalue_normal = function(a_b, chr, start, end){
  innerf = function(x, innermu=0., innersig=1.){
    vector_temp = na.omit(as.numeric(x))
    n = length(vector_temp)
    value = mean(vector_temp)
    dnorm(x=value, mean = innermu, sd = innersig/sqrt(n))
  }
  
  innerlog = function(x){
    log(x[!is.na(x)])
  }
  
  site_1 = grep('g1',names(a_b),value=T)
  site_2 = grep('g2',names(a_b),value=T)
  
  start_temp = start
  end_temp = end
  chr_temp = chr

  print(chr_temp)
  a_b_temp = a_b[(a_b$pos>=start_temp & a_b$pos<=end_temp & a_b$'chr'==chr_temp), ]

  vector_1 = unlist(a_b_temp[, site_1])
  vector_2 = unlist(a_b_temp[, site_2])
  vector = c(vector_1, vector_2)
  miu_1 = mean(vector_1, na.rm=T)
  sigma_1 = sd(vector_1, na.rm=T)
  miu_2 = mean(vector_2, na.rm=T)
  sigma_2 = sd(vector_2, na.rm=T)
  miu = mean(vector, na.rm=T)
  sigma = sd(vector, na.rm=T)
    
  norm_value_up_1 = apply(a_b_temp[, site_1], 2, function(x)innerf(x=x,innermu=miu_1, innersig=sigma_1))
  norm_value_down_1 = apply(a_b_temp[, site_1], 2, function(x)innerf(x=x, innermu=miu, innersig=sigma))
  norm_value_up_2 = apply(a_b_temp[, site_2], 2,  function(x)innerf(x=x, innermu=miu_2, innersig=sigma_2))
  norm_value_down_2 = apply(a_b_temp[, site_2], 2,  function(x)innerf(x=x, innermu=miu, innersig=sigma))

  e_value = exp(sum(innerlog(c(norm_value_up_1,norm_value_up_2))) -
                    sum(innerlog(c(norm_value_down_1,norm_value_down_2))))
  return(e_value)
}


