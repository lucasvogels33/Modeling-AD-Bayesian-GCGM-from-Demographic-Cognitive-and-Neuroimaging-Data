
#this R script does three things
#1) it loads the cleaned data
#2) it runs the Bayesian GCGM on these data
#3) it produces probability density plots for any selected variable pair


################################################
########1. load data########################
################################################

library(BDgraph)
setwd("C:/Users/lvogels/OneDrive - UvA/PhD-Project-Lucas/1_Projects/5. JAD/2. Data/Github")
filename = "Cleaned_data.Rdata"
load(file=filename)

#######################################################
###############2. run algorithm #######################
#######################################################

#define supporting functions
data_select_func = function(data,diagnosis_vec){
  
  #select subset of data according to diagnosis_vec
  index_diagnosis = which(data$Diagnosis %in% diagnosis_vec)
  data_1 = data[index_diagnosis,]
  
  #remove the diagnosis column
  data_1$Diagnosis = NULL #remove diagnosis column

  #output data
  return(data_1)
}
Gaussianize_func = function(data,not.cont,method){
  cont = 1 - not.cont
  index_cont = which(cont==1)
  data_cont = data[index_cont]
  if (method=="shrinkage"){data_cont_norm = bdgraph.npn(data=data_cont,npn="shrinkage")}
  if (method=="truncation"){data_cont_norm = bdgraph.npn(data=data_cont,npn="truncation")}
  if (method=="scaling"){data_cont_norm = scale(data_cont,center=TRUE,scale=TRUE)}
  data_1 = data
  data_1[index_cont] = data_cont_norm
  return(data_1)
}

#select subset of the data, with only the diagnosis that you want included, either 1,2,3,4 or a combination
diagnosis_vec = c(1,2,3,4) 
title = paste0("Diagnosis ",paste(diagnosis_vec,collapse=" "))
print(title)
data_select = data_select_func(data=data_5,diagnosis_vec=diagnosis_vec)

#normalize continuous data
not.cont=c(0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
matrix(c(colnames(data_select),not.cont),ncol=2) #check if non-continuous vector is correct (age assumed contninuos)
method = "truncation" 
data_norm = Gaussianize_func(data=data_select,not.cont=not.cont,method=method)

#select the burnin and MCMC iterations, and run the algorithm
burnin = 20000
iter = 120000

#run algorithm
sample_bd  = bdgraph( data = data_norm,algorithm = "bdmcmc",method="gcgm",not.cont=not.cont, iter = iter, burnin = burnin, jump = 1, save = TRUE,cores=1,g.start="empty",g.prior=0.2,verbose=TRUE)  

#save output
diag_text = paste(diagnosis_vec,collapse=" ")
filename = paste("Copula_output_saveTRUE_diag_",diag_text,".Rdata",sep="")
save(sample_bd,file=filename)


#######################################################################
############3. make probability plots #################################
#######################################################################

#supporting functions
prob_dens_plot = function(i=2,j=1,step_size=0.1,sample = sample_bd,hist=TRUE,max_ylim=0.75){
  
  #check if paramaters are of correct form 
  if (j >= i){ print("Please make sure j is smaller than i")}
  intervals = seq(-1,1,step_size)
  if (length(which(intervals==0))<1){
    print("interval vec needs to contain 0")
  }
  
  #save output of algorithm
  sample_graphs = sample$sample_graphs
  sample_K = sample$sample_K
  all_graphs = sample$all_graphs
  all_weights = sample$all_weights
  iter_tot = length(all_weights)
  
  #obtain series of partial correlation for a given connection
  ij_entry = 0.5*i*(i-1)+j
  diag_i_entry = 0.5*(i+1)*i
  diag_j_entry = 0.5*(j+1)*j
  C_ij_vec = rep(0,iter-burnin)
  for (iter_nr in (1:iter_tot)){
    
    #calculate precision matrix
    K_char = as.character(sample_K[iter_nr]) 
    K_vec = as.numeric(substring(K_char, seq(1, nchar(K_char), 5), seq(5, nchar(K_char), 5))) #split the character vector into strings of length 5
    K_ij = K_vec[ij_entry]
    K_ii = K_vec[diag_i_entry]
    K_jj = K_vec[diag_j_entry]
    
    #calculate partial correlations
    C_ij = -K_ij/sqrt(K_ii*K_jj)
    C_ij_vec[iter_nr] = C_ij
  }
  
  #create density vectors for the density plot
  len = length(intervals)
  weight = rep(0,len-1)
  mids = rep(0,len-1)
  for (q in 1:(len-1)){
    lower = intervals[q]
    upper = intervals[q+1]
    mids[q] = (lower+upper)/2
    if (upper<=0){
      indices = which(lower<=C_ij_vec & C_ij_vec<upper)
    }
    if (lower>=0){
      indices = which(lower<C_ij_vec & C_ij_vec<=upper)
    }
    weight[q] = sum(all_weights[indices])
  }
  density_vec = weight/sum(all_weights)
  
  #add the probabiltiy density for c_ij = 0
  weight_zero = sum(all_weights[which(C_ij_vec==0)])/sum(all_weights)
  mids = append(mids,0,after=length(mids)/2)
  density_vec = append(density_vec,weight_zero,after=length(density_vec)/2)
  
  var1 = colnames(sample$K_hat)[i]
  var2 = colnames(sample$K_hat)[j]
  if (hist==FALSE){
    #plot the density plot
    plot(NA,xlab="Partial correlation",ylab="Density",xlim=c(-1,1),ylim=c(0,max(density_vec)),main=title)
    points(x=mids,y=density_vec,type="l")
    
    #add a legend with mean, median, and sd
    mean_text = paste0("mean: ",round(mean(C_ij_vec),2))
    median_text = paste0("median: ",round(median(C_ij_vec),2))
    sd_text = paste0("sd: ",round(sd(C_ij_vec),2))
    legend("topleft",legend=c(mean_text,median_text,sd_text))
  }
  if (hist==TRUE){
    new_vec = c()
    for (i in 1:length(density_vec)){
      new_vec = c(new_vec,rep(mids[i],round(1000*density_vec[i])))
    }
    intervals_neg = seq(-1,-step_size,step_size)
    intervals_pos = seq(step_size,1,step_size)
    intervals = c(intervals_neg,-0.001,0.001,intervals_pos)
    out_zero = hist(new_vec,breaks=intervals,plot=FALSE)
    out = hist(new_vec,breaks=intervals,freq=TRUE,yaxt="n",xaxt="n",xlim=c(-0.5,0.5),ylab="density",xlab="partial correlation",ylim=c(0,sum(out_zero$counts)*max_ylim),col="azure3",border="azure3",main="")
    axis(2, at = seq(0,sum(out$counts),length=11),labels = seq(0,1,0.1))
    axis(1, at = seq(-0.5,0.5,0.1),labels = seq(-0.5,0.5,0.1))
    
    #add a legend with mean, median, and sd
    mean_text = paste0("mean: ",round(mean(C_ij_vec),2))
    median_text = paste0("median: ",round(median(C_ij_vec),2))
    sd_text = paste0("sd: ",round(sd(C_ij_vec),2))
    edge_prob = 100*(1-weight_zero)
    edge_prob_text = paste0("edge prob: ",round(edge_prob,0),"%")
    legend(x = c(-0.5, -0.1), y = c(0.75*sum(out$counts), 0.55*sum(out$counts)),legend=c(mean_text,median_text,sd_text,edge_prob_text))
  }
  
  return(list(mids=mids,density_vec=density_vec,C_ij_vec=C_ij_vec))
}

#load data
setwd("C:/Users/lvogels/OneDrive - UvA/PhD-Project-Lucas/1_Projects/5. JAD/2. Data/Github")
diagnosis_vec = c(1,2,3,4)
diag_text = paste(diagnosis_vec,collapse=" ")
filename = paste("Copula_output_saveTRUE_diag_",diag_text,".Rdata",sep="")
load(file=filename)
col_names = colnames(sample_bd$K_hat)

#make probability density plots for two selected variables
var1 = "Educ" #select the first variable
var2 = "Executive" #select the second variable
index_1 = which(col_names==var1)
index_2 = which(col_names==var2)
j = min(index_1,index_2)
i = max(index_1,index_2)
step_size = 0.05 #set the size of each bucket
max_ylim = 0.75 #set the limit of the y-axis
out_prob_dens = prob_dens_plot(i=i,j=j,step_size=step_size,sample = sample_bd,hist=TRUE,max_ylim=max_ylim)



