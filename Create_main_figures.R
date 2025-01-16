
#This file does three things
#1. It loads the cleaned data
#2. It runs the algorithm
#3. It creates all the figures in the article and the figures S2 and S3 in the supplementary material

####################################################
############# libraries ##########################
####################################################

#library(readxl) #for loading excel file to R
library(circlize)
library(gplots)
library(ggplot2)
library(BDgraph)
library(RColorBrewer) #for heatmaps
library(psych)
library(gridExtra)


################################################
########1. load data########################
################################################

library(BDgraph)
setwd("C:/Users/lvogels/OneDrive - UvA/PhD-Project-Lucas/1_Projects/5. JAD/2. Data/Github")
filename = "Cleaned_data.Rdata"
load(file=filename)

#######################################################
############2. run algorithms###########################
#######################################################

#define supporting functions
BD_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart,g.prior,not.cont){
  
  #obtain p and n
  p = ncol(data)
  n = nrow(data)
  
  #parameters
  burnin = 0
  jump = 1
  save = FALSE
  cores = 1
  verbose = FALSE
  
  #run the algorithm, produce output at every iteration in the iter_vec_thin vector
  len_iter = length(iter_vec_thin)
  len_burnin = length(burnin_iter_vec_thin)
  edge_vec = rep(0,len_iter+len_burnin)
  time_vec = rep(0,len_iter+len_burnin)
  
  for (j in 1:len_burnin){
    print(burnin_iter_vec_thin[j])
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = burnin_iter_vec_thin[j-1]
    }
    newiter = burnin_iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_bd  = bdgraph( data = data,algorithm = "bdmcmc",method="gcgm",not.cont=not.cont,iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #save metrics
    edge_vec[j] = sum(sample_bd$last_graph)/2
    
    #save data for next run
    MCMCstart = sample_bd
    
  }
  
  #run MCMC iterations (after burnin)
  weights_old = 0
  plinks_old = 0
  K_old = 0
  plinks_iter_matrix = matrix(0,p*(p-1)/2,len_iter)
  K_iter_matrix = matrix(0,p*(p-1)/2+p,len_iter)
  for (j in 1:len_iter){
    
    print(iter_vec_thin[j])
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_bd  = bdgraph( data = data,algorithm = "bdmcmc",method="gcgm",not.cont=not.cont, iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #calculate new p matrix
    weights_new = sample_bd$sum_weights
    plinks_new = (weights_old*plinks_old + weights_new*sample_bd$p_links)/(weights_old+weights_new)
    K_new = (weights_old*K_old + weights_new*sample_bd$K_hat)/(weights_old+weights_new)
    
    #save metrics
    edge_vec[len_burnin+j] = sum(sample_bd$last_graph)/2
    
    #save data for next run
    plinks_old = plinks_new
    K_old = K_new
    weights_old = weights_old+weights_new
    MCMCstart = sample_bd
    
    #save plinks in matrix
    plinks_vec = plinks_new[upper.tri(plinks_new)]
    plinks_iter_matrix[,j] = plinks_vec
    
    #save K in matrix
    K_vec = K_new[upper.tri(K_new,diag=TRUE)]
    K_iter_matrix[,j] = K_vec
    
  }
  
  plinks = plinks_new
  return(list(g.prior=g.prior,edge_vec=edge_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new,plinks_iter_matrix=plinks_iter_matrix,K_iter_matrix=K_iter_matrix))
} 
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
run_algorithm_func = function(data,burnin_iter_vec_thin,iter_vec_thin,not.cont){
  MCMCstart = "empty"
  g.prior = 0.2
  out = BD_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart,g.prior=g.prior,not.cont=not.cont)
  return(out)
}

#select subset of the data, with only the diagnosis that you want included, either 1,2,3,4 or a combination
diagnosis_vec = c(1,2,3,4) 
title = paste0("Diagnosis ",paste(diagnosis_vec,collapse=" "))
print(title)
data_select = data_select_func(data=data_5,diagnosis_vec=diagnosis_vec)

#normalize continuous data
not.cont=c(0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
matrix(c(colnames(data_select),not.cont),ncol=2) #check if non-continuous vector is correct (age assumed contninuos)
method = "truncation" #select either shrinkage, truncation or scaling
data_norm = Gaussianize_func(data=data_select,not.cont=not.cont,method=method)

#select the burnin and MCMC iterations, and run the algorithm
burnin_iter_vec_thin = c(10,50,100,500,1000,2500,5000,10000,20000)
iter_vec_thin = c(100,200,500,1000,2500,5000,7500,10000,15000,20000,30000,40000,50000,60000,70000,80000,90000,100000)
output_alg = run_algorithm_func(data=data_norm,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,not.cont=not.cont)

#save output
diag_text = paste(diagnosis_vec,collapse=" ")
filename = paste("Copula_output_diag_",diag_text,".Rdata",sep="")
save(output_alg,file=filename)

#######################################################
############make figures ########################
#######################################################

#introduce functions
conv_stat_func = function(out,converge_metric,title=title){
  
  #analyze convergence statistics
  if (convergence_metric == "P"){
    convergence_iter_matrix = out$plinks_iter_matrix
    ylab = "Edge inclusion probabilities"
  }
  if (convergence_metric == "K"){
    convergence_iter_matrix = out$K_iter_matrix
    ylab = "Precision matrix entry"
  }
  
  iter_vec_thin = out$iter_vec_thin
  convergence_final = convergence_iter_matrix[,length(iter_vec_thin)]
  links_needed = 10
  
  #find the link with the highest edge inclusion prob
  max_convergence = max(convergence_final)
  min_convergence = min(convergence_final)
  intervals = seq(min_convergence,max_convergence,(max_convergence-min_convergence)/(links_needed))
  index_links = c()
  for (i in 1:links_needed){
    lower = intervals[i]
    upper = intervals[i+1]
    index = which(convergence_final >= lower & convergence_final < upper)[1]
    if (length(index)>0){ #if an index is found within this interval, then add it to the list
      index_links = c(index_links,index)
    }
  }
  
  #plot convergence
  max_iter = iter_vec_thin[length(iter_vec_thin)]
  plot(NA,xlim=c(1,max_iter),xlab="MCMC iterations after burnin",ylab=ylab,ylim=c(min_convergence,max_convergence),main=title)
  for (i in 1:links_needed){
    index = index_links[i]
    points(x=iter_vec_thin,y=convergence_iter_matrix[index,],type="l")
  }
}
create_P_K_C = function(output_alg){
  
  #load p and K_iter_matrix
  plinks = output_alg$plinks
  K_iter_matrix = output_alg$K_iter_matrix
  
  #create the symmetric plinks matrix
  p = dim(plinks)[1]
  count = 0
  for (i in 2:p){ #fill lower matrix
    for (j in 1:(i-1)){
      count = count + 1
      plinks[i,j] = plinks[j,i]
    } 
  }
  
  #create the symmetric precision matrix
  final_iter = dim(K_iter_matrix)[2]
  K_vec = K_iter_matrix[,final_iter]
  K = matrix(0,p,p)
  K[upper.tri(K,diag=TRUE)] = K_vec
  count = 0
  for (i in 2:p){ #fill lower matrix
    for (j in 1:(i-1)){
      count = count + 1
      K[i,j] = K[j,i]
    } 
  }
  
  #create the partial correlations matrix
  K_diag_vec = 1/sqrt(diag(K))
  K_diag = diag(K_diag_vec)
  C = -K_diag%*%K%*%K_diag
  
  #set variable names
  old_names = colnames(plinks)
  new_names = new_names = c("Age","Memory","Executive","Sex","Educ","APOE4","Amy-stage",
                            "G Thal","V Thal",
                            "G PCC","V PCC",
                            "G Prec","V Prec",
                            "G Hipp","V Hipp",
                            "G Caud","V Caud",
                            "G Put","V Put")      
  matrix(c(old_names,new_names),byrow=FALSE,ncol=2)
  colnames(plinks) = new_names
  rownames(plinks) = new_names
  colnames(K) = new_names
  rownames(K) = new_names
  colnames(C) = new_names
  rownames(C) = new_names
  
  #reorder columns and rows
  order_vec=c(8,10,12,14,16,18,9,11,13,15,17,19,1,4,5,6,7,2,3)
  var_names = colnames(plinks)
  matrix(c(var_names,order_vec),ncol=2) #check if new labels and ordering is correct
  plinks = plinks[order_vec,order_vec]
  K = K[order_vec,order_vec]
  C = C[order_vec,order_vec]
  
  #save output
  return(list(plinks=plinks,K=K,C=C))
  
}
make_circle = function(plinks,C,highlight,highlight_type,cutoff){
  
  #create matrix with three columns (name of var1, name of var2, edge incl. prob.)
  var_names = colnames(plinks)
  len = length(var_names)
  pair_matrix = matrix(0,nrow=len*(len-1)/2,ncol=3)
  
  #define colors
  col_vec = rep(0,len*(len-1)/2)
  lightblue = "#0000FF10"
  lightred = "#FF000010"
  darkblue = "#0000FFA9"
  darkred = "#FF0000A9"
  
  #define edges in need of highlighting
  edge_incl_sex_cognition = c("SexMemory",
                              "SexExecutive",
                              "SexEduc",
                              "EducMemory",
                              "EducExecutive",
                              "SexAmy-stage",
                              "Amy-stageMemory",
                              "Amy-stageExecutive",
                              "SexV Hipp",
                              "V HippMemory",
                              "V HippExecutive",
                              "SexV PCC",
                              "V PCCMemory",
                              "V PCCExecutive")
  edge_incl_age_cognition = c("AgeMemory",
                              "AgeExecutive",
                              "AgeAmy-stage",
                              "Amy-stageMemory",
                              "Amy-stageExecutive",
                              "AgeV Hipp",
                              "V HippMemory",
                              "V HippExecutive",
                              "AgeV PCC",
                              "V PCCMemory",
                              "V PCCExecutive")
  edge_incl_biomarker_cognition = c("V HippMemory",
                                    "G PCCMemory",
                                    "V PCCExecutive",
                                    "Amy-stageMemory",
                                    "Amy-stageExecutive",
                                    "APOE4Amy-stage")
  edge_incl_demographic_biomarker = c("SexV Hipp",
                                      "SexV Prec",
                                      "SexV PCC",
                                      "SexV Put",
                                      "SexV Caud",
                                      "SexV Thal",
                                      "AgeV Put",
                                      "AgeV Caud",
                                      "AgeV Hipp",
                                      "AgeV Prec",
                                      "AgeV PCC",
                                      "AgeV Thal",
                                      "SexV Hipp",
                                      "SexG Prec",
                                      "SexG PCC",
                                      "SexG Put",
                                      "SexG Caud",
                                      "SexG Thal",
                                      "AgeG Put",
                                      "AgeG Caud",
                                      "AgeG Hipp",
                                      "AgeG Prec",
                                      "AgeG PCC",
                                      "AgeG Thal")
  edge_incl_gluc_vol = c("V HippG Hipp",
                         "V HippG Prec",
                         "V HippG Caud",
                         "V HippG Thal",
                         "V HippG PCC",
                         "V HippG Put",
                         "V PrecG Hipp",
                         "V PrecG Prec",
                         "V PrecG Caud",
                         "V PrecG Thal",
                         "V PrecG PCC",
                         "V PrecG Put",
                         "V CaudG Hipp",
                         "V CaudG Prec",
                         "V CaudG Caud",
                         "V CaudG Thal",
                         "V CaudG PCC",
                         "V CaudG Put",
                        "V ThalG Hipp",
                         "V ThalG Prec",
                         "V ThalG Caud",
                         "V ThalG Thal",
                         "V ThalG PCC",
                         "V ThalG Put",
                         "V PCCG Hipp",
                          "V PCCG Prec",
                          "V PCCG Caud",
                          "V PCCG Thal",
                          "V PCCG PCC",
                            "V PCCG Put",
                        "V PutG Hipp",
                        "V PutG Prec",
                        "V PutG Caud",
                        "V PutG Thal",
                        "V PutG PCC",
                        "V PutG Put")
  
  if (highlight){
    if (highlight_type=="sex_cog"){edge_incl = edge_incl_sex_cognition}
    if (highlight_type=="age_cog"){edge_incl = edge_incl_age_cognition}
    if (highlight_type=="bio_cog"){edge_incl = edge_incl_biomarker_cognition}
    if (highlight_type=="dem_bio"){edge_incl = edge_incl_demographic_biomarker}
    if (highlight_type=="glu_vol"){edge_incl = edge_incl_gluc_vol}
  }
  
  #fill matrix with edge inclusion prob and fill col_vec with colours
  count = 0
  for (i in 1:(len-1)){
    for (j in (i+1):len){
      count = count + 1
      
      #fill first two colmuns with var names
      pair_matrix[count,1] = var_names[i]
      pair_matrix[count,2] = var_names[j]
      
      #fill third column with edge inclusion prob
      if (plinks[i,j] > cutoff){
        pair_matrix[count,3] = plinks[i,j]
      }
      
      #fill colour vector with red (if part_cor < 0) and blue (if part_cor > 0)
      if (highlight){
        edgename1 = paste0(var_names[i],var_names[j])
        edgename2 = paste0(var_names[j],var_names[i])
        if ((edgename1 %in% edge_incl) | (edgename2 %in% edge_incl)){
          col_vec[count] = darkred
          if (C[i,j] > 0){
            col_vec[count] = darkblue
          }  
        } else {
          col_vec[count] = lightred
          if (C[i,j] > 0){
            col_vec[count] = lightblue
          }
        }
      } else {
        col_vec[count] = darkred
        if (C[i,j] > 0){
          col_vec[count] = darkblue
        } 
      }
      
    }
  }
  
  #transform matrix to data frame
  pairs_df = as.data.frame(pair_matrix)
  colnames(pairs_df) = c("var1","var2","edge_incl_prob")
  pairs_df$edge_incl_prob = as.numeric(pairs_df$edge_incl_prob)
  
  #set colour of the bording
  grid.col = c(Age = "#CCCCCC",
               Memory = "#CCCCCC",
               Executive = "#CCCCCC",
               Sex = "#CCCCCC",
               Educ = "#CCCCCC",
               APOE4 = "#CCCCCC",
               "Amy-stage" = "#CCCCCC",
               "G Thal" = "#999999",
               "V Thal" = "#666666",
               "G PCC" = "#999999",
               "V PCC" = "#666666",
               "G Prec" = "#999999",
               "V Prec"= "#666666",
               "G Hipp" = "#999999",
               "V Hipp"= "#666666",
               "G Caud" = "#999999",
               "V Caud"= "#666666",
               "G Put" = "#999999",
               "V Put"   = "#666666")
  
  #set gaps between borders
  circos.par(gap.after = 
               c(Age = 1,
                 Memory = 1,
                 Executive = 5,
                 Sex = 1,
                 Educ = 5,
                 APOE4 = 1,
                 "Amy-stage" = 5,
                 "G Thal" = 1,
                 "V Thal" = 1,
                 "G PCC" = 1,
                 "V PCC" = 1,
                 "G Prec" = 1,
                 "V Prec"= 1,
                 "G Hipp" = 1,
                 "V Hipp"= 1,
                 "G Caud" = 1,
                 "V Caud"= 1,
                 "G Put" = 5,
                 "V Put"   = 5)
  )
  
  #intialize starting degree
  circos.par(start.degree = 20)
  
  #make circle
  chordDiagram(pairs_df,grid.col=grid.col,col=col_vec,transparency = 0,
               ,annotationTrack = "grid", 
               ,preAllocateTracks = list(track.height = max(strwidth(var_names))))
  
  #add labels to circle
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) 
  
  
  #delete all settings (for next chart)
  circos.clear()
  
  
}
heatmap_part_func = function(plinks,C,K,cutoff){
  
  #set elements (i,j) of C to zero if plinks(i,j) < 0.5
  index_exclude = which(plinks<cutoff)
  C[index_exclude] = 0
  diag(C)=1
  
  #heatmap using gplot
  lmat = rbind(c(0,3,0),c(0,1,2),c(0,4,0))
  lwid = c(0.25,3.5,0.25)
  lhei = c(0.25,3.5,1.25)
  margins = c(6,6)
  out_heat2 = heatmap.2(x=C,margins=margins,lmat=lmat,lwid=lwid,lhei=lhei,dendrogram="none",Rowv=FALSE,Colv=FALSE,scale="none",trace="none",density.info = "none",
                        col=brewer.pal(11,"RdBu"),breaks=seq(-1,1,length.out=12),
                        key.title=NA,key.xlab="partial correlation",key.par=list(cex=1),ylab=NA,
                        labRow = rownames(C),labCol = colnames(C),cexRow=1.25,cexCol=1.25)
  
}
heatmap_pearson_func = function(data){
  
  #run hypothesis tests for each pair of variables
  cor_out = corr.test(data,adjust="none")
  pvalues = cor_out$p
  cor_matrix = cor_out$r
  
  #set all correlations to zero that have pvalues bigger than 0.05
  index_not_significant = which(pvalues>0.05)
  cor_matrix[index_not_significant] = 0
  
  #heatmap using gplot
  lmat = rbind(c(0,3,0),c(0,1,2),c(0,4,0))
  lwid = c(0.25,3.5,0.25)
  lhei = c(0.25,3.5,1.25)
  margins = c(6,6)
  out_heat2 = heatmap.2(x=cor_matrix,margins=margins,lmat=lmat,lwid=lwid,lhei=lhei,dendrogram="none",Rowv=FALSE,Colv=FALSE,scale="none",trace="none",density.info = "none",
                        col=brewer.pal(11,"RdBu"),breaks=seq(-1,1,length.out=12),
                        key.title=NA,key.xlab="Pearson correlation",key.par=list(cex=1),ylab=NA,
                        labRow = rownames(cor_matrix),labCol = colnames(cor_matrix),cexRow=1.25,cexCol=1.25)
  
}

#load output
setwd("C:/Users/lvogels/OneDrive - UvA/PhD-Project-Lucas/1_Projects/5. JAD/2. Data/Github")
diagnosis_vec = c(1,2,3,4)
diag_text = paste(diagnosis_vec,collapse=" ")
filename = paste("Copula_output_diag_",diag_text,".Rdata",sep="")
load(file=filename)

#perform convergence analysis
title = paste0("Diagnosis ",paste(diagnosis_vec,collapse=" "))
convergence_metric = "P" #select either P for edge inclusion convergence, or K for precision matrix convergence
conv_stat_func(out=output_alg,converge_metric = converge_metric,title=title)

#create plinks, K and C and data (all with correct names and ordering of the rows and columns)
out = create_P_K_C(output_alg=output_alg)
plinks = out$plinks
K = out$K
C = out$C

#create circle plot
highlight=TRUE #select TRUE (FALSE) if you (do not) want to highlight certain edges 
highlight_type="glu_vol" #select either sex_cog, age_cog, bio_cog, dem_bio, glu_vol
par(mar = c(0,0,0,0))
make_circle(plinks=plinks,C=C,highlight=highlight,highlight_type=highlight_type,cutoff=0.5)

#create heatmap partial correlation
cutoff = 0.5
heatmap_part_func(plinks=plinks,C=C,K=K,cutoff=cutoff)

#create heatmap Pearson correlation
heatmap_pearson_func(data_select)

#Pearson statistics
cor_out = corr.test(data_select,adjust="none")
pvalues = cor_out$p
cor_matrix = cor_out$r
index_not_significant = which(pvalues>0.05)
cor_matrix[index_not_significant] = 0
cor_matrix_vec = cor_matrix[upper.tri(cor_matrix)]
mean(abs(cor_matrix_vec))
sparsity = length(which(cor_matrix_vec==0))/(p*(p-1)/2)
sparsity
occ_big = (length(which(cor_matrix_vec>0.25))+length(which(cor_matrix_vec<(-0.25))))/(p*(p-1)/2)

#partial correlation statistics
p = dim(plinks)[1]
cutoff = 0.5
index_exclude = which(plinks<cutoff)
C[index_exclude] = 0
C_vec = C[upper.tri(C)]
mean(abs(C))
sparsity= length(which(C_vec==0))/(p*(p-1)/2)
sparsity
occ_big = (length(which(C_vec>0.25))+length(which(C_vec<(-0.25))))/(p*(p-1)/2)

########################################################
##get edge probs and part. corr of a certain link#######
########################################################
var1 = "Age"
var2 = "Memory"
index1 = which(colnames(plinks)==var1)
index2 = which(colnames(plinks)==var2)
plinks[index1,index2]
C[index1,index2]

#########################
###create barcharts #####
#########################

#select columns related to variable of interest
title = "Memory"
index = which(colnames(plinks)==title)
edge_incl_vec = plinks[,index]
part_cor_vec = C[,index]

#fill data frame with edge incl. prob and partial correlations
darkblue = "#0000FFA9"
darkred = "#FF0000A9"
data_fr = data.frame(name=colnames(plinks),edge_incl=edge_incl_vec,part_cor=part_cor_vec)
data_fr$Color = ifelse(data_fr$part_cor <0,darkblue,darkred)

#remove row containing main variable 
index_main = which(data_fr$name==title)
data_fr = data_fr[-c(index_main),]

#remove rows with a edge_incl_prob less than 0.1
remove_cutoff = 0.1
index_remove = which(data_fr$edge_incl<remove_cutoff)
data_fr = data_fr[-c(index_remove),]

#if only volume, remove all non-volume variables
only_volume = 0
only_meta = 0
if (only_volume==1){
  volume_index = grep("V ",data_fr$name)
  data_fr = data_fr[volume_index,]
}
if (only_meta==1){
  meta_index = grep("G ",data_fr$name)
  data_fr = data_fr[meta_index,]
}

#make barchart for edge inclusion probabilities
plot_edge_incl = ggplot(data_fr, aes(x = reorder(name,+edge_incl), y = edge_incl)) + coord_flip()
plot_edge_incl = plot_edge_incl + geom_bar(stat="identity", color='blue',fill='blue') + xlab("") +ylab("Edge incl. prob")
plot_edge_incl

#make barchart for partial correlations
plot_part_cor = ggplot(data_fr, aes(x = reorder(name,+edge_incl), y = part_cor,fill=Color)) + coord_flip()
plot_part_cor = plot_part_cor + geom_bar(stat="identity") +xlab("") + ylab("Partial correlation") + ylim(-1,1) + theme(legend.position="none")
plot_part_cor

grid.arrange(plot_edge_incl,plot_part_cor,ncol=2)

################################
######circle plot over time ####
################################

highlight_type="glu_vol" #select either sex_cog, age_cog, bio_cog, dem_bio, glu_vol
for (diagnosis in c(1,2,3,4)){
  #load file
  filename = paste("Copula_output_diag_",diagnosis,".Rdata",sep="")
  load(file=filename)
  
  #create plinks, K and C and data (all with correct names and ordering of the rows and columns)
  out = create_P_K_C(output_alg=output_alg,data=data_select)
  plinks = out$plinks
  K = out$K
  C = out$C
  data = out$data
  
  #create circle plot
  highlight=TRUE
  par(mar = c(0,0,0,0))
  make_circle(plinks=plinks,C=C,highlight=highlight,highlight_type=highlight_type,cutoff=0.5)
}
