
#This file does three things
#1. It loads the raw data from the excel file
#2. it cleans the data, but keeps the left and right brain regions separate
#3. It runs the algorithm
#4. It creates the figure S2 in the supplementary material

####################################################
############# libraries ##########################
####################################################

library(readxl) #for loading excel file to R
library(circlize)
library(BDgraph)


##################################################################
#####1. load data and clean data#####################################
##################################################################

#load data from excel file to R
setwd("C:/Users/lvogels/OneDrive - UvA/PhD-Project-Lucas/1_Projects/5. JAD/2. Data/Github")
data_orig <- read_excel("2020_master_data_merged.xlsx")
data_1 = data_orig

####################################################
#############2. clean data ########################
####################################################

#make a vector of all the variables of interested
brain_regions=c("L_Hippocampus","R_Hippocampus","L_Thalamus","R_Thalamus","L_Caudate","R_Caudate",
                "L_Putamen","R_Putamen","Cingulate_Gyrus_posterior_L","Cingulate_Gyrus_posterior_R",
                "Precuneus_Cortex_L","Precuneus_Cortex_R")
metrics = c("GM","FDGpvc") #select either GM (for volume), AV45pvc (for amyloid), or FDGpvc (for metabolism)
region_variables = c()
for (metric in metrics){
  region_variables = c(region_variables,paste0(metric,"_",brain_regions))
}
other_variables = c("Age","ADNI_MEM","ADNI_EF","diagnosis","Sex","EDUC","APOE4","Amy_Stage")

all_variables = c(region_variables,other_variables)

#select only the data of the variables of interest
data_2 = data_1[all_variables]

#change character columns to numeric columns
sapply(data_2, class) #observe that ADNI_MEM, ADNI_EF, APOE4,Amy_stage are string values
data_2$ADNI_MEM = as.numeric(data_2$ADNI_MEM)
data_2$ADNI_EF = as.numeric(data_2$ADNI_EF)
data_2$APOE4 = as.numeric(data_2$APOE4)
data_2$Amy_Stage = as.numeric(data_2$Amy_Stage)
sapply(data_2, class) #observe that all values are numeric now

#check where the non-applicable values are
na_vec = c()
for (variable_name in all_variables){
  na_vec = c(na_vec,sum(is.na(data_2[[variable_name]])))  
}
na_vec = matrix(c(all_variables,na_vec),ncol=2,byrow=FALSE)
na_vec

#remove all rows containing NAs
sum(is.na(data_2))
data_3 = na.omit(data_2) #remove rows containing non-applicable values
rows_removed = dim(data_2)[1] - dim(data_3)[1]
rows_removed
dim(data_3) 

#make the APOE4 variable binary
APOE4_vec = data_3$APOE4
index_0 = which(APOE4_vec==0)
length(index_0)
index_1 = which(APOE4_vec==1)
length(index_1)
index_2 = which(APOE4_vec==2)
length(index_2)
APOE4_vec[index_2] = 1 
data_3$APOE4 = APOE4_vec


#check number of individuals in each subgroup
diagnosis_encoding = c("cognitively normal","early mild cognitive impairment","late mild cognitive impairment","dementia due to Alzheimer's disease","Subjective cognitive decline")
for (i in 1:5){
  diagnosis_vec = data_3$diagnosis == i
  total_diag = sum(diagnosis_vec)
  text= paste0(diagnosis_encoding[i]," has ",total_diag, " observations")
  print(text)
}

#merge diagnosis 1 and 5
data_4 = data_3
diag_5_index = which(data_4$diagnosis==5)
data_4$diagnosis[diag_5_index] = 1

#merge brain regions
#brain regions are merged when they are positioned in the center of the brain and practically form one region
data_5= data_4

#rename variables
old_names = colnames(data_5)
new_names = c("V_Hipp_L","V_Hipp_R",
              "V_Thal_L","V_Thal_R",
              "V_Caud_L","V_Caud_R",
              "V_Put_L","V_Put_R",
              "V_PCC_L","V_PCC_R",
              "V_Prec_L","V_Prec_R",
              "G_Hipp_L","G_Hipp_R",
              "G_Thal_L","G_Thal_R",
              "G_Caud_L","G_Caud_R",
              "G_Put_L","G_Put_R",
              "G_PCC_L","G_PCC_R",
              "G_Prec_L","G_Prec_R",
              "Age","Memory","Executive","Diagnosis","Sex","Educ","APOE4","Amy-stage")
matrix(c(old_names,new_names),byrow=FALSE,ncol=2)
colnames(data_5) = new_names

#######################################################
###########3. run algorithms ############################
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
  last_K_iter_matrix = matrix(0,p*(p-1)/2+p,len_iter)
  weights_iter = rep(0,len_iter)
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
    
    #save last K in matrix
    last_K = sample_bd$last_K
    last_K_vec = last_K[upper.tri(last_K,diag=TRUE)]
    last_K_iter_matrix[,j] = last_K_vec
    
    #save weights
    weights_iter[j] = weights_new
  }
  
  plinks = plinks_new
  return(list(g.prior=g.prior,edge_vec=edge_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new,plinks_iter_matrix=plinks_iter_matrix,K_iter_matrix=K_iter_matrix,last_K_iter_matrix=last_K_iter_matrix,weights_iter=weights_iter))
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

#create all networks
diagnosis_vec = c(1,2,3,4) #select subset of diagnosis that you want included, either 1,2,3,4 or a combination
title = paste0("Diagnosis ",paste(diagnosis_vec,collapse=" "))
print(title)
data_select = data_select_func(data=data_5,diagnosis_vec=diagnosis_vec)

#normalize continuous data
not.cont=c(rep(0,24),0,0,0,1,1,1,1)
matrix(c(colnames(data_select),not.cont),ncol=2) #check if non-continuous vector is correct (age assumed contninuos)
method = "truncation" #select either shrinkage, truncation or scaling
data_norm = Gaussianize_func(data=data_select,not.cont=not.cont,method=method)

#select the burnin and MCMC iterations, and run the algorithm
burnin_iter_vec_thin = c(10,50,100,500,1000,2500,5000,10000,20000)
iter_vec_thin = c(100,200,500,1000,2500,5000,7500,10000,15000,20000,30000,40000,50000,60000,70000,80000,90000,100000)
output_alg = run_algorithm_func(data=data_norm,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,not.cont=not.cont)

#######################################################################
############4. create figures##########################################
#######################################################################

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
  new_names = c("V Hipp L","V Hipp R",
                "V Thal L","V Thal R",
                "V Caud L","V Caud R",
                "V Put L","V Put R",
                "V PCC L","V PCC R",
                "V Prec L","V Prec R",
                "G Hipp L","G Hipp R",
                "G Thal L","G Thal R",
                "G Caud L","G Caud R",
                "G Put L","G Put R",
                "G PCC L","G PCC R",
                "G Prec L","G Prec R",
                "Age","Memory","Executive","Sex","Educ","APOE4","Amy-stage")      
  matrix(c(old_names,new_names),byrow=FALSE,ncol=2)
  colnames(plinks) = new_names
  rownames(plinks) = new_names
  colnames(K) = new_names
  rownames(K) = new_names
  colnames(C) = new_names
  rownames(C) = new_names
  
  #reorder columns and rows
  order_vec=c(1:24,25,28,29,30,31,26,27)
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
               "G Thal L" = "#999999",
               "G Thal R" = "#999999",
               "V Thal L" = "#666666",
               "V Thal R" = "#666666",
               "G PCC L" = "#999999",
               "G PCC R" = "#999999",
               "V PCC L" = "#666666",
               "V PCC R" = "#666666",
               "G Prec L" = "#999999",
               "G Prec R" = "#999999",
               "V Prec L"= "#666666",
               "V Prec R"= "#666666",
               "G Hipp L" = "#999999",
               "G Hipp R" = "#999999",
               "V Hipp L"= "#666666",
               "V Hipp R"= "#666666",
               "G Caud L" = "#999999",
               "G Caud R" = "#999999",
               "V Caud L"= "#666666",
               "V Caud R"= "#666666",
               "G Put L" = "#999999",
               "G Put R" = "#999999",
               "V Put L"   = "#666666",
               "V Put R"   = "#666666")
  
  #set gaps between borders
  circos.par(gap.after = 
               c(Age = 1,
                 Memory = 1,
                 Executive = 5,
                 Sex = 1,
                 Educ = 5,
                 APOE4 = 1,
                 "Amy-stage" = 5,
                 "G Thal L" = 0.25,
                 "G Thal R" = 1,
                 "V Thal L" = 0.25,
                 "V Thal R" = 1,
                 "G PCC L" = 0.25,
                 "G PCC R" = 1,
                 "V PCC L" = 0.25,
                 "V PCC R" = 1,
                 "G Prec L" = 0.25,
                 "G Prec R" = 5,
                 "V Prec L"= 0.25,
                 "V Prec R"= 5,
                 "G Hipp L" = 0.25,
                 "G Hipp R" = 1,
                 "V Hipp L"= 0.25,
                 "V Hipp R"= 1,
                 "G Caud L" = 0.25,
                 "G Caud R" = 1,
                 "V Caud L"= 0.25,
                 "V Caud R"= 1,
                 "G Put L" = 0.25,
                 "G Put R" = 1,
                 "V Put L"   = 0.25,
                 "V Put R"   = 1)
  )
  
  #intialize starting degree
  circos.par(start.degree = 40)
  
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

#perform convergence analysis
title = paste0("Diagnosis ",paste(diagnosis_vec,collapse=" "))
convergence_metric = "P" #select either P or K
conv_stat_func(out=output_alg,converge_metric = converge_metric,title=title)

#create plinks, K and C and data (all with correct names and ordering of the rows and columns)
out = create_P_K_C(output_alg=output_alg)
plinks = out$plinks
K = out$K
C = out$C

#create circle plot
highlight=FALSE 
highlight_type= "NA" #select either sex_cog, age_cog, bio_cog, dem_bio, glu_vol
par(mar = c(0,0,0,0))
make_circle(plinks=plinks,C=C,highlight=highlight,highlight_type=highlight_type,cutoff=0.5)

#some statistics for the article
p_links_vec = plinks[upper.tri(plinks)]
qp = length(p_links_vec)
length(which(p_links_vec>0.5))/qp #% of edges that have an edge incl prob higher than 0.5
C_vec = C[upper.tri(C)]
index_exclude = which(p_links_vec<0.5)
C_vec[index_exclude] = 0
mean(abs(C_vec))
