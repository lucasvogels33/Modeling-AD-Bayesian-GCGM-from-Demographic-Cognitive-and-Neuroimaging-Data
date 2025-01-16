##################################################################################################################################################
######this file does three things: i) load the data from an excel file, ii) clean the data, and iii) save the cleaned data as .Rdata file########
##################################################################################################################################################

####################################################
############# i) load data#########################
####################################################

library(readxl) #for loading excel file to R
setwd("C:/Users/lvogels/OneDrive - UvA/PhD-Project-Lucas/1_Projects/5. JAD/2. Data/Github")
data_orig <- read_excel("2020_master_data_merged.xlsx")
data_1 = data_orig

####################################################
#############ii) clean data ########################
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

#merge diagnosis 1 (Cognitivly normal) and 5 (subjective decline)
data_4 = data_3
diag_5_index = which(data_4$diagnosis==5)
data_4$diagnosis[diag_5_index] = 1

#merge brain regions
merge_func = function(left,right,merge,metrics,df){
  for (metric in metrics){
    left_name = paste0(metric,"_",left)
    right_name = paste0(metric,"_",right)
    merge_name = paste0(metric,"_",merge)
    df[merge_name] = (df[left_name] + df[right_name])/2
    df[left_name] = NULL
    df[right_name] = NULL
  }
  return(df)
}
data_5 = data_4
metrics = c("FDGpvc","GM") # select either GM (for volume) or FDGpvc (for metabolism)

left = "L_Thalamus" 
right = "R_Thalamus"
merge = "Thalamus_merge"
data_5 = merge_func(left=left,right=right,merge=merge,metrics=metrics,df=data_5)

left = "Cingulate_Gyrus_posterior_L" 
right = "Cingulate_Gyrus_posterior_R"
merge = "Cingulate_Gyrus_posterior_merge"
data_5 = merge_func(left=left,right=right,merge=merge,metrics=metrics,df=data_5)

left = "Precuneus_Cortex_L" 
right = "Precuneus_Cortex_R"
merge = "Precuneus_Cortex_merge"
data_5 = merge_func(left=left,right=right,merge=merge,metrics=metrics,df=data_5)

left = "L_Hippocampus" 
right = "R_Hippocampus"
merge = "Hippocampus_merge"
data_5 = merge_func(left=left,right=right,merge=merge,metrics=metrics,df=data_5)

left = "L_Caudate"
right = "R_Caudate"
merge = "Caudate_merge"
data_5 = merge_func(left=left,right=right,merge=merge,metrics=metrics,df=data_5)

left = "L_Putamen"
right = "R_Putamen"
merge = "Putamen_merge"
data_5 = merge_func(left=left,right=right,merge=merge,metrics=metrics,df=data_5)

#rename variables
old_names = colnames(data_5)
new_names = c("Age","Memory","Executive","Diagnosis","Sex","Educ","APOE4","Amy-stage",
              "G Thal","V Thal",
              "G PCC","V PCC",
              "G Prec","V Prec",
              "G Hipp","V Hipp",
              "G Caud","V Caud",
              "G Put","V Put")              
matrix(c(old_names,new_names),byrow=FALSE,ncol=2)
colnames(data_5) = new_names


#make table 1 in the paper
diagnosis_groups = c(1,2,3,4)
variable_name = "ADNI-EF" #select one of the following variables: ADNI-EF, ADNI-MEM, Amyloid stage, Sex, APOE4, EDUC, Age
if (variable_name == "Amyloid stage"){var_vec = data_5$Amy-Stage}
if (variable_name == "Sex"){var_vec = data_5$Sex}
if (variable_name == "APOE4"){var_vec = data_5$APOE4}
if (variable_name == "EDUC"){var_vec = data_5$Educ}
if (variable_name == "Age"){var_vec = data_5$Age}
if (variable_name == "ADNI-MEM"){var_vec = data_5$Memory}
if (variable_name == "ADNI-EF"){var_vec = data_5$Executive}

for (diag in diagnosis_groups){
  index_diag = (data_5$Diagnosis==diag)
  average = mean(var_vec[index_diag])
  sd = sd(var_vec[index_diag])
  total_female = sum(var_vec[index_diag])
  female_perc = total_female/sum(index_diag)
  text = paste0("diagnosis group ",diag," has an average of ", variable_name," of ", average, "and a standard deviation of ", sd)
  if (variable_name=="Sex"){
    text = paste0("diagnosis group ",diag," has a total of ", total_female," women. That is ",female_perc)  
  }
  print(text)
}  


#######################
####save data #########
#######################
filename = "Cleaned_data.Rdata"
save(data_5,file=filename)
