####   Instructions   ####
# Interaction calculation requires robust linear regression, so we use rlm from the MASS package
# Species abundance and environmental data should both be data frames, where rows represent samples and columns represent species names or environmental parameter names
# The data should have no NA values and no columns with too many repeated values (e.g., zeros), thus data cleaning is required
library(utils)   # Provides txtProgressBar function
library(MASS)    # this package includes the rlm command 
library(reshape)
library(ggplot2)
library(robust)
library(openxlsx)

###############################  
#####         step 1      #####  
###############################  

load("Abundscale.RData")   # species abundance  
load("Parascale.rdata")    # environmental parameters
source("InteractionAnalysis_git.R")

###############################  
#####        step 2       #####  
###############################  
# Remove singular values
RelaAbund<-testsingula(Abundscale)
Parascale<-testsingula(Parascale)

###############################  
#####   step 3 Debug mode  ######  ## Requires at least around 200 inputs, time-consuming, temporarily disabled
###############################  
debug_mode <- FALSE    # Set to TRUE to enable debug mode (use subset of data), FALSE to run full dataset
debug_rows <- 150      # Number of rows to use in debug mode

if (debug_mode) {
  message(paste("Debug mode enabled, using only the first", debug_rows, "rows of data"))
  # Subset data for debugging
  if (nrow(Abundscale) > debug_rows) {
    Abundscale <- Abundscale[1:debug_rows, ]
  }
  if (nrow(Parascale) > debug_rows) {
    Parascale <- Parascale[1:debug_rows, ]
  }
}

###############################  
#####        step 4       #####  
###############################  
# Define parameters needed in the function
speciesprecision="high"
envprecision="high"
threshold01<-0.5   # used to determine the pattern in each environmental parameter
threshold<-0.3     # error level used in the robust test
type<-3            # use the workflow method to find the global interaction values
count<-15          # used to calculate the peak value of the distribution 
factor02<-2        # used to compare the difference between the positive and negative part

###############################  
#####        step 5       #####  
###############################  
# Use the modified function (add returnPerEnv=TRUE parameter)
result <- interMatrix02(
  RelaAbund, Parascale, 
  type = 3, count = 10, factor02 = 2,
  speciesprecision = "high", envprecision = "high", 
  threshold01 = 0.5,
  returnPerEnv = TRUE  # New parameter
)

# Extract the global matrix and results for each environmental factor
InterMatrix <- result$globalMatrix
perEnvMatrices <- result$perEnvMatrices
# Calculate Pearson correlation matrix
corrMatrix <- cor(Abundscale)
# Save the global matrix
save(InterMatrix, file = "GlobalMatrix.RData")

###############################  
#####        step 6       #####  
###############################  
# Generate filenames with timestamp (accurate to minute)
timestamp <- format(Sys.time(), "%Y%m%d_%H")

# Save global interaction matrix and Pearson correlation matrix
if (debug_mode) {
  global_result_file <- paste0("Global_Debug_", debug_rows, "rows_", timestamp, ".RData")
} else {
  global_result_file <- paste0("Global_Full_", timestamp, ".RData")
}

# Save global matrix and Pearson matrix
save(InterMatrix, corrMatrix, file = global_result_file)
message(paste("Global interaction matrix and Pearson correlation matrix saved to:", global_result_file))

# Save interaction matrix for each environmental factor
env_names <- names(perEnvMatrices)
for (env_name in env_names) {
  env_matrix <- perEnvMatrices[[env_name]]
  
  if (debug_mode) {
    env_result_file <- paste0("Env_", env_name, "_Debug_", debug_rows, "rows_", timestamp, ".RData")
  } else {
    env_result_file <- paste0("Env_", env_name, "_Full_", timestamp, ".RData")
  }
  
  save(env_matrix, file = env_result_file)
  message(paste("Interaction matrix for environmental factor", env_name, "saved to:", env_result_file))
}

###############################  
#####        Export Excel     #####  
###############################  

# Export global interaction matrix
global_inter_xlsx <- if (debug_mode) {
  paste0("Global_InteractionMatrix_Debug_", debug_rows, "rows_", timestamp, ".xlsx")
} else {
  paste0("Global_InteractionMatrix_Full_", timestamp, ".xlsx")
}
write.xlsx(InterMatrix, file = global_inter_xlsx)
message(paste("Global interaction matrix exported to:", global_inter_xlsx))

# Export global Pearson correlation matrix
global_corr_xlsx <- if (debug_mode) {
  paste0("Global_CorrelationMatrix_Debug_", debug_rows, "rows_", timestamp, ".xlsx")
} else {
  paste0("Global_CorrelationMatrix_Full_", timestamp, ".xlsx")
}
write.xlsx(corrMatrix, file = global_corr_xlsx)
message(paste("Global Pearson correlation matrix exported to:", global_corr_xlsx))

# Export matrix for each environmental factor
for (env_name in env_names) {
  env_matrix <- perEnvMatrices[[env_name]]
  env_xlsx <- if (debug_mode) {
    paste0("Env_", env_name, "_Debug_", debug_rows, "rows_", timestamp, ".xlsx")
  } else {
    paste0("Env_", env_name, "_Full_", timestamp, ".xlsx")
  }
  write.xlsx(env_matrix, file = env_xlsx)
  message(paste("Interaction matrix for environmental factor", env_name, "exported to:", env_xlsx))
}


###############################  
####         step 7      ######  
###############################   
#### Select six examples for the robust test: strong positive, median positive, low positive, ####
#### strong negative, median negative, low negative. #####  
#### The results (row and column) are the indices of the two species.

robust_test_rela<-six_test_robust(InterMatrix)  



###############################  
####         step 8      ######  
############################### 

#### Robust test on the six examples ###

#### k indicates the index of the six examples
for(k in 1:6){
  species_i<- robust_test_rela[k,1] ##### species i
  species_j<- robust_test_rela[k,2] ##### species j
  
  # ##### Remove one species #####
  # ###### Remove one environmental parameter ######
  #   permut_type<-"para"
  # 
  # permut_para<-interMatrix_ij_permut02 (Abundscale,Parascale,species_i,species_j,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
  # nam1<-paste("permut_para", rownames(robust_test_rela)[k],sep="_")  
  # assign(nam1, permut_para);
  # 
  
  ###### Remove samples #####
  permut_type<-"samples"
  #times<-500
  times<-5 ##### High frequency robust test required
  size<-15
  permut_sample<- interMatrix_ij_permut(RelaAbund,Parascale,species_i,species_j,times,size,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
  
  nam1<-paste("permut_sample", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, permut_sample);
  
  
  ###### Error robust on parameters ####
  #times<-500
  times<-5 ##### Only for testing
  robust_para<- interMatrix_ij_env_robust(RelaAbund,Parascale,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01)
  nam1<-paste("robust_para", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, robust_para);
  
  
  ###### Error robust on species ####
  
  #times<-500
  times<-5 ##### Only for testing
  robust_spe<- interMatrix_ij_spe_robust(RelaAbund,Parascale,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01)
  nam1<-paste("robust_spe", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, robust_spe);
}

###############################  
#####       Addition      #####  
###############################
# Identify dominant parameters: some environmental parameters play a leading role in determining interaction levels
# We calculate (1) frequency of parameters determining interaction values;
# (2) frequency of parameters with the smallest negative values;
# (3) frequency of parameters with the largest positive values.

out_Sample<-main_para(RelaAbund,Parascale,0,10,2,speciesprecision,envprecision,threshold01)
main_para_barplot(out_Sample)