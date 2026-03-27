## 在112仅做因子修改，每次复制最新的代码即可
####  使用说明   ####
# 相互作用的计算需要鲁棒线性回归，因此，我们在包MASS中使用rlm
# 物种丰度和环境数据都应该是数据框架，其中行代表样本，列代表物种名称或环境参数名称
# 数据中不应该有NA值，也不应该有包含太多重复值的列（例如0），因此，我们需要过滤掉数据

library(utils)   # 提供txtProgressBar函数
library(MASS)    # this package include the rlm command 
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
# 删除奇异值
RelaAbund<-testsingula(Abundscale)
Parascale<-testsingula(Parascale)

###############################  
#####   step 3 调试模式  ######  ## 最少需要200左右的输入，时间并不短，暂且关闭
###############################  
debug_mode <- FALSE    # 设为TRUE启用调试模式（使用部分数据），FALSE 则运行完整数据
debug_rows <- 150      # 调试时使用的数据行数

if (debug_mode) {
  message(paste("调试模式已启用，仅使用前", debug_rows, "行数据"))
  # 截取部分数据用于调试
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
# 定义函数中需要用到的参数
speciesprecision="high"
envprecision="high"
threshold01<-0.5   # used to determin the patern in each environmental parameter
threshold<-0.3     # error level used in the robust test
type<-3            # use the work flow method to find the global interaction values
count<-15          # used to calculate the peakvalue of the distribution 
factor02<-2      # used to compare the difference between the positive and negative part

###############################  
#####        step 5       #####  
###############################  
# 使用修改后的函数（增加 returnPerEnv=TRUE 参数）
result <- interMatrix02(
  RelaAbund, Parascale, 
  type = 3, count = 10, factor02 = 2,
  speciesprecision = "high", envprecision = "high", 
  threshold01 = 0.5,
  returnPerEnv = TRUE  # 新增参数
)

# 提取全局矩阵和每个环境因子的结果
InterMatrix <- result$globalMatrix
perEnvMatrices <- result$perEnvMatrices
# 计算皮尔逊相关矩阵
corrMatrix <- cor(Abundscale)
# 保存全局矩阵
save(InterMatrix, file = "GlobalMatrix.RData")

###############################  
#####        step 6       #####  
###############################  
# 生成包含时间戳的文件名（精确到分钟）
timestamp <- format(Sys.time(), "%Y%m%d_%H")

# 保存全局相互作用矩阵和皮尔逊相关矩阵
if (debug_mode) {
  global_result_file <- paste0("Global_Debug_", debug_rows, "rows_", timestamp, ".RData")
} else {
  global_result_file <- paste0("Global_Full_", timestamp, ".RData")
}

# 保存全局矩阵和皮尔逊矩阵
save(InterMatrix, corrMatrix, file = global_result_file)
message(paste("全局相互作用矩阵和皮尔逊相关矩阵已保存到:", global_result_file))

# 保存每个环境因子的相互作用矩阵
env_names <- names(perEnvMatrices)
for (env_name in env_names) {
  env_matrix <- perEnvMatrices[[env_name]]
  
  if (debug_mode) {
    env_result_file <- paste0("Env_", env_name, "_Debug_", debug_rows, "rows_", timestamp, ".RData")
  } else {
    env_result_file <- paste0("Env_", env_name, "_Full_", timestamp, ".RData")
  }
  
  save(env_matrix, file = env_result_file)
  message(paste("环境因子", env_name, "的相互作用矩阵已保存到:", env_result_file))
}

###############################  
#####       输出Excel     #####  
###############################  

# 导出全局相互作用矩阵
global_inter_xlsx <- if (debug_mode) {
  paste0("Global_InteractionMatrix_Debug_", debug_rows, "rows_", timestamp, ".xlsx")
} else {
  paste0("Global_InteractionMatrix_Full_", timestamp, ".xlsx")
}
write.xlsx(InterMatrix, file = global_inter_xlsx)
message(paste("全局相互作用矩阵已导出到:", global_inter_xlsx))

# 导出全局皮尔逊相关矩阵
global_corr_xlsx <- if (debug_mode) {
  paste0("Global_CorrelationMatrix_Debug_", debug_rows, "rows_", timestamp, ".xlsx")
} else {
  paste0("Global_CorrelationMatrix_Full_", timestamp, ".xlsx")
}
write.xlsx(corrMatrix, file = global_corr_xlsx)
message(paste("全局皮尔逊相关矩阵已导出到:", global_corr_xlsx))

# 导出每个环境因子的矩阵
for (env_name in env_names) {
  env_matrix <- perEnvMatrices[[env_name]]
  env_xlsx <- if (debug_mode) {
    paste0("Env_", env_name, "_Debug_", debug_rows, "rows_", timestamp, ".xlsx")
  } else {
    paste0("Env_", env_name, "_Full_", timestamp, ".xlsx")
  }
  write.xlsx(env_matrix, file = env_xlsx)
  message(paste("环境因子", env_name, "的相互作用矩阵已导出到:", env_xlsx))
}


###############################  
####         step 7      ######  
###############################   
#### chose the six examples for the robust test::  strong positive, median positive, low positive,####
####  strong negative, median negative, low negative.  #####  
#### the results (row and column) are the index of the two species, 

robust_test_rela<-six_test_robust(InterMatrix)  



###############################  
####         step 8      ######  
############################### 

#### robust test on the six examples ###

#### k indicates the index of the six examples
for(k in 1:6){
  species_i<- robust_test_rela[k,1] ##### species i
  species_j<- robust_test_rela[k,2] ##### species j
  
  # ##### 移除一个物种 #####
  # ###### remove one environmental parameters ######
  #   permut_type<-"para"
  # 
  # permut_para<-interMatrix_ij_permut02 (Abundscale,Parascale,species_i,species_j,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
  # nam1<-paste("permut_para", rownames(robust_test_rela)[k],sep="_")  
  # assign(nam1, permut_para);
  # 
  
  ###### remove  samples #####
  permut_type<-"samples"
  #times<-500
  times<-5 ##### 需要进行高频次的robust检验
  size<-15
  permut_sample<- interMatrix_ij_permut(RelaAbund,Parascale,species_i,species_j,times,size,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
  
  nam1<-paste("permut_sample", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, permut_sample);
  
  
  ###### error robust on para ####
  #times<-500
  times<-5 ##### 仅为测试
  robust_para<- interMatrix_ij_env_robust(RelaAbund,Parascale,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01)
  nam1<-paste("robust_para", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, robust_para);
  
  
  ###### error robust on species ####
  
  #times<-500
  times<-5 ##### 仅为测试
  robust_spe<- interMatrix_ij_spe_robust(RelaAbund,Parascale,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01)
  nam1<-paste("robust_spe", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, robust_spe);
}

# 下面这个暂时用不到不运行，运行时间很长
###############################  
#####       Addition      #####  
###############################
# 确定主导参数：某些环境参数在决定交互水平方面起主导作用 
# 我们计算(1)决定交互值的参数频率； 
# (2)负值最小的参数频率； 
# (3)正值最大的参数出现的频率。 

#out_Sample<-main_para(RelaAbund,Parascale,0,10,2,speciesprecision,envprecision,threshold01)
#main_para_barplot(out_Sample)
