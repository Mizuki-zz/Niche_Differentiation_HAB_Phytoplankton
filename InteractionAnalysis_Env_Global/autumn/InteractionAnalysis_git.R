######                       代码概览                         ######
##### 这个R文件包含交互分析中需要的所有函数
#####  part A --  过滤数据: FilterNA, FilterZero, Filterrepeat
#####  part B --  计算导数和原始相互作用值: rate_change, rate_change02, rate_change_low, interinf, interinf02, interinf03, interinf_low
#####  part C --  计算 interaction matirx:  interMatrix, interMatrix02, interMatrix_norm, interMatrix_norm_summ, interMatrix_low, interMatrix_ij, interMatrix_ij_low
#####  part D --  :  summ_inter_ij, summ_inter_ij_new, summ_interaction, anascatter
#####  part E --  main environmental parameters:  main_para, main_para_barplot
#####  part F --  robust test:  interMatrix_ij_permut, interMatrix_ij_env_spe_robust, interMatrix_ij_env_robust, interMatrix_ij_spe_robust, interMatrix_ij_permut02  
#####  part G --  test singularity: testsingula
#####  part H --  extract six examples from the interaction matrix: six_test_robust 
#####  part I --  make plot of change of interaciton value across the samples:  reshapedata, interinfPlot




####################################################################
######          part A：过滤丢失信息最少的数据                ######
####################################################################
FilterNA <- function (data,type){
  #### 目的：用于根据指定的方式（删除缺失值的列或行）清洗数据中的缺失值（NA）
  #### data is data frame of species abundance or environmental parameter 
  #### type = 0，删除包含缺失值的列或行，优先删除包含缺失值最多的列或行
  #### type = 1，删除包含缺失值的行
  #### type = other value, 删除包含缺失值的列
  if (type==0){
    while (anyNA(data)){
      NAcol_env<- apply(data, 2, function(x) length(which(is.na(x))))
      NArow_env<- apply(data, 1, function(x) length(which(is.na(x))))
      max_NAcol_env <- max(NAcol_env,na.rm=TRUE)
      max_NArow_env <-max(NArow_env,na.rm=TRUE)
      if (max_NAcol_env>=max_NArow_env){
        data <-data[,-which(NAcol_env==max_NAcol_env)]}  else{
          data <-data[-which(NArow_env==max_NArow_env),]
        }
    }
  } else if (type==1){
    NArow<- apply(data, 1, function(x) length(which(is.na(x))))
    data <- data[which(NArow==0),] 
  }else{
    NAcol<- apply(data, 2, function(x) length(which(is.na(x))))
    data <- data[,which(NAcol==0)]  
  }
  return(data)
}


FilterZero <- function (data,type){
  #### 清洗数据中的零值
  #### data is data frame of species abundance or environmental parameter 
  #### type = 0，删除包含零值的列或行，优先删除包含零值最多的列或行
  #### type = 1，删除包含零值的行
  #### type = other value, 删除包含零值的列
  if (type==0){
    while (any(data==0.0)){
      Zerocol_env<- apply(data, 2, function(x) length(which(x==0.0)))
      Zerorow_env<- apply(data, 1, function(x) length(which(x==0.0)))
      max_Zerocol_env <- max(Zerocol_env,na.rm=TRUE)
      max_Zerorow_env <-max(Zerorow_env,na.rm=TRUE)
      if (max_Zerocol_env>=max_Zerorow_env){
        data <-data[,-which(Zerocol_env==max_Zerocol_env)]}  else{
          data <-data[-which(Zerorow_env==max_Zerorow_env),]
        }
    }
  } else if (type==1){
    Zerorow<- apply(data, 1, function(x) length(which(x==0.0)))
    data <- data[which(Zerorow==0.0),] 
  }else{
    Zerocol<- apply(data, 2, function(x) length(which(x==0.0)))
    data <- data[,which(Zerocol==0.0)]  
  }
  return(data)
}


Filterrepeat <- function(data,threshold){
  #### remove the columns which contain too many repeat values
  N <- floor(threshold*nrow(data))  #### the threshold for the munber of repeat values
  reptimes <- vector()
  for (j in 1: ncol(data)){
    reptimes <- c(reptimes, max(table(data[,j])))
  }
  data <- data[,which(reptimes< N)]
  return(data)
}


#################################################################
######           part B：计算一阶导数和相互作用            ######
#################################################################  
#### 这部分包含了计算物种丰度变化率和原始四维相互作用值的所有函数，
#### 变化率的计算需要物种丰度和环境参数数据。交互还需要变化率的计算结果作为输入数据。
rate_change <- function(Y, X, i){
  #### 该函数将计算物种i丰度对所有样本中所有环境参数的一阶导数，
  #### 根据数据结构的不同，该函数将自动选择不同的分割方法。返回的结果是二维的
  #### dataframe: row is the samples, column is the environmental parameter
  ####   Y is the abundance, data frame 
  ####   X is the environmental parameters, data frame
  ####   i is the index of species i
  nsample <- nrow(Y)
  nSpe <-ncol(Y)
  nEnv <- ncol(X)
  derEnv <- matrix(, nrow = nsample-1, ncol = nEnv)
  #  主循环：对每一个样本 n（作为基准点），与其余样本比较
  for (n in 1:(nsample -1)){
    # change of  abundance
    # deltaSpe 是一个矩阵，表示每个物种在各个样本与样本 n 之间的差异
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,]))))
    # deltaSpe_i 是我们关心的第 i 个物种的变化量
    deltaSpe_i <- deltaSpe[,i]
    # change of environmental parameters
    # deltaEnv表示样本间环境参数的差异
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,]),1,function(x) x-as.matrix(X[n,]))))
    # 情况1：样本数大于物种数 + 环境变量数（推荐使用线性回归）
    if (nsample>=(nSpe+nEnv)){
      # 回归模型：将目标物种变化量（Spe_i）回归于所有环境变量 + 其他物种变化（多元回归）
      regdata <- data.frame(cbind( deltaEnv,deltaSpe[,-i]))
      fmal<- as.formula (paste("deltaSpe_i ~ -1 + ",paste(colnames(regdata),collapse="+")))
      fit <- lm (fmal,  regdata)    
      derEnv[n,]<- fit$coefficients[1:nEnv]
    } else if (nsample>=nEnv){
      # 情况2：样本数 >= 环境变量数，但不足以同时回归环境变量和其他物种变化
      # 只使用环境变量来拟合（robust regression：鲁棒回归）
      regdata <- data.frame( deltaEnv)
      fmal<- as.formula (paste("deltaSpe_i ~ -1 + ",paste(colnames(regdata),collapse="+")))
      fit <- rlm (fmal,  regdata)  
      derEnv[n,]<- fit$coefficients[1:nEnv]
    }else{
      # 情况3：样本太少，不能做有效回归，采用近似法（二点差分）
      temp<-apply(deltaEnv,2,function(x) deltaSpe_i/x)
      temp[which(!is.finite(temp))] <- NA
      derEnv[n,] <- apply(temp,2,function(x) mean(x,na.rm=TRUE))
    } 
  }
  derEnv <- data.frame(derEnv)
  rownames(derEnv) <-rownames(Y)[1:(nsample-1)]
  colnames(derEnv) <-colnames(X)
  return(derEnv)   #每一行表示一次样本差异的估计结果，每列是一个环境变量的“响应系数”或“导数”
}


rate_change02 <- function(Y, X, i){ 
  #### 该函数将计算物种i丰度对所有样本n中所有环境参数的一阶导数，
  #### 根据数据结构的不同，该函数只考虑环境参数的影响。返回的结果是二维的
  #### 该函数通过**鲁棒回归（robust linear regression）**来估计物种丰度对环境参数的响应程度，
  #### 其返回的结果是每个样本对每个环境参数的导数。
  #### data frame: row is the samples, column is the environmental parameter
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   i is the species i
  nsample <- nrow(Y)
  nSpe <-ncol(Y)
  nEnv <- ncol(X)
  # derEnv：用于保存每个样本的导数估计值，矩阵大小为 (nsample-1, nEnv)，即每行对应一个样本，列对应一个环境参数。
  derEnv <- matrix(, nrow = nsample-1, ncol = nEnv)
  # 计算每个样本 n 和其他样本之间物种丰度（Spe）和环境参数（Env）的变化量
  for (n in 1:(nsample -1)){
    # change of  abundance
    # deltaSpe_i 是第 i 个物种的丰度变化
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,]))))
    deltaSpe_i <- deltaSpe[,i]
    # change of environmental parameters
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,]),1,function(x) x-as.matrix(X[n,]))))
    
    
    regdata <- data.frame( deltaEnv)
    
    # use robust linear regression
    # lm一点也不健壮，rlm需要测试数据的奇异性
    # linear regression on   deltaSpe ~ deltaEnv
    fmal<- as.formula (paste("deltaSpe_i ~ -1 + ",paste(colnames(regdata),collapse="+")))
    fit <- rlm (fmal, regdata,  maxit = 500,acc=1e-3) 
    # 200 is  the running steps to have a convergence estimation, it can be changed based on the user?s test and experience
    # 该行代码提取鲁棒回归结果中的环境参数系数（即导数），并将其存储在 derEnv 中
    derEnv[n,]<- fit$coefficients[1:nEnv]
  }
  derEnv <- data.frame(derEnv)
  rownames(derEnv) <-rownames(Y)[1:(nsample-1)]
  colnames(derEnv) <-colnames(X)
  return(derEnv) 
}


rate_change_low <- function(Y, X, i){
  ####  该函数将计算物种i丰度对所有样本中所有环境参数的一阶导数，
  #### 根据数据结构的不同，该函数采用低精度方法。返回的结果是二维的
  #### data frame: row is the samples, column is the environmental parameter
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   i is the species i
  ####   this function will calculate the first derivative of species i abundace with respect to all the environmental parameters in the sample n
  nsample <- nrow(Y)
  nSpe <-ncol(Y)
  nEnv <- ncol(X)
  derEnv <- matrix(, nrow = nsample-1, ncol = nEnv)
  for (n in 1:(nsample -1)){
    # change of  abundance
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,]))))
    deltaSpe_i <- deltaSpe[,i]
    # change of environmental parameters
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,]),1,function(x) x-as.matrix(X[n,]))))
    # low precision, using two point to calculate the derivative, then choose the median 
    temp<-apply(deltaEnv,2,function(x) deltaSpe_i/x)
    temp[which(!is.finite(temp))] <- NA
    derEnv[n,] <- apply(temp,2,function(x) mean(x,na.rm=TRUE))
  }
  derEnv <- data.frame(derEnv)
  rownames(derEnv) <-rownames(Y)[1:(nsample-1)]
  colnames(derEnv) <-colnames(X)
  return(derEnv) 
}


interinf <- function(Y, X, DY, i){
  #### 该函数将计算所有样品中物种丰度的相互作用水平,
  #### 根据数据结构的不同，这个函数会自动选择不同的具体方法（鲁棒回归、普通回归、低精度方法）,
  #### 返回的结果是三维数组
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   DYY 对环境参数 X 的导数矩阵，表示每个样本中每个物种丰度相对于环境参数的变化率
  ####   i is the species i
  #####   this function will calculate the interaction influence  on species i 
  nsample <- nrow(DY)
  nEnv <- ncol(DY)
  nSpe <- ncol(Y)
  # 当 nsample >= (nSpe + nEnv),该模型考虑了物种丰度和环境参数的影响，采用泰勒展开
  if(nsample >=(nSpe+nEnv) ){
    # array: array[1] is the environmental parameter and the species, 
    #        array[2] is the environmental parameters, 
    #        array[3] is the samples
    # interaction 数组的维度为 (nEnv + nSpe, nEnv, nsample-1)，即包含所有物种和环境参数的交互作用
    interaction <- array (0, dim=c(nEnv+nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y),colnames(DY)),colnames(DY),rownames(Y[1:(nsample-1),])))
    # Spe_i <- Y[,i]    ### abundance of species i
    for (n in 1: (nsample-1)){
      # change of the DY
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      # change of  abundance
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      # change of environmental parameters
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      regdata <- data.frame(cbind(deltaSpe, deltaEnv))
      # linear regression on   deltaSpe ~ deltaEnv
      # 使用鲁棒回归来拟合模型 DY ~Spe +Env，并将系数存储到 interaction 中，进行交互作用的计算
      fmal<- as.formula (paste("deltaDY ~ -1 + ",paste(colnames(regdata),collapse="+")))
      fit <- rlm (fmal,  regdata[1:(nsample-1),])
      interaction [,,n] <- fit$coefficients/Y[n,i]
      if(Y[n,i]==0){
        interaction [,,n]<- NA
      } 
      missspecies <- which(Y[n,]==0,useNames = TRUE)
      interaction [missspecies,,n] <-NA
    }
  }else if (nsample >= nSpe){
    # 该模型只考虑物种丰度的影响，采用泰勒展开
    # array: array[1] is the species, 
    #        array[2] is the environmental parameters, 
    #        array[3] is the samples
    interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
    for (n in 1: (nsample-1)){
      # change of the DY
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      # change of  abundance
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      # change of environmental parameters
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      regdata <- data.frame(cbind(deltaSpe))
      # linear regression on   deltaSpe ~ deltaEnv
      # 使用 普通线性回归（lm）来估计物种丰度与环境参数变化之间的关系
      fmal<- as.formula (paste("deltaDY ~ -1 + ",paste(colnames(regdata),collapse="+")))
      fit <- lm (fmal,  regdata[1:(nsample-1),])
      interaction [,,n] <- fit$coefficients/Y[n,i]
      if(Y[n,i]==0){
        interaction [,,n]<- NA
      } 
      missspecies <- which(Y[n,]==0,useNames = TRUE)
      interaction [missspecies,,n] <-NA
    }
  }else{
    # the model use the low precision method
    # array: array[1] is the species, 
    #        array[2] is the environmental parameters, 
    #        array[3] is the samples
    interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
    for (n in 1: (nsample-1)){
      
      # change of the DY
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      # change of  abundance
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      # change of environmental parameters
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      temp<-apply(deltaSpe[1:(nsample-1),],c(1,2),function(x) 1/(x*nsample))
      temp[which(!is.finite(temp))] <- 0
      interaction [,,n] <-t(temp) %*% deltaDY    # this is equivalent to calculate the mean values
      if(Y[n,i]==0){
        interaction [,,n]<- NA
      } 
      missspecies <- which(Y[n,]==0,useNames = TRUE)
      interaction [missspecies,,n] <-NA  
    }
  }
  return(interaction)
}


interinf02 <- function(Y, X, DY, i){
  #### 目的：计算物种与环境之间的交互作用（使用线性回归分析环境变化对物种丰度变化的影响），并返回一个多维数组表示这些交互作用。 
  #### 输入参数：Y: 物种丰度数据（一个矩阵，行是样本，列是物种）。 
  ####           X: 环境因子数据（一个矩阵，行是样本，列是环境因子）。
  ####           DY: 物种丰度变化数据（每个物种在不同样本间的变化）。 
  ####           i: 当前物种（列索引）。 
  
  # 初始化变量
  nsample <- nrow(DY)   # 样本数量
  nEnv <- ncol(DY)      # 环境因子数量
  nSpe <- ncol(Y)       # 物种数量
  
  # 创建 interaction 数组：
  #interaction 数组是一个三维数组，用来存储每个物种、每个环境因子、每个样本的交互作用值（这里减去1是因为计算样本之间的变化）。
  interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
  
  #循环计算物种和环境的交互作用
  for (n in 1: (nsample-1)){
    
    # 计算 DY( 当前样本与其他样本之间环境因子的变化) 和 Spe(当前样本与其他样本之间物种丰度的变化)：
    deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
    
    #deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
    
    regdata <- data.frame(cbind(deltaSpe))
    for(k in 1:ncol(deltaDY)){
      deltaDY_k <- deltaDY[,k]
      
      
      # 线性回归模型
      # 通过 DY_k 与 Spe 的线性回归，分析环境因子的变化对物种丰度变化的影响
      fmal<- as.formula (paste("deltaDY_k ~ -1 + ",paste(colnames(regdata),collapse="+")))
      # 该模型使用了 rlm() 函数（鲁棒回归）来拟合数据，acc用来放宽收敛标准
      fit <- rlm (fmal,  regdata[1:(nsample-1),], maxit=500,acc=1e-3)
      # 存储交互作用系数：
      interaction [,k,n] <- fit$coefficients/Y[n,i]
    }
    
    # 处理零丰度情况
    # 如果当前物种在当前样本的丰度为零，则将该交互作用值设为 NA
    if(Y[n,i]==0){
      interaction [,,n]<- NA
    } 
    missspecies <- which(Y[n,]==0,useNames = TRUE)
    interaction [missspecies,,n] <-NA
  }
  #返回一个多维数组interaction，其形状为 nSpe x nEnv x (nsample - 1)，表示每个物种、每个环境因子、每个样本的交互作用系数。
  return(interaction)
}



interinf03 <- function(Y, X, DY, i){
  #### 目的： 通过线性回归计算物种丰度变化与环境因子变化之间的交互作用
  #### 输出： 一个三维数组，其中包含每个物种、每个环境因子与每个样本的交互作用系数
  #### 输入：
  ####       Y: 物种丰度矩阵（样本 × 物种）
  ####       X: 环境因子矩阵（样本 × 环境因子）
  ####       DY: 物种丰度变化矩阵（样本 × 物种的变化量）
  
  # 初始化变量：
  nsample <- nrow(DY)    # 样本数量
  nEnv <- ncol(DY)       # 样本数量
  nSpe <- ncol(Y)        # 物种数量（Y 的列数）
  
  # 创建一个三维数组 interaction 来存储每个物种、每个环境因子与每个样本的交互作用系数：
  interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
  
  # 循环遍历每个样本，（但排除了最后一个样本），目的是计算相邻样本间的变化。
  for (n in 1: (nsample-1)){
    
    # 计算物种和环境因子的变化：
    # DY: 当前样本与其他样本之间环境因子变化量
    deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
    # Spe: 当前样本与其他样本之间物种丰度变化量
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
    # Env: 当前样本与其他样本之间环境因子变化量
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
    
    # 线性回归模型   deltaSpe ~ deltaEnv
    # 创建一个数据框，包含物种丰度变化 Spe（自变量）
    regdata <- data.frame(cbind(deltaSpe))
    # fmal：构建回归公式，用环境因子变化 DY 作为因变量，物种丰度变化 Spe 作为自变量
    fmal<- as.formula (paste("deltaDY ~ -1 + ",paste(colnames(regdata),collapse="+")))
    # lm(fmal, ...)：用普通最小二乘法（OLS）进行线性回归分析，拟合 DY ~Spe 的模型
    fit <- lm (fmal,  regdata[1:(nsample-1),])
    # 存储回归系数
    interaction [,,n] <- fit$coefficients
    # 处理零丰度情况
    if(Y[n,i]==0){
      interaction [,,n]<- NA
    } 
    # 处理缺失物种
    missspecies <- which(Y[n,]==0,useNames = TRUE)
    interaction [missspecies,,n] <-NA
  }
  # 最终返回一个三维数组 interaction，它包含了每个物种、每个环境因子以及每个样本的交互作用系数
  return(interaction)
}



interinf_low <- function(Y, X, DY, i){
  #### 目的:计算物种丰度变化与环境因子变化之间的交互作用
  ####      与 interinf03 和 interinf02 类似，但使用了一种低精度的方法来估算交互作用
  #### 输入： 
  ####      Y: 物种丰度矩阵，维度为 样本 × 物种
  ####      X: 环境因子矩阵，维度为 样本 × 环境因子
  ####      DY: 物种丰度变化矩阵，维度为 样本 × 环境因子，表示每个样本的物种丰度变化
  ####      i: 当前分析的物种（列索引）。 
  # 初始化变量
  nsample <- nrow(DY)   # 样本数量
  nEnv <- ncol(DY)      # 环境因子数量
  nSpe <- ncol(Y)       # 物种数量
  #创建一个三维数组 interaction，用来存储每个物种、环境因子和样本之间的交互作用系数：
  interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
  # 循环遍历每个样本
  for (n in 1: (nsample-1)){
    #计算环境因子的变化(丰度对环境因子的影响)：
    deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
    #计算物种丰度变化：
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
    # 计算环境因子值的变化，主要计算环境因子本身在不同样本之间的变化
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
    
    # 低精度方法计算交互作用：
    temp<-apply(deltaSpe[1:(nsample-1),],c(1,2),function(x) 1/(x*nsample))
    temp[which(!is.finite(temp))] <- 0
    interaction [,,n] <-t(temp) %*% deltaDY    ###### this is equivalent to calculate the mean values
    
    if(Y[n,i]==0){
      interaction [,,n]<- NA
    } 
    
    missspecies <- which(Y[n,]==0,useNames = TRUE)
    interaction [missspecies,,n] <-NA  
  }
  return(interaction)
}


#################################################################
######            part C：计算全局交互矩阵                   ####
#################################################################

interMatrix <- function(Y,X,type){
  #### this function is used to calculate the global interaction matrix between species
  #### there are three methods to fix the global interaction value
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  ####   Y is the abundance 
  ####   X is the environmental parameter
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  results <-  data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  for (species_i in 1:nSpe){
    for (species_j in 1:nSpe){
      derEnv <- rate_change (Y, X, species_i)
      derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
      interaction <- interinf (Y, X[,colnames(derEnv)], derEnv, species_i)     
      interinf_ij <- interaction[species_j,,]
      if (type == 0){
        results[species_i,species_j] <- mean(interinf_ij,na.rm = TRUE)
      }
      if (type == 1){
        results[species_i,species_j] <-median(interinf_ij,na.rm = TRUE)
      }
      if (type == 2){
        paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
        results[species_i,species_j]<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
      }
    }
  }
  return(results) 
}



  #### 目的：用于计算物种间的 全局交互矩阵
  ####       根据不同的统计方法来确定物种间交互的强度（如均值、中位数或加权均值）
  #### 参数： 
  ####      Y：物种丰度矩阵，维度为 nsample × nSpe，表示每个样本中各物种的丰度
  ####      X：环境参数矩阵，维度为 nsample × nEnv，表示每个样本中不同环境因子的值
  #### type：指定计算交互值的方法：
  ####      type = 0：使用均值
  ####      type = 1：使用中位数
  ####      type = 2：使用加权均值
  #### count, integer number: 这用于根据最小值和最大值计算计数段中的直方图分布
  #### factor02, a factor used to compare the negative and positive parts
  #### threshold01: is the parameter used in function summ_inter_ij_new
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  interMatrix02 <- function(Y, X, type, count, factor02, speciesprecision, envprecision, threshold01, returnPerEnv = TRUE) {
    nsample <- nrow(Y)
    nEnv <- ncol(X)
    nSpe <- ncol(Y)
    results <- as.data.frame(matrix(vector(), nSpe, nSpe, dimnames = list(colnames(Y), colnames(Y))), stringsAsFactors = F)
    
    # 新增：存储每个环境因子的矩阵列表
    if (returnPerEnv) {
      perEnvResults <- vector("list", nEnv)
      names(perEnvResults) <- colnames(X)
      # 初始化每个环境因子的矩阵
      for (env in colnames(X)) {
        perEnvResults[[env]] <- matrix(NA, nrow = nSpe, ncol = nSpe, dimnames = list(colnames(Y), colnames(Y)))
      }
    }
    
    pb <- txtProgressBar(min = 0, max = nSpe, style = 3, char = "=")
    cat("\nCalculating interaction matrix:\n")
    
    for (species_i in 1:nSpe) {
      if (envprecision == "high") {
        derEnv <- rate_change02(Y, X, species_i)
      } else {
        derEnv <- rate_change_low(Y, X, species_i)
      }
      derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
      
      if (speciesprecision == "high") {
        interaction <- interinf02(Y, X[, colnames(derEnv)], derEnv, species_i)
      } else {
        interaction <- interinf_low(Y, X[, colnames(derEnv)], derEnv, species_i)
      }
      
      # 新增：为当前物种i初始化每个环境因子的矩阵（如果之前未初始化）
      if (returnPerEnv) {
        for (env in colnames(X)) {
          if (!env %in% names(perEnvResults)) {
            perEnvResults[[env]] <- matrix(NA, nrow = nSpe, ncol = nSpe, dimnames = list(colnames(Y), colnames(Y)))
          }
        }
      }
      
      for (species_j in 1:nSpe) {
        interinf_ij <- interaction[species_j, , ]
        
        # 新增：保存当前物种对(i,j)在每个环境因子的值
        if (returnPerEnv) {
          for (env in rownames(interinf_ij)) {
            if (env %in% colnames(X)) {
              perEnvResults[[env]][species_i, species_j] <- mean(interinf_ij[env, ], na.rm = TRUE)
            }
          }
        }
        
        if (type == 0) {
          results[species_i, species_j] <- mean(interinf_ij, na.rm = TRUE)
        } else if (type == 1) {
          results[species_i, species_j] <- median(interinf_ij, na.rm = TRUE)
        } else if (type == 2) {
          paratrans <- t(X)[rownames(interinf_ij), colnames(interinf_ij)]
          results[species_i, species_j] <- median(colSums(as.matrix(paratrans) * as.matrix(interinf_ij)), na.rm = TRUE)
        } else if (type == 3) {
          Summ_interinf_ij <- summ_inter_ij_new(interinf_ij, count, 0, threshold01)
          number_type_no <- length((which(Summ_interinf_ij$intertype == "no")))
          
          if (number_type_no > 5) {
            results[species_i, species_j] <- NA
          } else {
            Posimedian <- ifelse(
              length(which(Summ_interinf_ij$intertype == "+")) > 0,
              median(Summ_interinf_ij[which(Summ_interinf_ij$intertype == "+"), c("peakvalue")], na.rm = TRUE),
              NA
            )
            Negamedian <- ifelse(
              length(which(Summ_interinf_ij$intertype == "-")) > 0,
              median(Summ_interinf_ij[which(Summ_interinf_ij$intertype == "-"), c("peakvalue")], na.rm = TRUE),
              NA
            )
            
            if (is.na(Posimedian)) Posimedian <- 0
            if (is.na(Negamedian)) Negamedian <- 0
            
            posipara <- length(which(Summ_interinf_ij$intertype == "+"))
            negapara <- length(which(Summ_interinf_ij$intertype == "-"))
            
            if ((posipara == 0) | (negapara > 0.8 * nEnv)) {
              results[species_i, species_j] <- Negamedian
            } else if ((negapara == 0) | (posipara > 0.8 * nEnv)) {
              results[species_i, species_j] <- Posimedian
            } else if (Posimedian > (factor02 * abs(Negamedian))) {
              results[species_i, species_j] <- Posimedian
            } else if ((factor02 * Posimedian) < (abs(Negamedian))) {
              results[species_i, species_j] <- Negamedian
            } else {
              results[species_i, species_j] <- 0
            }
          }
        }
      }
      setTxtProgressBar(pb, species_i)
    }
    close(pb)
    cat("\nInteraction matrix calculation completed.\n")
    
    # 根据 returnPerEnv 决定返回值
    if (returnPerEnv) {
      return(list(globalMatrix = results, perEnvMatrices = perEnvResults))
    } else {
      return(results)
    }
  }


interMatrix_norm <- function(Y,X,speciesprecision,envprecision){
  #### 目的：根据物种丰度 (Y) 和环境参数 (X)，使用不同的精度设置计算物种之间的交互作用强度
  #### 参数：
  ####    Y：物种丰度矩阵（通常是一个物种-样本矩阵），行表示样本，列表示物种
  ####    X：环境参数矩阵，行表示样本，列表示环境参数
  ####    speciesprecision：物种精度
  ####    envprecision：环境精度
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  nsample <- nrow(Y)   # 获取样本数
  nEnv <- ncol(X)      # 获取环境参数的数量
  nSpe <- ncol(Y)      # 获取物种数量
  # 创建一个空的结果矩阵 nSpe x nSpe，行列分别为物种名称
  results <-  as.data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  # 外层循环是遍历所有物种 (species_i)，计算物种 i 与其他物种之间的交互作用
  for (species_i in 1:nSpe){
    if (envprecision=="high"){
      derEnv <- rate_change02 (Y, X, species_i)   # 计算环境变化率
    }else {
      derEnv <- rate_change_low (Y, X, species_i) # 计算环境变化率
    }
    # 从 derEnv 中去除含有 NA 的列。使用 apply 函数遍历每列，检查是否有 NA，并保留没有 NA 的列
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    
    if (speciesprecision=="high"){
      interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i)   # 计算物种交互作用
    }else{
      interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i) # 计算物种交互作用
      
    }
    # 计算交互作用矩阵并标准化：
    # 内层循环遍历所有物种 species_j，从 interaction 中获取物种 i 和 j 的交互作用值
    for (species_j in 1:nSpe){   
      interinf_ij <- interaction[species_j,,]
      temp<-interinf_ij
      temp[is.na(temp)] <- 0                       # 将 NA 值替换为 0
      results[species_i,species_j] <- norm(temp)   # 对交互作用进行标准化并存储
    }
  }    
  return(results) 
}


interMatrix_norm_summ <- function(Y,X,count,speciesprecision,envprecision){
  #### 目的：与前面的 interMatrix_norm 类似，但它在计算全局交互作用矩阵时做了一些修改
  ####     通过定义一个直方图分布来排除极值，并计算物种间的交互作用矩阵
  #### 参数：
  ####     Y：物种丰度矩阵。行表示样本，列表示物种
  ####     X：环境参数矩阵。行表示样本，列表示环境参数
  ####     count：用于计算直方图的段数。通过这个参数，代码会将交互作用值划分成多个区间来计算其分布
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  
  nsample <- nrow(Y)   # 样本的数量（即 Y 的行数）
  nEnv <- ncol(X)      # 环境参数的数量（即 X 的列数）
  nSpe <- ncol(Y)      # 物种的数量（即 Y 的列数）
  # 创建一个 nSpe x nSpe 的空矩阵，存储物种间的交互作用值，行列名为 Y 中的物种名称
  results <-  as.data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  # 外层循环遍历所有物种 species_i，计算物种 i 与其他物种之间的交互作用
  for (species_i in 1:nSpe){
    if (envprecision=="high"){
      derEnv <- rate_change02 (Y, X, species_i)   # 计算环境变化率
    }else {
      derEnv <- rate_change_low (Y, X, species_i) # 计算环境变化率
    }
    # 去除 derEnv 中包含 NA 的列。这是通过 apply 函数遍历每列，检查是否含有 NA，然后只保留没有 NA 的列
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    
    if (speciesprecision=="high"){
      interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i)   # 计算物种交互作用
    }else{
      interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i) # 计算物种交互作用
    }
    # 计算交互作用矩阵并处理极值
    # 内层循环遍历所有物种 species_j，从 interaction 中提取物种 i 和 j 的交互作用数据
    for (species_j in 1:nSpe){   
      interinf_ij <- interaction[species_j,,]
      temp<-c()   # 初始化一个空向量，存储计算结果
      for (para in 1:nrow(interinf_ij)){   # 遍历每一行交互作用数据
        a <-interinf_ij[para,]             # 取出物种 j 和所有样本的交互作用数据
        # 计算交互作用数据的直方图分布
        # cut 函数将交互作用数据划分为 count 个区间，计算每个区间的频率分布，生成直方图
        counttable<-table(cut(interinf_ij[para,],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        # which.max 查找直方图中最大频率对应的位置，这个位置对应着数据的峰值
        peakposition <- which.max(counttable)
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        # 提取位于峰值区间内的数据 (peakvalue)，并对其进行二范数规范化（norm(peakvalue, type="2")）
        # 将规范化后的结果存入 temp 向量
        peakvalue <- a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])]
        temp<- c(temp,norm(peakvalue,type="2"))
      }
      # 最后，对所有峰值的规范化结果再次进行二范数规范化，并存入 results[species_i, species_j]
      results[species_i,species_j]<- norm(temp,type="2")
    }    
  }    
  return(results)   # results 包含了物种间交互作用的规范化值
}


interMatrix_low <- function(Y,X,type,count,factor02,threshold01){
  #### 目的：基于环境变化和物种间交互作用，并通过不同的方法来处理交互作用的值
  #### 参数：
  ####     Y：物种丰度矩阵（样本 x 物种）。行表示样本，列表示物种
  ####     X：环境参数矩阵（样本 x 环境参数）。行表示样本，列表示环境参数
  ####     type：决定交互作用计算方法的参数。其值决定了使用哪种方法来汇总交互作用（如均值、中位数、加权和等）
  ####     factor02：一个常数因子，在 type == 3 的情况下，用于比较正负交互作用的中位数
  ####     threshold01：用于 summ_inter_ij_new 函数的阈值参数，影响交互作用的分组
  nsample <- nrow(Y)   # nsample：样本的数量（即 Y 的行数）
  nEnv <- ncol(X)      # nEnv：环境参数的数量（即 X 的列数）
  nSpe <- ncol(Y)      # nSpe：物种的数量（即 Y 的列数）
  # 创建一个 nSpe x nSpe 的空矩阵，存储物种间的交互作用值，行列名为 Y 中的物种名称
  results <-  as.data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  # 计算交互作用
  for (species_i in 1:nSpe){
    # 计算环境变化率
    derEnv <- rate_change_low (Y, X, species_i)
    # 去除包含NA的列
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    # 计算物种间的交互作用
    interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i) 
    # 根据 type 计算交互作用：
    for (species_j in 1:nSpe){   
      interinf_ij <- interaction[species_j,,]   #提取物种 j 的交互作用数据
      if (type == 0){   # 使用交互作用数据的均值
        results[species_i,species_j] <- mean(interinf_ij,na.rm = TRUE)
      }
      if (type == 1){   # 使用交互作用数据的中位数
        results[species_i,species_j] <-median(interinf_ij,na.rm = TRUE)
      }
      if (type == 2){   # 计算加权和的中位数，其中加权系数是环境参数 X 的转置
        paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
        results[species_i,species_j]<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
      }
      if (type == 3){   # 特殊处理：通过更复杂的逻辑，结合正负交互作用的中位数和阈值来决定交互作用值
        #factor02<-5
        # 使用summ_inter_ij_new函数来处理交互作用数据
        Summ_interinf_ij<- summ_inter_ij_new(interinf_ij,count,0,threshold01)
        number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
        # 如果"no"类型的交互作用超过5个，则设置为NA
        if (number_type_no>5){
          results[species_i,species_j] <- NA
        }else{
          # 计算正负交互作用的中位数
          Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
          Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
          # 计算正负交互作用的数量
          posipara<-length(which(Summ_interinf_ij$intertype=="+"))
          negapara<-length(which(Summ_interinf_ij$intertype=="-"))
          
          if(posipara==0){
            results[species_i,species_j] <-  Negamedian   # 只有负交互作用时，结果为负中位数
          }else if (negapara==0){
            results[species_i,species_j] <- Posimedian    # 只有正交互作用时，结果为正中位数
          }else if (Posimedian>(factor02*abs(Negamedian))){
            results[species_i,species_j]<- Posimedian     # 如果正交互作用大于负交互作用的阈值，使用正中位数
          }else if ((factor02*Posimedian)<(abs(Negamedian))){
            results[species_i,species_j] <- Negamedian    # 如果负交互作用大于正交互作用的阈值，使用负中位数
          }else {
            results[species_i,species_j] <- 0             # 否则设为0
          }           
        }         
      }
    }    
  }
  
  return(results)   # 函数返回计算得到的交互作用矩阵 results，其中包含了物种间交互作用的值
}


interMatrix_ij <- function(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01){
  #### 目的：用于计算物种之间的单一交互作用值（即 β_ij，表示物种 i 对物种 j 的影响），
  ####       它是构建整个物种交互作用矩阵（interaction matrix）中某一元素的一部分
  #### 参数：
  ####     Y ：物种丰度矩阵（样本 x 物种）
  ####     X ：环境因子矩阵（样本 x 环境变量）
  ####     type	：计算交互作用方式（0=均值，1=中位数，2=加权和，3=采用工作流法进行统计估计）
  ####     count ：整数：用于根据最小值和最大值计算计数段中的直方图分布 
  ####     factor02 ：正负交互作用比较的调节因子
  ####     threshold01: summ_inter_ij_new函数中使用的参数
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  # 环境变化率计算
  if (envprecision=="high"){
    derEnv <- rate_change02 (Y, X, species_i)
  }else {
    derEnv <- rate_change_low (Y, X, species_i)
  }
  # 去除含 NA 的列，确保数据完整性
  derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
  # 计算物种间交互作用矩阵（对所有物种）
  if (speciesprecision=="high"){
    interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i) 
  }else{
    interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)    
  }
  # 提取物种 i 对物种 j 的交互值数据，是一个矩阵或向量，代表不同环境或时间条件下的影响值
  interinf_ij <- interaction[species_j,,]
  if (type == 0){   # 求 均值
    results <- mean(interinf_ij,na.rm = TRUE)
  }
  if (type == 1){   # 求 中位数
    results <-median(interinf_ij,na.rm = TRUE)
  } 
  if (type == 2){   # 加权求和（考虑环境参数）
    paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
    results<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
  }
  if (type == 3){   # 采用工作流法进行统计估计
    # 调用 summ_inter_ij_new，将交互作用值分类为："+"：正交互作用；"-"：负交互作用"no"；：不显著交互作用
    Summ_interinf_ij<- summ_inter_ij_new(interinf_ij,count,0,threshold01)
    number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
    # 判断“无效交互”是否太多
    if (number_type_no>5){
      results<- NA   # 如果 no 的数量超过 5，认为数据太弱，不可靠，返回 NA
    }else{
      # 否则分析正负交互作用：
      Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
      Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
      # 统计正负数量
      posipara<-length(which(Summ_interinf_ij$intertype=="+"))
      negapara<-length(which(Summ_interinf_ij$intertype=="-"))
      # 比较正负交互的中位数，并结合 factor02：
      # 目的：根据正/负交互强度决定最终的交互作用方向和大小
      if(posipara==0){
        results <-  Negamedian
      }else if (negapara==0){
        results <- Posimedian
      }else if (Posimedian>(factor02*abs(Negamedian))){
        results<- Posimedian
      }else if ((factor02*Posimedian)<(abs(Negamedian))){
        results <- Negamedian
      }else {
        results <- 0
      }
    }
  }
  return(results) 
}


interMatrix_ij_low <- function(Y,X,species_i,species_j,type,count,factor02,threshold01){
  # 目的：计算物种 i 对物种 j 的交互作用值，类似前面的 interMatrix_ij 函数，
  #       其用低精度的环境变化率计算和交互作用推断方法，利用了低精度函数（如 rate_change_low 和 interinf_low）来估算交互作用
  # 参数：同上
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  # 计算环境变化率
  derEnv <- rate_change_low (Y, X, species_i)
  # 去除包含 NA 的列
  derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
  # 计算交互作用
  interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)
  # 提取物种 j 的交互作用数据
  interinf_ij <- interaction[species_j,,]
  if (type == 0){
    results <- mean(interinf_ij,na.rm = TRUE)
  }
  if (type == 1){
    results <-median(interinf_ij,na.rm = TRUE)
  }
  if (type == 2){
    paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
    results<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
  }
  if (type==3){   # 调用 summ_inter_ij_new 函数，对交互作用值进行进一步的统计推断
    Summ_interinf_ij<- summ_inter_ij_new(interinf_ij,count,0,threshold01)
    number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
    if (number_type_no>5){
      results<- NA
    }else{
      Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
      Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
      
      posipara<-length(which(Summ_interinf_ij$intertype=="+"))
      negapara<-length(which(Summ_interinf_ij$intertype=="-"))
      # 根据正负交互作用的数量和大小做出决策
      if(posipara==0){
        results <-  Negamedian
      }else if (negapara==0){
        results <- Posimedian
      }else if (Posimedian>(factor02*abs(Negamedian))){
        results<- Posimedian
      }else if ((factor02*Posimedian)<(abs(Negamedian))){
        results <- Negamedian
      }else {
        results <- 0
      }
    }
  }
  return(results) 
}


################################################################# 
######         part D：交互结果的统计汇总结果              ######
#################################################################
# 在计算出原始的四维交互水平值后，我们需要做一些统计汇总来确定全局交互矩阵，
# 根据模式，我们可以确定全局交互是正的、负的，还是没有明确的模式

summ_inter_ij <- function(Y,X,species_i,species_j,count,type,speciesprecision,envprecision){
  #### 目的：该函数用于计算原始交互值\beta ^k＿ij的统计汇总\alpha
  #### 参数：
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belong to the types of "posi", "nega" or "no"
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  # 计算环境变化率（derEnv）
  if (envprecision=="high"){
    derEnv <- rate_change02 (Y, X, species_i)
  }else {
    derEnv <- rate_change_low (Y, X, species_i)
  }
  # 去除包含缺失值的环境变化率
  derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
  # 计算物种交互作用
  if (speciesprecision=="high"){
    interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i) 
  }else{
    interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)
  }
  # 提取物种 i 对物种 j 的交互作用数据，这是一个三维数组
  # 创建统计结果的列名
  interinf_ij <- interaction[species_j,,]
  b<- as.character(c(1:count))
  b<- t(as.matrix(as.character(c(1:(count-1)))))
  b<- apply(b,1,function(x) paste("count",x,sep=""))
  # 创建一个包含统计结果的列名列表，统计结果包括：
  # intertype：交互作用类型（正、负、无）。
  # min, 1stQu, median, mean, 3rdQu, max：数据的基本统计量。
  # peakrange：交互作用数据的峰值范围。
  # peakvalue：峰值的中位数。
  # variance：数据的方差。
  # numpositive, numnegative：正交互作用和负交互作用的数量。
  # count1, count2, ..., countN：根据分段数量（count）计算的直方图的各个区间。
  namelist <- c("intertype","min","1stQu","median","mean","3rdQu","max", "peakrange", "peakvalue", "variance",  "numpositive","numnegative",b)
  # 计算交互作用的统计特性（type == 0，按样本计算）
  if (type==0){
    sta_inter <-  data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
    for (para in 1:nrow(interinf_ij)){
      a <-interinf_ij[para,]
      temp <- as.numeric(summary(a))[1:6]   # summary(a)：计算交互作用数据的基本统计量
      variance <- var(a,na.rm = TRUE)       # var(a)：计算交互作用数据的方差
      nposi <- sum(  a > 0,na.rm = TRUE )   # nposi：计算正交互作用的数量
      nnega <- nsample-nposi-sum(is.na(a))  # nnega：计算负交互作用的数量
      n_non_NA<-nsample- sum(is.na(a))
      if (nposi/nsample>0.8) {
        # 正交互作用占比超过 80%，则判定为正交互作用（+）
        intertype <- "+"    
      }else if (nnega/nsample>0.8){
        # 负交互作用占比超过 80%，则判定为负交互作用（-）
        intertype <-"-"
      }else {
        # 否则，判定为无明显交互作用（no）
        intertype <-"no"
      }
      
      if (n_non_NA== 1 | n_non_NA== 0){
        sta_inter[para,1:6]
      }else{
        # 通过 cut 函数根据交互作用值的范围，将其划分为 count 个区间
        counttable<-table(cut(interinf_ij[para,],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        peakposition <- which.max(counttable)   # peakposition 表示交互作用数据的峰值所在位置
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        peakvalue <- median(a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])])
        sta_inter[para,c(-8,-1)]<-c(temp,peakvalue,variance,nposi,nnega,as.numeric(counttable))
        sta_inter[para,1]<-intertype
        sta_inter[para,8]<-peakrange
      }
    }
  }else {
    # 计算交互作用的统计特性（type != 0，按环境参数计算）           
    sta_inter <-  data.frame(matrix(vector(), ncol(interinf_ij), length(namelist),dimnames=list(colnames(interinf_ij), namelist)))
    for (n in 1:nrow(interinf_ij)){
      a <-interinf_ij[,n]
      temp <- as.numeric(summary(a))[1:6]
      variance <- var(a,na.rm = TRUE) 
      nposi <- sum(  a > 0 ,na.rm = TRUE )
      nnega <- nsample-2-nposi-sum(is.na(a))
      n_non_NA<-nsample- sum(is.na(a))
      if (n_non_NA== 1 | n_non_NA== 0){
        sta_inter[para,1:6]
      }else{
        counttable<-table(cut(interinf_ij[,n],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        peakposition <- which.max(counttable)
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        peakvalue <- median(a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])])
        sta_inter[n,-8]<-c(temp,peakvalue,variance,nposi,nnega,as.numeric(counttable))
        sta_inter[n,8]<-peakrange
      }
    }
  }
  return(sta_inter)
}


summ_inter_ij_new <- function(interinf_ij,count,type,threshold01){
  #### 该函数用于计算原始交互值\beta ^k＿ij的统计汇总\alpha
  ####  interinf_ij is the interaction data frame of  beta_ij
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  #### threshold01: 用于判断交互作用的阈值，决定交互作用是正向、负向还是无明显模式。
  # 初始化和列名设置
  nsample <- ncol(interinf_ij)   # nsample：获取样本数量，即 interinf_ij 数据框的列数
  # b：生成包含 count 个区间名称的向量，用于描述交互作用数据的各个分段区间
  b<- as.character(c(1:count))
  b<- t(as.matrix(as.character(c(1:(count-1)))))
  b<- apply(b,1,function(x) paste("count",x,sep=""))
  #namelist：设置统计结果的列名，包括基本统计量（如最小值、四分位数、中位数等），
  # 峰值信息（peakrange 和 peakvalue），交互作用类型（intertype）以及分段统计数据
  namelist <- c("intertype","min","1stQu","median","mean","3rdQu","max", "peakrange", "peakvalue", "variance",  "numpositive","numnegative",b)
  # 统计摘要计算：按样本（type = 0）
  if (type==0){
    sta_inter <-  data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
    for (para in 1:nrow(interinf_ij)){
      a <-interinf_ij[para,]
      temp <- as.numeric(summary(a))[1:6]   # summary(a)计算数据的基本统计量
      variance <- var(a,na.rm = TRUE)       # var(a)：计算数据的方差
      nposi <- sum(  a > 0,na.rm = TRUE )   # nposi：正值的个数，表示正向交互作用
      nnega <- nsample-nposi-sum(is.na(a))  # nnega：负值的个数，表示负向交互作用
      n_non_NA<-nsample- sum(is.na(a))      # n_non_NA：非缺失值的样本数
      if (nposi/nsample>threshold01) {
        # 如果正向交互作用占比大于 threshold01，则标记为正向（+）
        intertype <- "+"
      }else if (nnega/nsample>threshold01){
        # 如果负向交互作用占比大于 threshold01，则标记为负向（-）
        intertype <-"-"
      }else {
        # 否则，标记为没有明确模式（no）
        intertype <-"no"
      }
      # 数据分段和峰值计算
      if (n_non_NA== 1 | n_non_NA== 0){
        sta_inter[para,1:6]
      }else{
        counttable<-table(cut(interinf_ij[para,],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        peakposition <- which.max(counttable)
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        peakvalue <- median(a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])])
        sta_inter[para,c(-8,-1)]<-c(temp,peakvalue,variance,nposi,nnega,as.numeric(counttable))
        sta_inter[para,1]<-intertype
        sta_inter[para,8]<-peakrange
      }
    }
  }else {
    # 统计摘要计算：按环境参数（type != 0）                
    sta_inter <-  data.frame(matrix(vector(), ncol(interinf_ij), length(namelist),dimnames=list(colnames(interinf_ij), namelist)))
    for (n in 1:nrow(interinf_ij)){
      a <-interinf_ij[,n]
      temp <- as.numeric(summary(a))[1:6]
      variance <- var(a,na.rm = TRUE) 
      nposi <- sum(  a > 0 ,na.rm = TRUE )
      nnega <- nsample-2-nposi-sum(is.na(a))
      n_non_NA<-nsample- sum(is.na(a))
      if (n_non_NA== 1 | n_non_NA== 0){
        sta_inter[para,1:6]
      }else{
        counttable<-table(cut(interinf_ij[,n],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        peakposition <- which.max(counttable)
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        peakvalue <- median(a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])])
        
        sta_inter[n,-8]<-c(temp,peakvalue,variance,nposi,nnega,as.numeric(counttable))
        sta_inter[n,8]<-peakrange
      }
    }
  }
  return(sta_inter) # 包含了所有计算的统计信息，如交互作用类型、峰值范围、方差、正负交互作用的数量等
}


summ_interaction <-function(Y,X,count,type){
  # 目的：计算物种间的交互作用统计量
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  # b：创建一个字符串向量，表示交互作用的分段
  b<- t(as.matrix(as.character(c(1:(count-1)))))
  b<- apply(b,1,function(x) paste("count",x,sep=""))
  # namelist：包含计算的统计量的名称，如最小值、四分位数、中位数、均值、方差、皮尔逊相关系数、斯皮尔曼相关系数等
  namelist <- c("min","1stQu","median","mean","3rdQu","max","variance", "PearsonCorr", "PearsonPvalue","SpearmanCorr","SpearmanPvalue",  "numpositive","numnegative",b)
  # type == 0，计算的统计量是按物种的交互作用
  if (type==0){
    # sta_inter：初始化一个四维数组，用来存储每一对物种间（species_i 和 species_j）的交互作用统计数据，
    # 维度是 (nSpe, nSpe, nEnv, length(namelist))，对应物种对、环境因子和统计量名称
    sta_inter <-  array(NA, dim=c(nSpe,nSpe, nEnv, length(namelist)),dimnames=list(colnames(Y),colnames(Y),colnames(X),namelist))
    # 计算每对物种的交互作用
    for(species_i in 1:nSpe){
      for(species_j in 1:nSpe){
        # 调用 rate_change 函数，计算物种响应的环境变化
        derEnv <- rate_change (Y, X, species_i)
        derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
        # 通过 interinf 函数计算交互作用的矩阵
        interaction <- interinf (Y, X[,colnames(derEnv)], derEnv, species_i)     
        # 从交互作用矩阵中提取 species_j 对 species_i 的交互作用数据
        interinf_ij <- interaction[species_j,,]
        # 计算统计量：
        for (para in 1:nrow(interinf_ij)){
          a <-interinf_ij[para,]
          temp <- as.numeric(summary(a))[1:6]
          variance <- var(a,na.rm = TRUE) 
          Pear_cor <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a)$estimate
          Pear_p <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a)$p.value
          Spear_cor <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a,method="spearman")$estimate
          Spear_p <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a,method="spearman")$p.value
          nposi <- sum(  a > 0 ,na.rm = TRUE)
          nnega <- nsample-2-nposi-sum(is.na(a))
          counttable<-table(cut(interinf_ij[para,],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
          sta_inter[species_i,species_j,rownames(interinf_ij)[para],]<-c(temp,variance,Pear_cor,Pear_p,Spear_cor,Spear_p,nposi,nnega,counttable)
        }
      }
    }
  }else {
    # type != 0：按样本计算统计量
    sta_inter <-  array(NA, dim=c(nSpe,nSpe, nsample-2, length(namelist)),dimnames=list(colnames(Y),colnames(Y),rownames(Y)[1:(nsample-2)],namelist))
    for(species_i in 1:nSpe){
      for(species_j in 1:nSpe){  
        derEnv <- rate_change (Y, X, species_i)
        derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
        interaction <- interinf (Y, X[,colnames(derEnv)], derEnv, species_i)     
        interinf_ij <- interaction[species_j,,]
        for (n in 1:nrow(interinf_ij)){
          a <-interinf_ij[,n]
          temp <- as.numeric(summary(a))[1:6]
          variance <- var(a,na.rm = TRUE) 
          nposi <- sum(  a > 0,na.rm = TRUE )
          nnega <- nsample-2-nposi-sum(is.na(a))
          counttable<-table(cut(interinf_ij[,n],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
          sta_inter[species_i,species_j,colnames(interinf_ij)[para],]<-c(temp,variance,nposi,nnega,counttable)
        }
      }
    }
  }
  return(sta_inter)
}


anascatter <-function(Y,X,interaction_ij,para,count){
  #### 目的：计算并返回物种交互作用的极值差异以及与环境参数的关系差异
  #### 交互作用interaction_ij为一对物种间的相互作用强度二维矩阵，rownames为环境参数，colnames为样本###
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  position <- which(rownames(interaction_ij)==para)   #查找 para 所在的行号（即该环境参数在交互作用矩阵中的位置）
  window <-5   # window：窗口大小，指的是每次分析时考虑的样本数量（默认值是 5）
  N<- floor(nsample/window)   # N：样本数除以窗口大小，表示一共有多少个窗口
  # 创建名列表和数据框
  b<- t(as.matrix(as.character(c(1:(nsample-1-window)))))
  b<- apply(b,1,function(x) paste("window",x,sep=""))
  namelist <- c("extreme",b)
  difference<-  data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
  paradiffer <- data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
  # 循环计算交互作用差异：
  for (n_para in (1:nrow(interaction_ij))){
    paraorig <-X[colnames(interaction_ij),n_para]
    paraorder<- order(paraorig)
    counttable<-table(cut(paraorig,breaks=seq(min(paraorig,na.rm = TRUE),max(paraorig,na.rm = TRUE),length=count))) #### check the data extrem properties 
    diff_temp <- vector(mode="numeric", length=0)
    paradiff_temp <- vector(mode="numeric", length=0)
    # 如果频数的方差大于 20，则认为数据具有较大的波动，
    # 在这种情况下，程序会重新划分 paraorig 的区间，使得数据集中在最大频数的区间
    if (var(counttable)>20){
      breakpoint <- seq(min(paraorig,na.rm = TRUE),max(paraorig,na.rm = TRUE),length=count)
      pp <- which(counttable == max(counttable,na.rm = TRUE))
      paraorig <- paraorig[which(paraorig >breakpoint[pp] & paraorig <= breakpoint[pp+1])]
      paraorder<- order(paraorig)
    }
    # 计算交互作用的极值差异和环境参数的差异：
    for (k in 1:(nsample-1-window)){
      # section_k：为当前窗口的交互作用数据
      section_k <- interaction_ij[position,paraorder[k:(k+window-1)]]
      # diff_temp：计算每个窗口的交互作用的最大值与最小值的差异
      diff_temp<-c(diff_temp,max(section_k,na.rm = TRUE)-min(section_k,na.rm = TRUE))
      # paradiff_temp：计算每个窗口的环境参数（paraorig）的最大值与最小值的差异
      Parasection_k <- paraorig[paraorder[k:(k+window-1)]]
      paradiff_tem<-c(paradiff_temp,max(Parasection_k,na.rm = TRUE)-min(Parasection_k,na.rm = TRUE))
      #}
      #      else{
      #       section_k <- interaction_ij[position,paraorder[((k-1)*window+1):(nsample-2)]]
      #       diff <-c(diff,max(section_k)-min(section_k))
      #     }
      
    }
    difference[n_para,] <- c(max(counttable),diff_temp)
    paradiffer[n_para,] <- c(max(counttable),paradiff_temp)  
  }
  return(difference)
  #return(cbind(difference,paradiffer))
}


####################################################################
######                part E：计算主要环境参数                ######
####################################################################
main_para<-function(Y,X,type,count,factor02,speciesprecision,envprecision,threshold01){
  ##### this function is used to find the main parameters
  ##### Y species abundance
  ####  X environmental parameters
  
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belong to the types of "posi", "nega" or "no"
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  
  #InterMatrix_Summ<-  data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(RelaAbund), colnames(RelaAbund))),stringsAsFactors=F)
  
  sig_para <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
  
  posi_domi_para <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
  nega_domi_para <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
  
  for (i in 1:nSpe){
    if (envprecision=="high"){
      derEnv <- rate_change02 (Y, X, i)
    }else {
      derEnv <- rate_change_low (Y, X, i)
    }
    
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    
    if (speciesprecision=="high"){
      interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, i) 
    }else{
      interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, i)
      
    }
    
    for(j in 1:nSpe){
      interinf_ij <- interaction[j,,]
      #Summ_interinf_ij <- summ_inter_ij(Y,X,i,j,count,0)
      Summ_interinf_ij <- summ_inter_ij_new(interinf_ij,count,0,threshold01)
      
      number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
      
      temp <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
      temp_domi_posi <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
      temp_domi_nega <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
      
      if (number_type_no>5){
        #InterMatrix_Summ_test[i,j]<- NA
      }else{
        #       x1<-Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")]
        #       x1<-as.matrix(x1)
        if("+" %in% Summ_interinf_ij$intertype ){
          Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
          posipara<-rownames(Summ_interinf_ij)[which(Summ_interinf_ij$intertype=="+")]
          maxi_posi <- max(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
          domi_posi <- posipara[which(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")]>0.7*maxi_posi)]
          temp_domi_posi[1,domi_posi]<-1
        }
        
        if("-" %in% Summ_interinf_ij$intertype ){
          Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
          negapara<-rownames(Summ_interinf_ij)[which(Summ_interinf_ij$intertype=="-")]
          mini_nega <- min(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
          domi_nega <- negapara[which(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")]<0.7*mini_nega)]
          temp_domi_nega[1,domi_nega]<-1
        }
        
        
        
        if ((length(posipara)==nEnv)| (Posimedian>(factor02*abs(Negamedian))) ) {
          
          #InterMatrix_Summ[i,j] <- Posimedian
          
          temp[1,posipara]<-1
        } else if ((length(negapara)==nEnv)|((factor02*Posimedian)<(abs(Negamedian)))){
          #InterMatrix_Summ[i,j] <- Negamedian
          
          temp[1,negapara]<--1
          
        }else {
          #InterMatrix_Summ[i,j] <- 0
        }    
        
      }
      rownames(temp)<- paste("interinf",i,j,sep="_")  
      sig_para<- rbind(sig_para, temp)
      rownames(temp_domi_posi)<- paste("interinf",i,j,sep="_")  
      posi_domi_para<- rbind(posi_domi_para, temp_domi_posi)
      rownames(temp_domi_nega)<- paste("interinf",i,j,sep="_")  
      nega_domi_para<- rbind(nega_domi_para, temp_domi_nega)
      
    }
  }
  out<-list(posi_domi_para,nega_domi_para,sig_para)
  return(out)
}


main_para_barplot<- function(out){
  
  posi_domi_para<-out[[1]]
  nega_domi_para<-out[[2]]
  sig_para<-out[[3]]
  ###### make bar plot #### 
  p1<-apply(posi_domi_para, 2, function(x) sum(x,na.rm=TRUE) )
  p2<-apply(nega_domi_para, 2, function(x) sum(x,na.rm=TRUE) )
  p3<-apply(sig_para, 2, function(x) sum(!is.na(x)) )
  
  infor_domi_para<-as.data.frame(rbind(p1,p2,p3))
  
  rownames(infor_domi_para)<- c("positive","negative","determining interaction")
  
  
  infor_domi_para$statinfor <- rownames(infor_domi_para)
  sss<- melt(infor_domi_para,id="statinfor")
  colnames(sss)[1]<-c("dominant_parameters")
  
  
  term<-ggplot(sss, aes(variable, value))+
    geom_bar(stat = "identity",aes(fill=dominant_parameters),position = "dodge",width=0.5)+
    coord_flip()+xlab("environmental parameters")+ylab("number of frequency")+
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
    theme(axis.title.x = element_text(size = rel(1.5), angle = 0))+
    theme(axis.text.x = element_text(size = rel(1.5), angle = 0))+
    theme(axis.text.y = element_text(size = rel(1.5), angle = 0))+
    theme(legend.title = element_text(size = rel(1.0), angle = 0))+
    theme(legend.text = element_text(size = rel(1.0), angle = 0))+
    theme(legend.position = "bottom")
  
  pdf(paste("Bar_Plot_statistic_Dominant_Parameters.pdf",sep=""),height = 10, width = 10)
  print(term)
  dev.off()
  return(1)
}



#################################################################
######                 part F：robust test                 ######
#################################################################
interMatrix_ij_permut <- function(Y,X,species_i,species_j,times,size,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01){
  
  
  ###### randomly remove  size  species for times, save the 2 dimension interaction table into the data frame #####
  #####  permut_type is "samples",  "species", "para"  
  ##### 
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  ####  time, how many times in  the random test
  #### size,  control the number of movement of species or parameters
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  #### permutation type:  "species", "para", "samples"
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belongs to the types of "posi", "nega" or "no"
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  species_i_name<-colnames(Y)[species_i]
  species_j_name<-colnames(Y)[species_j]
  
  interinf_ij <-  interMatrix_ij(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
  
  
  #interinf_ij_permut <- vector(mode="numeric", length=0)
  testpermut <- c(1:times)
  removenames<-vector()
  for(k in 1:size){
    removenames<- c(removenames,paste(permut_type,k,"remove",sep="_"))
  }
  
  interinf_ij_permut <-  data.frame(matrix(vector(),  times+1,size+1,dimnames=list( c(testpermut,"original"),c("values",removenames))),stringsAsFactors=F)
  
  
  
  if(permut_type=="species"){
    testspecies <- colnames(Y)[c(-species_i,-species_j)]
    
    for (n in 1:times){
      kickoutspe<- sample(testspecies,size,replace = FALSE)
      
      Y_n <- Y[,-which(colnames(Y) %in% kickoutspe)]
      
      newposition_i<- which(colnames(Y_n)==species_i_name)
      newposition_j<- which(colnames(Y_n)==species_j_name)
      
      
      interinf_ij_n <-interMatrix_ij(Y_n,X,newposition_i,newposition_j,type,count,factor02,speciesprecision,envprecision,threshold01)
      
      #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij_n )
      interinf_ij_permut[n,1]<-interinf_ij_n 
      interinf_ij_permut[n,-1]<- kickoutspe
      
    }
    #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij)
    interinf_ij_permut[times+1,1]<-interinf_ij
  }else if (permut_type=="samples"){
    samplesnames <- rownames(Y)
    
    for (n in 1:times){
      kickoutsample<- sample(samplesnames,size,replace = FALSE)
      
      Y_n <- Y[-which(rownames(Y) %in% kickoutsample),]
      X_n<- X[-which(rownames(Y) %in% kickoutsample),]
      interinf_ij_n <-interMatrix_ij(Y_n,X_n,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
      #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij_n ) 
      interinf_ij_permut[n,1]<-interinf_ij_n 
      interinf_ij_permut[n,-1]<- kickoutsample
    }
    #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij)
    interinf_ij_permut[times+1,1]<-interinf_ij
  }else{
    paranames <- colnames(X)
    
    for (n in 1:times){
      kickoutpara<- sample(paranames,size,replace = FALSE)
      
      
      X_n<- X[,-which(rownames(X) %in% kickoutpara)]
      interinf_ij_n <-interMatrix_ij(Y,X_n,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
      #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij_n ) 
      interinf_ij_permut[n,1]<-interinf_ij_n 
      interinf_ij_permut[n,-1]<- kickoutsample
    }
    #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij)
    interinf_ij_permut[times+1,1]<-interinf_ij
    
  }
  
  return(interinf_ij_permut)
  
}



######### add the random error to the original data ####
interMatrix_ij_env_spe_robust <- function(Y,X,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01){
  ##### this function is used to calculate the robust test on species abundance and environmental parameter,  add the random error to the original data
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  #### times, how many times of the random test
  #### threshold, error bar level
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belongs to the types of "posi", "nega" or "no"
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  ### threshold means 5% or 10% or etc. 
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  N_para<-nsample*nEnv
  N_spe<-nsample*nSpe
  coef <- 1+threshold
  
  interinf_ij <-  interMatrix_ij(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
  
  testpermut <- c(1:times)
  
  
  interinf_ij_robust <-  data.frame(matrix(vector(),  times+1,1,dimnames=list( c(as.vector(testpermut),"original"),"values")),stringsAsFactors=F)
  
  for (n in 1:times){
    randerror_para<- matrix((1+runif(N_para,-1,1)*threshold),nsample,nEnv)
    X_n <- X*randerror_para
    
    randerror_spe<- matrix((1+runif(N_spe,-1,1)*threshold),nsample,nSpe)
    Y_n <- Y*randerror_spe
    
    
    interinf_ij_n <-interMatrix_ij(Y_n,X_n,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
    
    interinf_ij_robust[n,1]<-interinf_ij_n 
    
  }
  
  interinf_ij_robust[times+1,1]<-interinf_ij
  
  return(interinf_ij_robust)
}


interMatrix_ij_env_robust <- function(Y,X,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01){
  
  ##### this function is used to calculate the robust test on environmental parameter,  add the random error to the original data
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  #### times, how many times of the random test
  #### threshold, error bar level
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belongs to the types of "posi", "nega" or "no"
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  
  ### threshold means 5% or 10% or etc. 
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  N<-nsample*nEnv
  coef <- 1+threshold
  
  interinf_ij <-  interMatrix_ij(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
  
  testpermut <- c(1:times)
  
  
  interinf_ij_robust <-  data.frame(matrix(vector(),  times+1,1,dimnames=list( c(as.vector(testpermut),"original"),"values")),stringsAsFactors=F)
  
  for (n in 1:times){
    randerror<- matrix((1+runif(N,-1,1)*threshold),nsample,nEnv)
    X_n <- X*randerror
    
    interinf_ij_n <-interMatrix_ij(Y,X_n,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
    
    interinf_ij_robust[n,1]<-interinf_ij_n 
    
  }
  
  interinf_ij_robust[times+1,1]<-interinf_ij
  
  return(interinf_ij_robust)
}

interMatrix_ij_spe_robust <- function(Y,X,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01){
  
  ### threshold means 5% or 10% or etc. 
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  N<-nsample*nSpe
  coef <- 1+threshold
  
  interinf_ij <-  interMatrix_ij(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
  
  testpermut <- c(1:times)
  
  
  interinf_ij_robust <-  data.frame(matrix(vector(),  times+1,1,dimnames=list( c(as.vector(testpermut),"original"),"values")),stringsAsFactors=F)
  
  for (n in 1:times){
    randerror<- matrix((1+runif(N,-1,1)*threshold),nsample,nSpe)
    Y_n <- Y*randerror
    
    interinf_ij_n <-interMatrix_ij(Y_n,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
    
    interinf_ij_robust[n,1]<-interinf_ij_n 
    
  }
  
  interinf_ij_robust[times+1,1]<-interinf_ij
  
  return(interinf_ij_robust)
}



interMatrix_ij_permut02 <- function(Y,X,species_i,species_j,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01){
  ##### randomly remove one species  or parameter , save the 2 dimension interaction table into the array
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belong to the types of "posi", "nega" or "no"
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  ######  remove one  species  or one parameter #####
  #####  permut_type  is species  or para  #####
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  species_i_name<-colnames(Y)[species_i]
  species_j_name<-colnames(Y)[species_j]
  interinf_ij <-  interMatrix_ij(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
  
  if(permut_type=="species"){
    testspecies <- colnames(Y)[c(-species_i,-species_j)]
    interinf_ij_permut <-  data.frame(matrix(vector(),  length(testspecies)+1,1,dimnames=list( c(testspecies,"original"),c("values"))),stringsAsFactors=F)
    for (n in 1:length(testspecies)){
      Y_n <- Y[,-(which(colnames(Y)==testspecies[n]))]
      newposition_i<- which(colnames(Y_n)==species_i_name)
      newposition_j<- which(colnames(Y_n)==species_j_name)
      interinf_ij_n <-interMatrix_ij(Y_n,X,newposition_i,newposition_j,type,count,factor02,speciesprecision,envprecision,threshold01)
      interinf_ij_permut[n,1]<-interinf_ij_n 
    }
    #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij)
    interinf_ij_permut[length(testspecies)+1,1]<-interinf_ij
  }else{
    testpara <- colnames(X)
    interinf_ij_permut <-  data.frame(matrix(vector(),  length(testpara)+1,1,dimnames=list( c(testpara,"original"),c("values"))),stringsAsFactors=F)
    for (n in 1:length(testpara)){
      X_n <- X[,-(which(colnames(X)==testpara[n]))]
      interinf_ij_n <-interMatrix_ij(Y,X_n,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
      interinf_ij_permut[n,1]<-interinf_ij_n 
    }
    #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij)
    interinf_ij_permut[length(testpara)+1,1]<-interinf_ij
  }
  return(interinf_ij_permut)
}


####################################################################
######                    part G：试验数据奇点                ######
####################################################################
testsingula <- function(data) {
  #### 检验数据的奇异性检验数据是否有奇异性，如果有，去掉引起奇异性的列
  N1 <- ncol(data)     # 获取数据的列数，即变量的数量
  N2 <- qr(data)$rank  # 计算数据的秩（rank），即矩阵的线性独立列数
  names <- colnames(data)  # 获取数据的列名
  N3 <- N1 - N2        # 计算奇异性的数量，奇异性是列数减去秩（rank）
  
  if (N1 == N2) {      # 如果列数等于秩，说明矩阵没有奇异性
    results <- data    # 没有奇异性，返回原始数据
  } else {             # 如果列数不等于秩，说明存在奇异性
    temp <- data       # 创建一个数据副本用于处理
    
    while (qr(temp)$rank < ncol(temp)) {            # 检查数据的秩是否小于列数
      remove <- sample(names, N3, replace = FALSE)  # 随机选择要删除的列名（删除数量为 N3，即奇异性列的数量）
      temp <- data[, -which(names == remove)]       # 删除对应列的列名，使用 -which 来排除这些列
    }
    results <- temp   # 返回处理后的数据
  }
  return(results)     # 返回结果数据
}




####################################################################
######                part I：make some plots                 ######
####################################################################
# reshape the data for ggplot
reshapedata<- function(interactiondata,Sample_Info){
  #### 该函数的目的是对输入的交互作用数据进行重新整理，并合并样本信息数据，
  #### 它还会添加环境参数的数值，并将数据转换为适合绘图的格式
  x <- as.data.frame(t(interactiondata))   # interactiondata：包含物种之间交互作用的矩阵
  needed_colnames <- colnames(x)
  # 确保列名称清楚地显示交互值
  colnames(x) <- paste(colnames(x), "_IAV", sep = "")
  x$plot <- rownames(x)
  # merge：将 x 与 Sample_Info 中的 Plot 和环境参数列进行合并，确保每个 plot 对应其环境参数
  # Sample_Info：包含与样本相关的环境信息
  x <- merge(x, Sample_Info[, c("Plot", needed_colnames)], by.x = "plot", by.y = "Plot", sort = FALSE)
  # melt：将数据框 x 转换为长格式，每个样本、环境参数和交互作用值在一行，
  # 新生成的列名修改为 value_IAV 和 envpara_IAV，分别表示交互作用值和环境参数。
  xx <- melt(x, id = c("plot", needed_colnames))
  colnames(xx)[colnames(xx) == "value"] <- "value_IAV"
  colnames(xx)[colnames(xx) == "variable"] <- "envpara_IAV"
  # melt：再次将数据转换为长格式，形成每个交互作用和环境参数的一对多关系
  xxx <- melt(xx, id = c("plot", "envpara_IAV", "value_IAV"))
  # 通过去掉 _IAV 后缀并进行比较，检查环境参数名称是否匹配。如果匹配，则保留该行数据
  xxx$check <- gsub("_IAV", "", xxx$envpara_IAV)
  xxx$checktrue <- xxx$check == xxx$variable
  xxx <- droplevels(subset(xxx, checktrue == TRUE))
  return(xxx)   # 返回整理后的数据框 xxx，其包含每个样本、物种交互作用和环境参数的相关信息
}  



interinfPlot <- function(xxx,i,j,specieslist){
  #### 目的：用于绘制物种交互作用的图表，并保存为 PDF 文件，
  ####    它根据传入的物种对 i 和 j，绘制交互作用关系图，并为每个环境参数生成子图。
  # 这段代码添加了谁与谁交互的信息
  # remove any row with NA
  # xxx：由 reshapedata 函数返回的数据框，包含交互作用数据
  xxx <- xxx[which(apply(xxx, 1, function(x) !any(is.na(x)))), ]
  species_i<-specieslist[i]
  species_j<-specieslist[j]
  species_j_on_species_i <- xxx
  information <-paste("Interaction of", species_j, "on", species_i, sep=" ")
  species_j_on_species_i$comp <- rep(information,  nrow(species_j_on_species_i))
  
  # 下面的代码生成图，并写入一个单独的对象
  species_j_on_species_i_P <- ggplot(species_j_on_species_i, aes(value, value_IAV)) +
    theme_bw() +
    facet_wrap( ~ envpara_IAV, scales = "free", nrow = 1) +  # compare to: scales = "free_x"
    geom_hline(yintercept = 0, color = "green", size = 1) +
    geom_point() +
    ggtitle(species_j_on_species_i[1, "comp"])
  
  filenames <-paste(information,".pdf",sep="")
  setwd("./figures")
  # This code compiles all the figures into one plot.
  pdf(filenames, height = 15, width = 25)
  grid.arrange( species_j_on_species_i_P, ncol = 1)
  
  dev.off()
  setwd("..")
  
  return(species_j_on_species_i_P)
}

