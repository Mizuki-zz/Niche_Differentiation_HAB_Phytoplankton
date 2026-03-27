######                        Code Overview                         ######
##### This R file contains all functions required for interaction analysis
#####  part A --  Data filtering: FilterNA, FilterZero, Filterrepeat
#####  part B --  Calculation of derivatives and raw interaction values: rate_change, rate_change02, rate_change_low, interinf, interinf02, interinf03, interinf_low
#####  part C --  Calculation of interaction matrix: interMatrix, interMatrix02, interMatrix_norm, interMatrix_norm_summ, interMatrix_low, interMatrix_ij, interMatrix_ij_low
#####  part D --  : summ_inter_ij, summ_inter_ij_new, summ_interaction, anascatter
#####  part E --  Main environmental parameters: main_para, main_para_barplot
#####  part F --  Robust test: interMatrix_ij_permut, interMatrix_ij_env_spe_robust, interMatrix_ij_env_robust, interMatrix_ij_spe_robust, interMatrix_ij_permut02  
#####  part G --  Test singularity: testsingula
#####  part H --  Extract six examples from the interaction matrix: six_test_robust 
#####  part I --  Make plot of change of interaction value across the samples: reshapedata, interinfPlot




####################################################################
######          part A: Filter data with minimal missing info                ######
####################################################################
FilterNA <- function (data,type){
  #### Purpose: Clean missing values (NA) from data based on specified mode (remove columns or rows)
  #### data is data frame of species abundance or environmental parameter 
  #### type = 0, remove columns or rows containing missing values, prioritize removing the column or row with the most missing values
  #### type = 1, remove rows containing missing values
  #### type = other value, remove columns containing missing values
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
  #### Clean zero values from data
  #### data is data frame of species abundance or environmental parameter 
  #### type = 0, remove columns or rows containing zeros, prioritize removing the column or row with the most zeros
  #### type = 1, remove rows containing zeros
  #### type = other value, remove columns containing zeros
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
  N <- floor(threshold*nrow(data))  #### the threshold for the number of repeat values
  reptimes <- vector()
  for (j in 1: ncol(data)){
    reptimes <- c(reptimes, max(table(data[,j])))
  }
  data <- data[,which(reptimes< N)]
  return(data)
}


#################################################################
######           part B: Calculate first derivatives and interactions            ######
#################################################################  
#### This part contains all functions for calculating species abundance change rates and raw four‑dimensional interaction values.
#### Rate change calculations require species abundance and environmental parameter data. Interactions also require the results of rate change as input.
rate_change <- function(Y, X, i){
  #### This function calculates the first derivative of species i abundance with respect to all environmental parameters across all samples.
  #### Based on the data structure, the function automatically selects different splitting methods. Returns a two‑dimensional
  #### data frame: rows are samples, columns are environmental parameters.
  ####   Y is the abundance, data frame 
  ####   X is the environmental parameters, data frame
  ####   i is the index of species i
  nsample <- nrow(Y)
  nSpe <-ncol(Y)
  nEnv <- ncol(X)
  derEnv <- matrix(, nrow = nsample-1, ncol = nEnv)
  #   Main loop: for each sample n (as reference point), compare with other samples
  for (n in 1:(nsample -1)){
    # change of abundance
    # deltaSpe is a matrix showing the differences of each species between each other sample and sample n
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,]))))
    # deltaSpe_i is the change of the i‑th species we care about
    deltaSpe_i <- deltaSpe[,i]
    # change of environmental parameters
    # deltaEnv represents differences in environmental parameters between samples
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,]),1,function(x) x-as.matrix(X[n,]))))
    # Case 1: number of samples >= (number of species + number of environmental variables) (recommended to use linear regression)
    if (nsample>=(nSpe+nEnv)){
      # Regression model: regress the change of target species (Spe_i) on all environmental variables + changes of other species (multiple regression)
      regdata <- data.frame(cbind( deltaEnv,deltaSpe[,-i]))
      fmal<- as.formula (paste("deltaSpe_i ~ -1 + ",paste(colnames(regdata),collapse="+")))
      fit <- lm (fmal,  regdata)    
      derEnv[n,]<- fit$coefficients[1:nEnv]
    } else if (nsample>=nEnv){
      # Case 2: number of samples >= number of environmental variables, but insufficient to regress both environmental variables and other species changes
      # Use only environmental variables to fit (robust regression)
      regdata <- data.frame( deltaEnv)
      fmal<- as.formula (paste("deltaSpe_i ~ -1 + ",paste(colnames(regdata),collapse="+")))
      fit <- rlm (fmal,  regdata)  
      derEnv[n,]<- fit$coefficients[1:nEnv]
    }else{
      # Case 3: too few samples for effective regression, use approximation (two‑point difference)
      temp<-apply(deltaEnv,2,function(x) deltaSpe_i/x)
      temp[which(!is.finite(temp))] <- NA
      derEnv[n,] <- apply(temp,2,function(x) mean(x,na.rm=TRUE))
    } 
  }
  derEnv <- data.frame(derEnv)
  rownames(derEnv) <-rownames(Y)[1:(nsample-1)]
  colnames(derEnv) <-colnames(X)
  return(derEnv)   # each row represents an estimate from one sample difference, each column is a "response coefficient" or "derivative" for an environmental variable
}


rate_change02 <- function(Y, X, i){ 
  #### This function calculates the first derivative of species i abundance with respect to all environmental parameters across all samples n.
  #### Based on the data structure, it only considers the influence of environmental parameters. Returns a two‑dimensional
  #### data frame: rows are samples, columns are environmental parameters.
  #### It estimates the response of species abundance to environmental parameters via **robust linear regression**,
  #### and the returned result is the derivative for each sample for each environmental parameter.
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   i is the species i
  nsample <- nrow(Y)
  nSpe <-ncol(Y)
  nEnv <- ncol(X)
  # derEnv: matrix to store derivative estimates for each sample, size (nsample-1, nEnv)
  derEnv <- matrix(, nrow = nsample-1, ncol = nEnv)
  # Calculate changes in species abundance (Spe) and environmental parameters (Env) between each sample n and other samples
  for (n in 1:(nsample -1)){
    # change of abundance
    # deltaSpe_i is the change in abundance of the i‑th species
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,]))))
    deltaSpe_i <- deltaSpe[,i]
    # change of environmental parameters
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,]),1,function(x) x-as.matrix(X[n,]))))
    
    
    regdata <- data.frame( deltaEnv)
    
    # use robust linear regression
    # lm is not robust at all; rlm needs testing for data singularity
    # linear regression on deltaSpe ~ deltaEnv
    fmal<- as.formula (paste("deltaSpe_i ~ -1 + ",paste(colnames(regdata),collapse="+")))
    fit <- rlm (fmal, regdata,  maxit = 500,acc=1e-3) 
    # 200 is the number of iterations to achieve convergence; can be changed based on user testing and experience
    # This line extracts the coefficients for environmental parameters from the robust regression results (derivatives) and stores them in derEnv
    derEnv[n,]<- fit$coefficients[1:nEnv]
  }
  derEnv <- data.frame(derEnv)
  rownames(derEnv) <-rownames(Y)[1:(nsample-1)]
  colnames(derEnv) <-colnames(X)
  return(derEnv) 
}


rate_change_low <- function(Y, X, i){
  #### This function calculates the first derivative of species i abundance with respect to all environmental parameters across all samples.
  #### Based on the data structure, it uses a low‑precision method. Returns a two‑dimensional data frame.
  #### data frame: row is the samples, column is the environmental parameter
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   i is the species i
  ####   this function will calculate the first derivative of species i abundance with respect to all the environmental parameters in the sample n
  nsample <- nrow(Y)
  nSpe <-ncol(Y)
  nEnv <- ncol(X)
  derEnv <- matrix(, nrow = nsample-1, ncol = nEnv)
  for (n in 1:(nsample -1)){
    # change of abundance
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,]))))
    deltaSpe_i <- deltaSpe[,i]
    # change of environmental parameters
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,]),1,function(x) x-as.matrix(X[n,]))))
    # low precision, using two‑point difference to calculate the derivative, then take the median 
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
  #### This function calculates the interaction influence on species i across all samples.
  #### Based on the data structure, the function automatically selects different specific methods (robust regression, ordinary regression, low‑precision method).
  #### Returns a three‑dimensional array.
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   DYY is the derivative matrix with respect to environmental parameters X, representing the rate of change of each species abundance relative to environmental parameters for each sample
  ####   i is the species i
  #####   this function will calculate the interaction influence on species i 
  nsample <- nrow(DY)
  nEnv <- ncol(DY)
  nSpe <- ncol(Y)
  # When nsample >= (nSpe + nEnv), the model considers both species abundance and environmental parameters, using Taylor expansion
  if(nsample >=(nSpe+nEnv) ){
    # array: array[1] is the environmental parameter and the species, 
    #        array[2] is the environmental parameters, 
    #        array[3] is the samples
    # interaction array dimensions (nEnv + nSpe, nEnv, nsample-1), containing interactions for all species and environmental parameters
    interaction <- array (0, dim=c(nEnv+nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y),colnames(DY)),colnames(DY),rownames(Y[1:(nsample-1),])))
    # Spe_i <- Y[,i]    ### abundance of species i
    for (n in 1: (nsample-1)){
      # change of DY
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      # change of abundance
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      # change of environmental parameters
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      regdata <- data.frame(cbind(deltaSpe, deltaEnv))
      # linear regression on deltaDY ~ Spe + Env
      # Use robust regression to fit the model DY ~ Spe + Env and store coefficients in interaction, performing interaction calculation
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
    # This model only considers the influence of species abundance, using Taylor expansion
    # array: array[1] is the species, 
    #        array[2] is the environmental parameters, 
    #        array[3] is the samples
    interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
    for (n in 1: (nsample-1)){
      # change of DY
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      # change of abundance
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      # change of environmental parameters
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      regdata <- data.frame(cbind(deltaSpe))
      # linear regression on deltaDY ~ deltaSpe
      # Use ordinary linear regression (lm) to estimate the relationship between species abundance change and environmental parameter change
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
    # the model uses the low precision method
    # array: array[1] is the species, 
    #        array[2] is the environmental parameters, 
    #        array[3] is the samples
    interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
    for (n in 1: (nsample-1)){
      
      # change of DY
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      # change of abundance
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      # change of environmental parameters
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      temp<-apply(deltaSpe[1:(nsample-1),],c(1,2),function(x) 1/(x*nsample))
      temp[which(!is.finite(temp))] <- 0
      interaction [,,n] <-t(temp) %*% deltaDY    # this is equivalent to calculating the mean values
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
  #### Purpose: Calculate interactions between species and environment (using linear regression to analyze the influence of environmental changes on species abundance changes), and return a multi‑dimensional array representing these interactions. 
  #### Input: Y: species abundance data (a matrix, rows are samples, columns are species). 
  ####        X: environmental factor data (a matrix, rows are samples, columns are environmental factors).
  ####        DY: species abundance change data (changes of each species across samples). 
  ####        i: current species (column index). 
  
  # Initialize variables
  nsample <- nrow(DY)   # number of samples
  nEnv <- ncol(DY)      # number of environmental factors
  nSpe <- ncol(Y)       # number of species
  
  # Create interaction array:
  # interaction is a three‑dimensional array storing the interaction value for each species, each environmental factor, and each sample (minus one because changes are between samples).
  interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
  
  # Loop to calculate interactions for each species and environment
  for (n in 1: (nsample-1)){
    
    # Calculate changes in DY (changes in environmental factors between current sample and other samples) and Spe (changes in species abundance between current sample and other samples):
    deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
    
    #deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
    
    regdata <- data.frame(cbind(deltaSpe))
    for(k in 1:ncol(deltaDY)){
      deltaDY_k <- deltaDY[,k]
      
      
      # Linear regression model
      # By regressing DY_k on Spe, analyze the impact of environmental factor changes on species abundance changes
      fmal<- as.formula (paste("deltaDY_k ~ -1 + ",paste(colnames(regdata),collapse="+")))
      # This model uses rlm() (robust regression) to fit the data; acc relaxes convergence criteria
      fit <- rlm (fmal,  regdata[1:(nsample-1),], maxit=500,acc=1e-3)
      # Store interaction coefficients:
      interaction [,k,n] <- fit$coefficients/Y[n,i]
    }
    
    # Handle zero abundance cases
    # If the current species has zero abundance in the current sample, set the interaction value to NA
    if(Y[n,i]==0){
      interaction [,,n]<- NA
    } 
    missspecies <- which(Y[n,]==0,useNames = TRUE)
    interaction [missspecies,,n] <-NA
  }
  # Return a multi‑dimensional array interaction with dimensions nSpe x nEnv x (nsample - 1), representing interaction coefficients for each species, environmental factor, and sample.
  return(interaction)
}



interinf03 <- function(Y, X, DY, i){
  #### Purpose: Calculate interactions between species abundance changes and environmental factor changes via linear regression.
  #### Output: A three‑dimensional array containing interaction coefficients for each species, environmental factor, and sample.
  #### Input:
  ####       Y: species abundance matrix (samples × species)
  ####       X: environmental factor matrix (samples × environmental factors)
  ####       DY: species abundance change matrix (samples × species changes)
  
  # Initialize variables:
  nsample <- nrow(DY)    # number of samples
  nEnv <- ncol(DY)       # number of environmental factors (should be ncol(DY) actually, but note variable name)
  nSpe <- ncol(Y)        # number of species (columns of Y)
  
  # Create a three‑dimensional array interaction to store interaction coefficients for each species, environmental factor, and sample:
  interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
  
  # Loop over each sample (excluding the last one) to compute changes between adjacent samples.
  for (n in 1: (nsample-1)){
    
    # Calculate changes in species and environmental factors:
    # DY: change in environmental factors between current sample and other samples
    deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
    # Spe: change in species abundance between current sample and other samples
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
    # Env: change in environmental factor values between current sample and other samples
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
    
    # Linear regression model: deltaSpe ~ deltaEnv
    # Create a data frame containing species abundance changes Spe (predictors)
    regdata <- data.frame(cbind(deltaSpe))
    # fmal: build regression formula, using environmental factor changes DY as the response variable and species abundance changes Spe as predictors
    fmal<- as.formula (paste("deltaDY ~ -1 + ",paste(colnames(regdata),collapse="+")))
    # lm(fmal, ...): use ordinary least squares (OLS) for linear regression, fitting DY ~ Spe model
    fit <- lm (fmal,  regdata[1:(nsample-1),])
    # Store regression coefficients
    interaction [,,n] <- fit$coefficients
    # Handle zero abundance cases
    if(Y[n,i]==0){
      interaction [,,n]<- NA
    } 
    # Handle missing species
    missspecies <- which(Y[n,]==0,useNames = TRUE)
    interaction [missspecies,,n] <-NA
  }
  # Finally return a three‑dimensional array interaction containing interaction coefficients for each species, environmental factor, and sample.
  return(interaction)
}



interinf_low <- function(Y, X, DY, i){
  #### Purpose: Calculate interactions between species abundance changes and environmental factor changes.
  ####           Similar to interinf03 and interinf02, but uses a low‑precision method to estimate interactions.
  #### Input: 
  ####      Y: species abundance matrix, dimensions samples × species
  ####      X: environmental factor matrix, dimensions samples × environmental factors
  ####      DY: species abundance change matrix, dimensions samples × environmental factors, representing changes in species abundance for each sample
  ####      i: current species index (column index). 
  # Initialize variables
  nsample <- nrow(DY)   # number of samples
  nEnv <- ncol(DY)      # number of environmental factors
  nSpe <- ncol(Y)       # number of species
  # Create a three‑dimensional array interaction to store interaction coefficients for each species, environmental factor, and sample:
  interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
  # Loop over each sample
  for (n in 1: (nsample-1)){
    # Calculate changes in environmental factors (impact of abundance on environmental factors):
    deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
    # Calculate changes in species abundance:
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
    # Calculate changes in environmental factor values, primarily changes in environmental factors themselves across samples
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
    
    # Low‑precision method to calculate interactions:
    temp<-apply(deltaSpe[1:(nsample-1),],c(1,2),function(x) 1/(x*nsample))
    temp[which(!is.finite(temp))] <- 0
    interaction [,,n] <-t(temp) %*% deltaDY    ###### this is equivalent to calculating the mean values
    
    if(Y[n,i]==0){
      interaction [,,n]<- NA
    } 
    
    missspecies <- which(Y[n,]==0,useNames = TRUE)
    interaction [missspecies,,n] <-NA  
  }
  return(interaction)
}


#################################################################
######            part C: Calculate global interaction matrix                   ####
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



  #### Purpose: Calculate the global interaction matrix between species.
  ####          Determine the strength of interactions between species based on different statistical methods (e.g., mean, median, or weighted mean).
  #### Parameters: 
  ####      Y: species abundance matrix, dimensions nsample × nSpe, representing the abundance of each species in each sample.
  ####      X: environmental parameter matrix, dimensions nsample × nEnv, representing the values of different environmental factors in each sample.
  #### type: method for calculating interaction values:
  ####      type = 0: use mean
  ####      type = 1: use median
  ####      type = 2: use weighted mean
  #### count, integer: used to calculate histogram distribution in count segments based on minimum and maximum values.
  #### factor02, a factor used to compare the negative and positive parts.
  #### threshold01: parameter used in function summ_inter_ij_new.
  #### speciesprecision, or envprecision == "high", use rate_change02 and interinf02.
  #### speciesprecision, or envprecision == "low", use rate_change_low and interinf_low.
  
  interMatrix02 <- function(Y, X, type, count, factor02, speciesprecision, envprecision, threshold01, returnPerEnv = TRUE) {
    nsample <- nrow(Y)
    nEnv <- ncol(X)
    nSpe <- ncol(Y)
    results <- as.data.frame(matrix(vector(), nSpe, nSpe, dimnames = list(colnames(Y), colnames(Y))), stringsAsFactors = F)
    
    # New: store list of matrices for each environmental factor
    if (returnPerEnv) {
      perEnvResults <- vector("list", nEnv)
      names(perEnvResults) <- colnames(X)
      # Initialize matrix for each environmental factor
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
      
      # New: initialize matrices for each environmental factor for current species i if not already done
      if (returnPerEnv) {
        for (env in colnames(X)) {
          if (!env %in% names(perEnvResults)) {
            perEnvResults[[env]] <- matrix(NA, nrow = nSpe, ncol = nSpe, dimnames = list(colnames(Y), colnames(Y)))
          }
        }
      }
      
      for (species_j in 1:nSpe) {
        interinf_ij <- interaction[species_j, , ]
        
        # New: save value for current species pair (i,j) for each environmental factor
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
    
    # Return value based on returnPerEnv
    if (returnPerEnv) {
      return(list(globalMatrix = results, perEnvMatrices = perEnvResults))
    } else {
      return(results)
    }
  }


interMatrix_norm <- function(Y,X,speciesprecision,envprecision){
  #### Purpose: Calculate the strength of interactions between species based on species abundance (Y) and environmental parameters (X) using different precision settings.
  #### Parameters:
  ####    Y: species abundance matrix (usually a species‑sample matrix), rows are samples, columns are species
  ####    X: environmental parameter matrix, rows are samples, columns are environmental parameters
  ####    speciesprecision: species precision
  ####    envprecision: environmental precision
  #### speciesprecision, or envprecision == "high", use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low", use rate_change_low and interinf_low
  
  nsample <- nrow(Y)   # get number of samples
  nEnv <- ncol(X)      # get number of environmental parameters
  nSpe <- ncol(Y)      # get number of species
  # Create an empty result matrix nSpe x nSpe with row and column names as species names
  results <-  as.data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  # Outer loop iterates over all species (species_i), calculating interactions between species i and others
  for (species_i in 1:nSpe){
    if (envprecision=="high"){
      derEnv <- rate_change02 (Y, X, species_i)   # calculate environmental change rates
    }else {
      derEnv <- rate_change_low (Y, X, species_i) # calculate environmental change rates
    }
    # Remove columns from derEnv that contain NA. Use apply to check each column for NA and keep those without NA.
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    
    if (speciesprecision=="high"){
      interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i)   # calculate species interactions
    }else{
      interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i) # calculate species interactions
      
    }
    # Calculate interaction matrix and normalize:
    # Inner loop over all species species_j, extract interaction values for species i and j from interaction
    for (species_j in 1:nSpe){   
      interinf_ij <- interaction[species_j,,]
      temp<-interinf_ij
      temp[is.na(temp)] <- 0                       # replace NA with 0
      results[species_i,species_j] <- norm(temp)   # normalize the interaction and store
    }
  }    
  return(results) 
}


interMatrix_norm_summ <- function(Y,X,count,speciesprecision,envprecision){
  #### Purpose: Similar to interMatrix_norm, but with modifications when calculating the global interaction matrix.
  ####           It excludes extreme values by defining a histogram distribution and calculates the interaction matrix between species.
  #### Parameters:
  ####     Y: species abundance matrix. Rows are samples, columns are species.
  ####     X: environmental parameter matrix. Rows are samples, columns are environmental parameters.
  ####     count: number of bins for histogram. This parameter divides interaction values into intervals to compute their distribution.
  #### speciesprecision, or envprecision == "high", use rate_change02 and interinf02.
  #### speciesprecision, or envprecision == "low", use rate_change_low and interinf_low.
  
  
  nsample <- nrow(Y)   # number of samples (rows of Y)
  nEnv <- ncol(X)      # number of environmental parameters (columns of X)
  nSpe <- ncol(Y)      # number of species (columns of Y)
  # Create an nSpe x nSpe empty matrix to store interaction values between species, row and column names are species names from Y
  results <-  as.data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  # Outer loop over all species species_i, calculating interactions between species i and others
  for (species_i in 1:nSpe){
    if (envprecision=="high"){
      derEnv <- rate_change02 (Y, X, species_i)   # calculate environmental change rates
    }else {
      derEnv <- rate_change_low (Y, X, species_i) # calculate environmental change rates
    }
    # Remove columns from derEnv containing NA. Use apply to check each column for NA and keep only those without NA.
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    
    if (speciesprecision=="high"){
      interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i)   # calculate species interactions
    }else{
      interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i) # calculate species interactions
    }
    # Calculate interaction matrix and handle extreme values
    # Inner loop over all species species_j, extract interaction data for species i and j from interaction
    for (species_j in 1:nSpe){   
      interinf_ij <- interaction[species_j,,]
      temp<-c()   # initialize an empty vector to store computed results
      for (para in 1:nrow(interinf_ij)){   # iterate over each row of interaction data
        a <-interinf_ij[para,]             # extract interaction data for species j and all samples
        # Compute histogram distribution of interaction data
        # cut divides interaction data into count intervals, counts frequencies, generating a histogram
        counttable<-table(cut(interinf_ij[para,],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        # which.max finds the position of the maximum frequency in the histogram, corresponding to the peak of the data
        peakposition <- which.max(counttable)
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        # Extract data within the peak interval (peakvalue) and normalize using L2 norm (norm(peakvalue, type="2"))
        # Store the normalized result in the temp vector
        peakvalue <- a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])]
        temp<- c(temp,norm(peakvalue,type="2"))
      }
      # Finally, apply L2 normalization to the vector of normalized peak values and store in results[species_i, species_j]
      results[species_i,species_j]<- norm(temp,type="2")
    }    
  }    
  return(results)   # results contains the normalized interaction values between species
}


interMatrix_low <- function(Y,X,type,count,factor02,threshold01){
  #### Purpose: Based on environmental changes and species interactions, and using different methods to process interaction values.
  #### Parameters:
  ####     Y: species abundance matrix (samples x species). Rows are samples, columns are species.
  ####     X: environmental parameter matrix (samples x environmental parameters). Rows are samples, columns are environmental parameters.
  ####     type: parameter determining the method for calculating interactions. Its value determines which method to use to summarize interactions (e.g., mean, median, weighted sum, etc.).
  ####     factor02: a constant factor, used when type == 3 to compare medians of positive and negative interactions.
  ####     threshold01: threshold parameter used in summ_inter_ij_new, affecting interaction grouping.
  nsample <- nrow(Y)   # nsample: number of samples (rows of Y)
  nEnv <- ncol(X)      # nEnv: number of environmental parameters (columns of X)
  nSpe <- ncol(Y)      # nSpe: number of species (columns of Y)
  # Create an nSpe x nSpe empty matrix to store interaction values between species, row and column names are species names from Y
  results <-  as.data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  # Calculate interactions
  for (species_i in 1:nSpe){
    # Calculate environmental change rates
    derEnv <- rate_change_low (Y, X, species_i)
    # Remove columns containing NA
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    # Calculate species interactions
    interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i) 
    # Compute interaction based on type:
    for (species_j in 1:nSpe){   
      interinf_ij <- interaction[species_j,,]   # extract interaction data for species j
      if (type == 0){   # use mean of interaction data
        results[species_i,species_j] <- mean(interinf_ij,na.rm = TRUE)
      }
      if (type == 1){   # use median of interaction data
        results[species_i,species_j] <-median(interinf_ij,na.rm = TRUE)
      }
      if (type == 2){   # calculate median of weighted sum, where weights are the transpose of environmental parameters X
        paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
        results[species_i,species_j]<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
      }
      if (type == 3){   # special handling: determine interaction value using a more complex logic combining medians of positive and negative interactions and thresholds
        #factor02<-5
        # Use summ_inter_ij_new to process interaction data
        Summ_interinf_ij<- summ_inter_ij_new(interinf_ij,count,0,threshold01)
        number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
        # If there are more than 5 "no" type interactions, set to NA
        if (number_type_no>5){
          results[species_i,species_j] <- NA
        }else{
          # Calculate medians of positive and negative interactions
          Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
          Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
          # Count numbers of positive and negative interactions
          posipara<-length(which(Summ_interinf_ij$intertype=="+"))
          negapara<-length(which(Summ_interinf_ij$intertype=="-"))
          
          if(posipara==0){
            results[species_i,species_j] <-  Negamedian   # only negative interactions, result = negative median
          }else if (negapara==0){
            results[species_i,species_j] <- Posimedian    # only positive interactions, result = positive median
          }else if (Posimedian>(factor02*abs(Negamedian))){
            results[species_i,species_j]<- Posimedian     # if positive interaction exceeds threshold times negative, use positive median
          }else if ((factor02*Posimedian)<(abs(Negamedian))){
            results[species_i,species_j] <- Negamedian    # if negative interaction exceeds threshold times positive, use negative median
          }else {
            results[species_i,species_j] <- 0             # otherwise set to 0
          }           
        }         
      }
    }    
  }
  
  return(results)   # return the interaction matrix results, containing interaction values between species
}


interMatrix_ij <- function(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01){
  #### Purpose: Calculate a single interaction value between species (i.e., β_ij, representing the effect of species i on species j).
  ####          This is part of constructing an element of the overall species interaction matrix.
  #### Parameters:
  ####     Y : species abundance matrix (samples x species)
  ####     X : environmental factor matrix (samples x environmental variables)
  ####     type : method for calculating interaction (0=mean, 1=median, 2=weighted sum, 3=workflow method for statistical estimation)
  ####     count : integer: used to calculate histogram distribution in count segments based on minimum and maximum values
  ####     factor02 : factor for comparing positive and negative interactions
  ####     threshold01: parameter used in summ_inter_ij_new
  #### speciesprecision, or envprecision == "high", use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low", use rate_change_low and interinf_low
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  # Calculate environmental change rates
  if (envprecision=="high"){
    derEnv <- rate_change02 (Y, X, species_i)
  }else {
    derEnv <- rate_change_low (Y, X, species_i)
  }
  # Remove columns containing NA to ensure data integrity
  derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
  # Calculate species interaction matrix (for all species)
  if (speciesprecision=="high"){
    interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i) 
  }else{
    interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)    
  }
  # Extract interaction data for species j on species i, which is a matrix or vector representing influence values under different environments or times
  interinf_ij <- interaction[species_j,,]
  if (type == 0){   # use mean
    results <- mean(interinf_ij,na.rm = TRUE)
  }
  if (type == 1){   # use median
    results <-median(interinf_ij,na.rm = TRUE)
  } 
  if (type == 2){   # weighted sum (considering environmental parameters)
    paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
    results<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
  }
  if (type == 3){   # use workflow method for statistical estimation
    # Call summ_inter_ij_new to classify interaction values into: "+": positive interactions; "-": negative interactions; "no": non‑significant interactions
    Summ_interinf_ij<- summ_inter_ij_new(interinf_ij,count,0,threshold01)
    number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
    # Determine if too many "no" interactions
    if (number_type_no>5){
      results<- NA   # if number of "no" exceeds 5, consider data too weak and return NA
    }else{
      # Otherwise analyze positive and negative interactions:
      Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
      Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
      # Count numbers of positive and negative
      posipara<-length(which(Summ_interinf_ij$intertype=="+"))
      negapara<-length(which(Summ_interinf_ij$intertype=="-"))
      # Compare medians of positive and negative interactions, using factor02:
      # Purpose: Determine final interaction direction and magnitude based on strength of positive/negative interactions
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
  # Purpose: Calculate the interaction value of species i on species j, similar to interMatrix_ij,
  #          but using low‑precision environmental change rate and interaction inference methods, utilizing low‑precision functions (e.g., rate_change_low and interinf_low) to estimate interactions.
  # Parameters: same as above
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  # Calculate environmental change rates
  derEnv <- rate_change_low (Y, X, species_i)
  # Remove columns containing NA
  derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
  # Calculate interactions
  interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)
  # Extract interaction data for species j
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
  if (type==3){   # Call summ_inter_ij_new for further statistical inference on interaction values
    Summ_interinf_ij<- summ_inter_ij_new(interinf_ij,count,0,threshold01)
    number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
    if (number_type_no>5){
      results<- NA
    }else{
      Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
      Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
      
      posipara<-length(which(Summ_interinf_ij$intertype=="+"))
      negapara<-length(which(Summ_interinf_ij$intertype=="-"))
      # Make decision based on counts and magnitudes of positive/negative interactions
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
######         part D: Statistical summary of interaction results              ######
#################################################################
# After calculating the raw four‑dimensional interaction values, we need to perform some statistical summaries to determine the global interaction matrix.
# Based on the pattern, we can determine whether the global interaction is positive, negative, or has no clear pattern.

summ_inter_ij <- function(Y,X,species_i,species_j,count,type,speciesprecision,envprecision){
  #### Purpose: This function calculates the statistical summary α of the raw interaction values β^k_ij.
  #### Parameters:
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
  # Calculate environmental change rates (derEnv)
  if (envprecision=="high"){
    derEnv <- rate_change02 (Y, X, species_i)
  }else {
    derEnv <- rate_change_low (Y, X, species_i)
  }
  # Remove environmental change rates with missing values
  derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
  # Calculate species interactions
  if (speciesprecision=="high"){
    interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i) 
  }else{
    interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)
  }
  # Extract interaction data for species i on species j, a three‑dimensional array
  # Create column names for statistical results
  interinf_ij <- interaction[species_j,,]
  b<- as.character(c(1:count))
  b<- t(as.matrix(as.character(c(1:(count-1)))))
  b<- apply(b,1,function(x) paste("count",x,sep=""))
  # Create a list of column names for statistical results, including:
  # intertype: interaction type (positive, negative, none).
  # min, 1stQu, median, mean, 3rdQu, max: basic statistics of the data.
  # peakrange: peak range of interaction data.
  # peakvalue: median of the peak.
  # variance: variance of the data.
  # numpositive, numnegative: counts of positive and negative interactions.
  # count1, count2, ..., countN: histogram bins according to count.
  namelist <- c("intertype","min","1stQu","median","mean","3rdQu","max", "peakrange", "peakvalue", "variance",  "numpositive","numnegative",b)
  # Calculate statistical properties of interactions (type == 0, by sample)
  if (type==0){
    sta_inter <-  data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
    for (para in 1:nrow(interinf_ij)){
      a <-interinf_ij[para,]
      temp <- as.numeric(summary(a))[1:6]   # summary(a): calculate basic statistics of interaction data
      variance <- var(a,na.rm = TRUE)       # var(a): calculate variance of interaction data
      nposi <- sum(  a > 0,na.rm = TRUE )   # nposi: count of positive interactions
      nnega <- nsample-nposi-sum(is.na(a))  # nnega: count of negative interactions
      n_non_NA<-nsample- sum(is.na(a))
      if (nposi/nsample>0.8) {
        # If proportion of positive interactions exceeds 80%, classify as positive (+)
        intertype <- "+"    
      }else if (nnega/nsample>0.8){
        # If proportion of negative interactions exceeds 80%, classify as negative (-)
        intertype <-"-"
      }else {
        # Otherwise, classify as no clear interaction (no)
        intertype <-"no"
      }
      
      if (n_non_NA== 1 | n_non_NA== 0){
        sta_inter[para,1:6]
      }else{
        # Use cut to divide interaction values into count intervals
        counttable<-table(cut(interinf_ij[para,],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        peakposition <- which.max(counttable)   # peakposition indicates the location of the peak in interaction data
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        peakvalue <- median(a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])])
        sta_inter[para,c(-8,-1)]<-c(temp,peakvalue,variance,nposi,nnega,as.numeric(counttable))
        sta_inter[para,1]<-intertype
        sta_inter[para,8]<-peakrange
      }
    }
  }else {
    # Calculate statistical properties of interactions (type != 0, by environmental parameter)           
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
  #### This function calculates the statistical summary α of the raw interaction values β^k_ij.
  ####  interinf_ij is the interaction data frame of beta_ij
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  #### threshold01: threshold for judging interactions, determines whether interaction is positive, negative, or has no clear pattern.
  # Initialization and column name setup
  nsample <- ncol(interinf_ij)   # nsample: get number of samples, i.e., number of columns of interinf_ij
  # b: generate a vector of count bin names for describing segmentation intervals of interaction data
  b<- as.character(c(1:count))
  b<- t(as.matrix(as.character(c(1:(count-1)))))
  b<- apply(b,1,function(x) paste("count",x,sep=""))
  #namelist: set column names of statistical results, including basic statistics (min, quartiles, median, etc.),
  # peak information (peakrange and peakvalue), interaction type (intertype), and bin statistics
  namelist <- c("intertype","min","1stQu","median","mean","3rdQu","max", "peakrange", "peakvalue", "variance",  "numpositive","numnegative",b)
  # Statistical summary calculation: by sample (type == 0)
  if (type==0){
    sta_inter <-  data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
    for (para in 1:nrow(interinf_ij)){
      a <-interinf_ij[para,]
      temp <- as.numeric(summary(a))[1:6]   # summary(a) calculates basic statistics of data
      variance <- var(a,na.rm = TRUE)       # var(a): calculate variance of data
      nposi <- sum(  a > 0,na.rm = TRUE )   # nposi: count of positive values, representing positive interactions
      nnega <- nsample-nposi-sum(is.na(a))  # nnega: count of negative values, representing negative interactions
      n_non_NA<-nsample- sum(is.na(a))      # n_non_NA: number of non‑missing samples
      if (nposi/nsample>threshold01) {
        # If proportion of positive interactions exceeds threshold01, mark as positive (+)
        intertype <- "+"
      }else if (nnega/nsample>threshold01){
        # If proportion of negative interactions exceeds threshold01, mark as negative (-)
        intertype <-"-"
      }else {
        # Otherwise, mark as no clear pattern (no)
        intertype <-"no"
      }
      # Data binning and peak calculation
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
    # Statistical summary calculation: by environmental parameter (type != 0)                
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
  return(sta_inter) # contains all computed statistical information, such as interaction type, peak range, variance, counts of positive/negative interactions, etc.
}


summ_interaction <-function(Y,X,count,type){
  # Purpose: Calculate interaction statistics between species
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  # b: create a character vector representing interaction bins
  b<- t(as.matrix(as.character(c(1:(count-1)))))
  b<- apply(b,1,function(x) paste("count",x,sep=""))
  # namelist: names of statistics to compute, such as min, quartiles, median, mean, variance, Pearson correlation, Spearman correlation, etc.
  namelist <- c("min","1stQu","median","mean","3rdQu","max","variance", "PearsonCorr", "PearsonPvalue","SpearmanCorr","SpearmanPvalue",  "numpositive","numnegative",b)
  # type == 0, statistics are calculated by species interaction
  if (type==0){
    # sta_inter: initialize a four‑dimensional array to store interaction statistics for each pair of species (species_i and species_j),
    # dimensions (nSpe, nSpe, nEnv, length(namelist)), corresponding to species pairs, environmental factors, and statistic names
    sta_inter <-  array(NA, dim=c(nSpe,nSpe, nEnv, length(namelist)),dimnames=list(colnames(Y),colnames(Y),colnames(X),namelist))
    # Calculate interactions for each pair of species
    for(species_i in 1:nSpe){
      for(species_j in 1:nSpe){
        # Call rate_change to calculate environmental changes in response to species
        derEnv <- rate_change (Y, X, species_i)
        derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
        # Calculate interaction matrix via interinf
        interaction <- interinf (Y, X[,colnames(derEnv)], derEnv, species_i)     
        # Extract interaction data for species_j on species_i from the interaction matrix
        interinf_ij <- interaction[species_j,,]
        # Calculate statistics:
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
    # type != 0: calculate statistics by sample
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
  #### Purpose: Calculate and return extreme value differences of species interactions and their relationship differences with environmental parameters.
  #### interaction_ij is a two‑dimensional matrix of interaction strengths between a pair of species, rownames are environmental parameters, colnames are samples ###
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  position <- which(rownames(interaction_ij)==para)   # find the row number for para (i.e., the position of this environmental parameter in the interaction matrix)
  window <-5   # window: size of the window, number of samples considered in each analysis (default = 5)
  N<- floor(nsample/window)   # N: number of samples divided by window size, giving the number of windows
  # Create name list and data frame
  b<- t(as.matrix(as.character(c(1:(nsample-1-window)))))
  b<- apply(b,1,function(x) paste("window",x,sep=""))
  namelist <- c("extreme",b)
  difference<-  data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
  paradiffer <- data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
  # Loop to calculate interaction differences:
  for (n_para in (1:nrow(interaction_ij))){
    paraorig <-X[colnames(interaction_ij),n_para]
    paraorder<- order(paraorig)
    counttable<-table(cut(paraorig,breaks=seq(min(paraorig,na.rm = TRUE),max(paraorig,na.rm = TRUE),length=count))) #### check data extreme properties 
    diff_temp <- vector(mode="numeric", length=0)
    paradiff_temp <- vector(mode="numeric", length=0)
    # If variance of frequencies >20, data is considered to have large fluctuations.
    # In this case, the program re‑bins paraorig so that data are concentrated in the interval with the highest frequency.
    if (var(counttable)>20){
      breakpoint <- seq(min(paraorig,na.rm = TRUE),max(paraorig,na.rm = TRUE),length=count)
      pp <- which(counttable == max(counttable,na.rm = TRUE))
      paraorig <- paraorig[which(paraorig >breakpoint[pp] & paraorig <= breakpoint[pp+1])]
      paraorder<- order(paraorig)
    }
    # Calculate extreme value differences of interactions and differences of environmental parameters:
    for (k in 1:(nsample-1-window)){
      # section_k: interaction data for the current window
      section_k <- interaction_ij[position,paraorder[k:(k+window-1)]]
      # diff_temp: difference between max and min interaction in each window
      diff_temp<-c(diff_temp,max(section_k,na.rm = TRUE)-min(section_k,na.rm = TRUE))
      # paradiff_temp: difference between max and min environmental parameter (paraorig) in each window
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
######                part E: Calculate main environmental parameters                ######
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
######                 part F: robust test                 ######
#################################################################
interMatrix_ij_permut <- function(Y,X,species_i,species_j,times,size,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01){
  
  
  ###### randomly remove size species for times, save the 2 dimension interaction table into the data frame #####
  #####  permut_type is "samples", "species", "para"  
  ##### 
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  ####  time, how many times in the random test
  #### size, control the number of movement of species or parameters
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



######### add random error to the original data ####
interMatrix_ij_env_spe_robust <- function(Y,X,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01){
  ##### this function is used to calculate the robust test on species abundance and environmental parameter, add random error to the original data
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
  
  ##### this function is used to calculate the robust test on environmental parameter, add random error to the original data
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
  ##### randomly remove one species or parameter, save the 2 dimension interaction table into the array
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
  
  ######  remove one species or one parameter #####
  #####  permut_type is species or para #####
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
######                    part G: Test data singularity                ######
####################################################################
testsingula <- function(data) {
  #### Test data singularity: check if data has singularities; if so, remove the columns causing singularity.
  N1 <- ncol(data)     # number of columns, i.e., number of variables
  N2 <- qr(data)$rank  # rank of the matrix (number of linearly independent columns)
  names <- colnames(data)  # column names of the data
  N3 <- N1 - N2        # number of singularities = number of columns minus rank
  
  if (N1 == N2) {      # if number of columns equals rank, matrix has no singularity
    results <- data    # no singularity, return original data
  } else {             # if columns ≠ rank, singularity exists
    temp <- data       # create a copy of the data for processing
    
    while (qr(temp)$rank < ncol(temp)) {            # check if rank is less than number of columns
      remove <- sample(names, N3, replace = FALSE)  # randomly select column names to remove (remove N3 columns, i.e., number of singular columns)
      temp <- data[, -which(names == remove)]       # remove corresponding columns using -which
    }
    results <- temp   # return processed data
  }
  return(results)     # return the result
}




####################################################################
######                part I: make some plots                 ######
####################################################################
# reshape the data for ggplot
reshapedata<- function(interactiondata,Sample_Info){
  #### The purpose of this function is to reshape the input interaction data and merge with sample information,
  #### also adding environmental parameter values and converting the data into a format suitable for plotting.
  x <- as.data.frame(t(interactiondata))   # interactiondata: matrix containing species interactions
  needed_colnames <- colnames(x)
  # Ensure column names clearly show interaction values
  colnames(x) <- paste(colnames(x), "_IAV", sep = "")
  x$plot <- rownames(x)
  # merge: merge x with Sample_Info columns Plot and environmental parameters, ensuring each plot corresponds to its environmental parameters
  # Sample_Info: contains environmental information related to samples
  x <- merge(x, Sample_Info[, c("Plot", needed_colnames)], by.x = "plot", by.y = "Plot", sort = FALSE)
  # melt: convert data frame x to long format, each sample, environmental parameter, and interaction value in one row,
  # new column names are value_IAV and envpara_IAV, representing interaction value and environmental parameter respectively.
  xx <- melt(x, id = c("plot", needed_colnames))
  colnames(xx)[colnames(xx) == "value"] <- "value_IAV"
  colnames(xx)[colnames(xx) == "variable"] <- "envpara_IAV"
  # melt: again convert to long format, creating a one‑to‑many relationship between each interaction and environmental parameter
  xxx <- melt(xx, id = c("plot", "envpara_IAV", "value_IAV"))
  # Remove "_IAV" suffix and compare to check if environmental parameter names match. If they match, keep the row.
  xxx$check <- gsub("_IAV", "", xxx$envpara_IAV)
  xxx$checktrue <- xxx$check == xxx$variable
  xxx <- droplevels(subset(xxx, checktrue == TRUE))
  return(xxx)   # return reshaped data frame xxx, containing information on each sample, species interaction, and environmental parameter
}  



interinfPlot <- function(xxx,i,j,specieslist){
  #### Purpose: Plot species interaction graphs and save as PDF files.
  ####          Based on the input species pair i and j, it draws interaction relationship plots, creating subplots for each environmental parameter.
  # This code adds information about which species interact with which.
  # remove any row with NA
  # xxx: data frame returned by reshapedata, containing interaction data
  xxx <- xxx[which(apply(xxx, 1, function(x) !any(is.na(x)))), ]
  species_i<-specieslist[i]
  species_j<-specieslist[j]
  species_j_on_species_i <- xxx
  information <-paste("Interaction of", species_j, "on", species_i, sep=" ")
  species_j_on_species_i$comp <- rep(information,  nrow(species_j_on_species_i))
  
  # The following code generates the plot and writes it into a separate object
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