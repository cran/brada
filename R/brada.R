# brada - Bayesian response adaptive trial design analysis
#     Copyright (C) 2022  Riko Kelter
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' function for simulating and analyzing a Bayesian response-adaptive trial design for binary endpoints
#'
#' @method 
#' @export
#' @examples
brada <- function(a0=1,b0=1,Nmax=40,batchsize=5,nInit,p_true,p0,p1,theta_T=0.90,theta_L=0.1,theta_U=1,nsim=100,seed=42,method="PP",refFunc="flat",nu=0,shape1=1,shape2=1,truncation=1,cores=2){
  # WE MUST LEAVE THESE TWO SINGLE-TRIAL SIMULATION FUNCTIONS INSIDE AND AT THE TOP OF THIS FUNCTION TO ENSURE THAT PARALLEL COMPUTING RECOGNIZES THEM; 
  # OTHERWISE THEY MUST BE INCLUDED IN A SEPARATE R PACKAGE WHICH CAN BE LOADED VIA THE .packages ARGUMENT OF THE foreach FUNCTION
  
  ##################################################################################################################################
  singleTrial_PP <- function(s,n,responseMatrix,nInit,Nmax,batchsize,a0,b0,p0,p1){
    n=nInit # first n patients
    x <- sum(responseMatrix[s,1:n]) # number of successes in first n patients
    stop = FALSE
    PPs <- NULL # vector for later line plot
    tempMatrix = matrix(-1,nrow=1,ncol=2+length(seq(nInit,Nmax,by=batchsize)))
    while(stop == FALSE){
      if(!(n == nInit)){
        sumNextBatch = sum(responseMatrix[s,(n-batchsize+1):n])
        x = x+sumNextBatch # add batchsize new observations if n != nInit (e.g. 10)
      }
      PP = 0 # initialize PP
      for(i in 0:(Nmax-n)){ # iterate over possible results for successes
        marg=extraDistr::dbbinom(i,size=Nmax-n,a0+x,b0+n-x) # P(Y=i|x)
        postProb=1-pbeta(p1,a0+x+i,b0+Nmax-x-i) # P(p>p1|X=x,Y=i) = evidence for H1
        indicator = 0
        if(postProb > theta_T){
          indicator = 1
        }
        PP = PP + (marg*indicator)
      }
      print(PP); print(n)
      PPs = c(PPs,PP) # append PP value for line plot later
      if(PP <= theta_L){ # check for futility
        # stop for futility (futility is coded as 1)
        tempMatrix[1,1:2] = c(1,n) # temporary matrix which stores the result
        #trials<-rbind(trials,data.frame(futility=1,n=n,PP=PP))
        stop = TRUE # stop trial
        #lines(seq(nInit,nInit+(length(PPs)-1)*batchsize,by=batchsize),PPs,ty="l",col=rgb(255, 0, 0, max = 255, alpha = 125)) # add line to plot
      }
      if(PP >= theta_U){ # check for efficacy
        # stop for efficacy (efficacy is coded as 2)
        tempMatrix[1,1:2] = c(2,n) # temporary matrix which stores the result
        stop = TRUE # stop trial
      }
      if(n==Nmax && stop == FALSE){
        stop = TRUE # stop trial
        # trials<-rbind(trials,data.frame(futility=0,n=n,PP=PP)) # store data
        tempMatrix[1,1:2] = c(0,n) # temporary matrix which stores the result
        #lines(seq(nInit,Nmax,by=batchsize),PPs,ty="l",col=rgb(100, 149, 237, max = 255, alpha = 125)) # add line to plot
      }
      if(n<Nmax && stop == FALSE){
        n=n+batchsize # increment n by batchsize
      }
    }
    for(j in 1:length(PPs)){
      tempMatrix[1,2+j] = PPs[j]
    }
    # Set all values inside tempMatrix to the last PP value where PP either reached the futility or efficacy threshold
    if(min(tempMatrix[1,]) == -1){
      tempMatrix[1,min(which(tempMatrix[1,] == -1)):ncol(tempMatrix)] = tempMatrix[1,min(which(tempMatrix[1,] == -1))-1]
    }
    tempMatrix
  }
  
  
  singleTrial_PPe <- function(s,n,responseMatrix,nInit,Nmax,batchsize,a0,b0,refFunc,shape1,shape2,truncation=1,p0,p1){
    n=nInit # first n patients
    x <- sum(responseMatrix[s,1:n]) # number of successes in first n patients
    stop = FALSE
    PPes <- NULL # vector for later line plot
    tempMatrix = matrix(-1,nrow=1,ncol=2+length(seq(nInit,Nmax,by=batchsize)))
    while(stop == FALSE){
      if(!(n == nInit)){
        sumNextBatch = sum(responseMatrix[s,(n-batchsize+1):n])
        x = x+sumNextBatch # add batchsize new observations if n != nInit (e.g. 10)
      }
      PPe = 0 # initialize PP
      for(i in 0:(Nmax-n)){ # iterate over possible results for successes
        marg=extraDistr::dbbinom(i,size=Nmax-n,a0+x,b0+n-x) # P(Y=i|x)
        # postDensDraws = rbeta(500000,a0+x+i,b0+Nmax-x-i)
        if(refFunc=="beta"){
          # FBET=fbet(posteriorDensityDraws = postDensDraws,
          #          interval=c(p1,1),
          #          FUN = dbeta,
          #          par = list(shape1=shape1,shape2=shape2),
          #          nu=nu) # evalue for alternative
          FBET_H1 = fbst::fbet(posterior=dbeta, # analytic posterior version; increases computational speed
                         par_posterior=list(shape1=a0+x+i,shape2=b0+Nmax-x-i),
                         interval=c(p1,1),
                         FUN = dbeta,
                         par = list(shape1=shape1,shape2=shape2),
                         nu=nu)
        }
        if(refFunc=="binaryStep"){
          FBET_H1 = fbst::fbet(posterior=dbeta, # analytic posterior version; increases computational speed
                               par_posterior=list(shape1=a0+x+i,shape2=b0+Nmax-x-i),
                               interval=c(p1,1),
                               FUN = binaryStep_truncated,
                               par = list(p0=p0,truncation=truncation),
                               nu=nu)
        }
        if(refFunc=="relu"){
          FBET_H1 = fbst::fbet(posterior=dbeta, # analytic posterior version; increases computational speed
                               par_posterior=list(shape1=a0+x+i,shape2=b0+Nmax-x-i),
                               interval=c(p1,1),
                               FUN = genShiftedReLU_truncated,
                               par = list(p0=p0,truncation=truncation),
                               nu=nu)
        }
        if(refFunc=="palu"){
          FBET_H1 = fbst::fbet(posterior=dbeta, # analytic posterior version; increases computational speed
                               par_posterior=list(shape1=a0+x+i,shape2=b0+Nmax-x-i),
                               interval=c(p1,1),
                               FUN = PaLU_truncated,
                               par = list(p0=p0,truncation=truncation),
                               nu=nu)
        }
        if(refFunc=="lolu"){
          FBET_H1 = fbst::fbet(posterior=dbeta, # analytic posterior version; increases computational speed
                               par_posterior=list(shape1=a0+x+i,shape2=b0+Nmax-x-i),
                               interval=c(p1,1),
                               FUN = LoLU_truncated,
                               par = list(p0=p0,truncation=truncation),
                               nu=nu)
        }
        if(refFunc=="flat"){
          #FBET=fbet(posteriorDensityDraws = postDensDraws,
          #          interval=c(p1,1),
          #          nu=nu) # evalue for alternative
          FBET_H1 = fbst::fbet(posterior=dbeta, # analytic posterior version; increases computational speed
                         par_posterior=list(shape1=a0+x+i,shape2=b0+Nmax-x-i),
                         interval=c(p1,1),
                         FUN = dbeta, # beta(1,1) ref function is equal to flat reference function
                         par = list(shape1=1,shape2=1),
                         nu=nu)
        }
        eValue = FBET_H1 # this is the evidence for H1:p>p1
        indicator = 0
        if(eValue > theta_T){
          indicator = 1
        }
        PPe = PPe + (marg*indicator)
      }
      PPes = c(PPes,PPe) # append PPe value for line plot later
      if(PPe <= theta_L){ # check for futility
        # stop for futility (futility is coded as 1)
        tempMatrix[1,1:2] = c(1,n) # temporary matrix which stores the result
        stop = TRUE # stop trial
      }
      if(PPe >= theta_U){ # check for efficacy
        # stop for efficacy (efficacy is coded as 2)
        tempMatrix[1,1:2] = c(2,n) # temporary matrix which stores the result
        stop = TRUE # stop trial
      }
      if(n==Nmax && stop == FALSE){
        stop = TRUE # stop trial
        tempMatrix[1,1:2] = c(0,n) # temporary matrix which stores the result
      }
      if(n<Nmax && stop == FALSE){
        n=n+batchsize # increment n by batchsize
      }
    }
    for(j in 1:length(PPes)){
      tempMatrix[1,2+j] = PPes[j]
    }
    # Set all values inside tempMatrix to the last PP value where PP either reached the futility or efficacy threshold
    if(min(tempMatrix[1,]) == -1){
      tempMatrix[1,min(which(tempMatrix[1,] == -1)):ncol(tempMatrix)] = tempMatrix[1,min(which(tempMatrix[1,] == -1))-1]
    }
    tempMatrix
  }
  ##################################################################################################################################
  
  # ANN reference functions
  binaryStep_truncated <- function(x,p0,truncation=1){
    ifelse(x<p0 | x>truncation, 0, 1/(1-p0))
  }
  
  genShiftedReLU_truncated <- function(x,p0,truncation=1){
    ifelse(x<p0 | x>truncation, 0, (x-p0)*2/((1-p0)^2))
  }
  
  PaLU_truncated <- function(x,p0,truncation=1){
    ifelse(x<p0 | x>truncation,0,((1/3-p0*(1/3*p0^2+1-p0))^(-1))*(x-p0)^2)
  }
  
  LoLU_truncated <- function(x,p0,truncation=1){
    ifelse(x<p0 | x>truncation,0,(1/(p0-(p0-2)*log(2-p0)-1))*log(1+x-p0))
  }
  # plot(seq(0,1,by=0.01),genShiftedReLU_truncated(x=seq(0,1,by=0.01),p0=0.2,truncation=0.9),ty="l")
  
  ##################################################################################################################################
  
  # Set up cluster
  
  # num_workers=parallel::detectCores()  # setup parallel backend to use many processors
  # CRAN limits number of cores to 2 due to performance reasons, so this is checked first
  # chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  # if (nzchar(chk) && (chk != "false")){  # then limit the workers
  #   num_workers <- 2L
  # } else {
  #   # use all cores
  #   num_workers <- parallel::detectCores()
  # }
  #   
  # chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  cl <- parallel::makeCluster(cores)
  #doParallel::registerDoParallel(cl)
  doSNOW::registerDoSNOW(cl)
  ##################################################################################################################################
  
  # Parameters
  # library(extraDistr)
  #a0=a0 # beta prior parameters
  #b0=b0
  #Nmax=Nmax # maximum number of patients
  p=p_true # true response probability
  #p0 = p0 # H0 boundary
  #p1 = p1 # H1 boundary
  #theta_T = theta_T # threshold for posterior prob
  #theta_U = theta_U # threshold for efficacy for PP
  #theta_L = theta_L # threshold for futility for PP
  set.seed(seed)
  #nsim=nsim
  ##################################################################################################################################
  
  # Generate data
  responseMatrix = generateData(p=p,Nmax=Nmax,nsim=nsim,seed=seed)
  
  ##################################################################################################################################
  
  # Plot
  # graphics.off()
  # df <- data.frame()
  # resultPlot <- ggplot(df) + geom_point() + xlim(nInit, Nmax) + ylim(0, 1) + 
  #   theme_bw() + 
  #   geom_hline(yintercept=theta_L) + 
  #   geom_hline(yintercept=theta_U) +
  #   geom_vline(xintercept = Nmax, linetype="dotted", 
  #               color = "black", size=1) +
  #   labs(x = "Patients included in the trial", ylab="PPe")
  # graphics.off()
  
  ##################################################################################################################################
  
  # Main simulation
  
  # Setup progress bar
  # progress bar ------------------------------------------------------------
  # requireNamespace(progress)
  
  pb <- progress_bar$new(
    format = "Iteration = :letter [:bar] :elapsed | expected time till finish: :eta",
    total = nsim,    # 100
    width = 120)

  progress_letter <- seq(1,nsim)  # token reported in progress bar

  # allowing progress bar to be used in foreach -----------------------------
  progress <- function(n){
    pb$tick(tokens = list(letter = progress_letter[n]))
  }

  opts <- list(progress = progress)
  
  
  s <- NULL
  ######
  if(method=="PP"){
    finalMatrix <- foreach::foreach(s=1:nsim, .combine=rbind, .packages = c("extraDistr", "fbst"), .options.snow = opts) %dopar% {
      tempMatrix = singleTrial_PP(s = s, n=nInit, responseMatrix = responseMatrix, nInit = nInit, Nmax = Nmax, batchsize = batchsize, a0 = a0, b0 = b0, p0 = p0, p1 = p1)
      
      tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    }
  }
  # First column stores binary variable stopped for futility (=1) or not (=0). Second column stores sample size n at which the trial stopped. Last column stores the PP value at stopping.
  # Columns 4+ include the PP values for nInit to NMax/batchsize interim analysis points. Are used for plotting later.
  
  if(method=="PPe"){
    refFunc = refFunc
    nu = nu
    shape1 = shape1
    shape2 = shape2 # use the parameters for reference function, evidence threshold and the beta-shapes for the beta reference function used in the simTrial function
    if(refFunc == "flat"){
      finalMatrix <- foreach::foreach(s=1:nsim, .combine=rbind, .packages = c("extraDistr", "fbst"), .options.snow = opts) %dopar% {
        tempMatrix = singleTrial_PPe(s = s, n=nInit, responseMatrix = responseMatrix, nInit = nInit, Nmax = Nmax, batchsize = batchsize, a0 = a0, b0 = b0, refFunc = "flat", p0 = p0, p1 = p1)
        
        tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
      }
    }
    if(refFunc == "beta"){
      finalMatrix <- foreach::foreach(s=1:nsim, .combine=rbind, .packages = c("extraDistr", "fbst"), .options.snow = opts) %dopar% {
        tempMatrix = singleTrial_PPe(s = s, n=nInit, responseMatrix = responseMatrix, nInit = nInit, Nmax = Nmax, batchsize = batchsize, a0 = a0, b0 = b0, refFunc = "beta", 
                                     shape1 = shape1, shape2 = shape2, p0 = p0, p1 = p1)
        
        tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
      }
    }
    if(refFunc == "binaryStep"){
      finalMatrix <- foreach::foreach(s=1:nsim, .combine=rbind, .packages = c("extraDistr", "fbst"), .options.snow = opts) %dopar% {
        tempMatrix = singleTrial_PPe(s = s, n=nInit, responseMatrix = responseMatrix, nInit = nInit, Nmax = Nmax, batchsize = batchsize, a0 = a0, b0 = b0, refFunc = "binaryStep", 
                                     shape1 = shape1, shape2 = shape2, truncation = truncation, p0 = p0, p1 = p1)
        
        tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
      }
    }
    if(refFunc == "relu"){
      finalMatrix <- foreach::foreach(s=1:nsim, .combine=rbind, .packages = c("extraDistr", "fbst"), .options.snow = opts) %dopar% {
        tempMatrix = singleTrial_PPe(s = s, n=nInit, responseMatrix = responseMatrix, nInit = nInit, Nmax = Nmax, batchsize = batchsize, a0 = a0, b0 = b0, refFunc = "relu", 
                                     shape1 = shape1, shape2 = shape2, truncation = truncation, p0 = p0, p1 = p1)
        
        tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
      }
    }
    if(refFunc == "palu"){
      finalMatrix <- foreach::foreach(s=1:nsim, .combine=rbind, .packages = c("extraDistr", "fbst"), .options.snow = opts) %dopar% {
        tempMatrix = singleTrial_PPe(s = s, n=nInit, responseMatrix = responseMatrix, nInit = nInit, Nmax = Nmax, batchsize = batchsize, a0 = a0, b0 = b0, refFunc = "palu", 
                                     shape1 = shape1, shape2 = shape2, truncation = truncation, p0 = p0, p1 = p1)
        
        tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
      }
    }
    if(refFunc == "lolu"){
      finalMatrix <- foreach::foreach(s=1:nsim, .combine=rbind, .packages = c("extraDistr", "fbst"), .options.snow = opts) %dopar% {
        tempMatrix = singleTrial_PPe(s = s, n=nInit, responseMatrix = responseMatrix, nInit = nInit, Nmax = Nmax, batchsize = batchsize, a0 = a0, b0 = b0, refFunc = "lolu", 
                                     shape1 = shape1, shape2 = shape2, truncation = truncation, p0 = p0, p1 = p1)
        
        tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
      }
    }
  }
  parallel::stopCluster(cl) #stop cluster
  # First column stores binary variable stopped for futility (=1) or not (=0). Second column stores sample size n at which the trial stopped. Last column stores the PP value at stopping.
  # Columns 4+ include the PP values for nInit to NMax/batchsize interim analysis points. Are used for plotting later.
  
  ##################################################################################################################################
  

  
  # for(t in 1:nrow(finalMatrix)){
  #   x = seq(nInit,Nmax,by=batchsize)
  #   y = finalMatrix[t,4:ncol(finalMatrix)]
  #   if(finalMatrix[t,3] < theta_L){
  #     resultPlot = resultPlot + geom_line(data=data.frame(x,y),aes(x=x,y=y),col=rgb(255, 0, 0, max = 255, alpha = 125))
  #   }
  #   if(finalMatrix[t,3] >= theta_U){
  #     resultPlot = resultPlot + geom_line(data=data.frame(x,y),aes(x=x,y=y),col=rgb(100, 149, 237, max = 255, alpha = 125))
  #   }
  #   if(finalMatrix[t,3] < theta_U && finalMatrix[t,ncol(finalMatrix)] > theta_L){
  #     resultPlot = resultPlot + geom_line(data=data.frame(x,y),aes(x=x,y=y),col=rgb(190, 190, 190, max = 255, alpha = 125))
  #   }
  # }
  
  
  ##################################################################################################################################

  
  #resultPlot = resultPlot +
  #labs(tag = paste0("Stopped for futility: ", percentageFutility*100, " %")) +    
  #  theme(plot.tag.position = c(0.25, 0.25)) 
  #resultPlot = resultPlot +
  #  annotate("text", x = 0.1, y = 1.05, label = paste0("Stopped for efficacy: ", percentageEfficacy*100, " %")) +
  #  coord_cartesian(ylim = c(4, 8), clip = "off")

  # dt = data.frame(x=rep(1,nrow(finalMatrix)),y=finalMatrix[,2])
  # #stripchart(dt$y ~ dt$x, vertical = FALSE, data = dt, 
  # #method = "jitter", add = TRUE, pch = 20, col = 'black', cex=1.5)
  # 
  #   #boxplot(dt$y ~ dt$x, col="skyblue",main="",horizontal=TRUE,ylab="",xlab="Trial size")
  # boxpl <- ggplot(dt, aes(x=x, y=y)) + 
  #   geom_boxplot(outlier.colour="black", outlier.shape=16,
  #            outlier.size=2, notch=FALSE, fill="cornflowerblue", color="black") +
  #   coord_flip() +
  #   geom_jitter(height = 0.25) +
  #   theme_bw() +
  #   labs(
  #     title = "",
  #     subtitle = "",
  #     caption = "",
  #     tag = "",
  #     x = "",
  #     y = "Trial sample size"
  #   ) +
  #   theme(axis.title.y=element_blank(),
  #       axis.text.y=element_blank())
  # 
  # library(ggpubr)
  # library(gridExtra)
  # grid.arrange(boxpl,resultPlot,nrow=2,heights=c(0.5,1))
  #finalPlot <- ggarrange(
  #  boxpl,                # First row with line plot
  # Second row with box and dot plots
  #  resultPlot, 
  #  nrow = 2,
  #  ncol = 1,
  #  widths = c(0.6, 1), heights = c(0.5, 1)
  #)
  #finalPlot
  
  ##################################################################################################################################
  
  # futility and efficacy percentages
  #percentageFutility = sum(finalMatrix[,ncol(finalMatrix)-1] < theta_L)/nsim
  #percentageEfficacy = sum(finalMatrix[,ncol(finalMatrix)-1] >= theta_U)/nsim # CAREFUL: We must use the second last column as otherwise the PP / PPe sum is empty! Thus, we check at the last interim analysis time point Nmax-batchsize whether the threshold for efficacy has been passed.
  percentageFutility = sum(finalMatrix[,1] == 1)/nsim 
  percentageEfficacy = sum(finalMatrix[,1] == 2)/nsim
  percentageInconclusive = sum(finalMatrix[,1] == 0)/nsim
  
  # cond1 = finalMatrix[,ncol(finalMatrix)-1] < theta_U
  # cond2 = finalMatrix[,ncol(finalMatrix)-1] > theta_L
  # cond = cond1 == cond2
  # if(is.null(nrow(finalMatrix[cond,]))){
  #   percentageInconclusive = 0
  # } else {
  #   percentageInconclusive = nrow(finalMatrix[cond,])/nsim
  # }

  
  #percentageInconclusive = sum(as.numeric(finalMatrix[,ncol(finalMatrix)-1] > theta_L & finalMatrix[,ncol(finalMatrix)-1] <theta_U))/nsim
  
  # return results as brada object
  if (is.null(refFunc)){
    refString = "Flat"
  } else {
    refString = "User-defined"
  }
  # return fbst object
  res = new("brada", data = list(a0 = a0,
                                b0 = b0,
                                Nmax = Nmax,
                                batchsize = batchsize,
                                nInit = nInit,
                                p_true = p_true,
                                p0 = p0,
                                p1 = p1,
                                theta_T = theta_T,
                                theta_L = theta_L,
                                theta_U = theta_U,
                                nsim = nsim,
                                seed = seed,
                                method = method,
                                refFunc = refFunc,
                                nu = nu,
                                shape1 = shape1,
                                shape2 = shape2,
                                truncation = truncation,
                                trials = finalMatrix,
                                futility = percentageFutility,
                                efficacy = percentageEfficacy,
                                inconclusive = percentageInconclusive,
                                expectedSampleSize = mean(finalMatrix[,2])))
  
  res
}

#' Generates trial data
#'
#' @method 
#' @export
#' @examples
## Function for data generation
generateData <- function(p,Nmax,nsim,seed=420){
  # Generate nsim*Nmax response values
  responses <- matrix(data = NA,nrow=nsim,ncol=Nmax)
  for(s in 1:nsim){
    for(n in 1:Nmax){
      responses[s,n] = rbinom(1,1,prob=p)
    }
  }
  responses
}


#' Computes power of the design
#'
#' @method 
#' @export
#' @examples
## Function for power calculation
power <- function(brada_object, p_true, nsim=100, cores = 2){
  power = brada(a0 = brada_object$a0,
                b0 = brada_object$b0,
                Nmax = brada_object$Nmax,
                batchsize = brada_object$batchsize,
                nInit = brada_object$nInit,
                p_true = p_true,
                p0 = brada_object$p0,
                p1 = brada_object$p1,
                theta_T = brada_object$theta_T,
                theta_L = brada_object$theta_L,
                theta_U = brada_object$theta_U,
                nsim = nsim,
                seed = brada_object$seed,
                method = brada_object$method,
                refFunc = brada_object$refFunc,
                nu = brada_object$nu,
                shape1 = brada_object$shape1,
                shape2 = brada_object$shape2,
                truncation = brada_object$truncation,
                cores = cores)
  power
}


#' Monitors a design
#'
#' @method 
#' @export
#' @examples
## Function for trial monitoring
monitor <- function(brada_object, obs){
  
  n=brada_object$nInit # first n patients
  x <- sum(obs) # number of successes observed
  n = length(obs)
  
  if(brada_object$method == "PP"){
    PP = 0 # initialize PP
    for(i in 0:(brada_object$Nmax-n)){ # iterate over possible results for successes
      marg=extraDistr::dbbinom(i,size=brada_object$Nmax-n,brada_object$a0+x,brada_object$b0+n-x) # P(Y=i|x)
      postProb=1-pbeta(brada_object$p1,brada_object$a0+x+i,brada_object$b0+brada_object$Nmax-x-i) # P(p>p1|X=x,Y=i) = evidence for H1
      indicator = 0
      if(postProb > brada_object$theta_T){
        indicator = 1
      }
      PP = PP + (marg*indicator)
    }
  }
  
  if(brada_object$method == "PPe"){
    PPe = 0 # initialize PP
    for(i in 0:(brada_object$Nmax-n)){ # iterate over possible results for successes
      marg=extraDistr::dbbinom(i,size=brada_object$Nmax-n,brada_object$a0+x,brada_object$b0+n-x) # P(Y=i|x)
      if(brada_object$refFunc=="beta"){
        FBET_H1 = fbst::fbet(posterior=dbeta, # analytic posterior version; increases computational speed
                             par_posterior=list(shape1=brada_object$a0+x+i,shape2=brada_object$b0+brada_object$Nmax-x-i),
                             interval=c(brada_object$p1,1),
                             FUN = dbeta,
                             par = list(shape1=brada_object$shape1,shape2=brada_object$shape2),
                             nu=brada_object$nu)
      }
      if(brada_object$refFunc=="flat"){
        FBET_H1 = fbst::fbet(posterior=dbeta, # analytic posterior version; increases computational speed
                             par_posterior=list(shape1=brada_object$a0+x+i,shape2=brada_object$b0+brada_object$Nmax-x-i),
                             interval=c(brada_object$p1,1),
                             FUN = dbeta, # beta(1,1) ref function is equal to flat reference function
                             par = list(shape1=1,shape2=1),
                             nu=brada_object$nu)
      }
      eValue = FBET_H1 # this is the evidence for H1:p>p1
      indicator = 0
      if(eValue > brada_object$theta_T){
        indicator = 1
      }
      PPe = PPe + (marg*indicator)
    }
  }
   
  if (brada_object$method == "PP"){
    designString = "Trial design: Predictive probability design"
  }
  if (brada_object$method == "PPe"){
    designString = "Trial design: Predictive evidence value design"
  }
  
  message("--------- BRADA TRIAL MONITORING ---------\nPrimary endpoint: binary\nTest of H_0: p <= ",brada_object$p0, " against H_1: p > ", brada_object$p1, "\n", designString, "\nMaximum sample size: ", brada_object$Nmax, "\nFirst interim analysis at: ", brada_object$nInit, "\nInterim analyses after each ", brada_object$batchsize, " observations \nLast interim analysis at: ", brada_object$Nmax-brada_object$batchsize, " observations \n-----------------------------------------\nCurrent trial size: ", n, " patients \n--------------- RESULTS -----------------\n")

  if(!(n %in% seq(brada_object$nInit,brada_object$Nmax,by=brada_object$batchsize))){
    message("WARNING: YOU ARE VIOLATING THE TRIAL PROTOCOL BY RUNNING AN INTERIM ANALYSIS NOW \n")
  }
  
  if(brada_object$method == "PP"){
    if(PP < brada_object$theta_L){ # check for futility
        message("Predictive probability of trial success: ", round(PP,5), "\n Futility threshold: ", brada_object$theta_L, "\nDecision: Stop for futility")
    }
    if(PP > brada_object$theta_U){ # check for efficacy
      message("Predictive probability of trial success: ", round(PP,5), "\nEfficacy threshold: ", brada_object$theta_U, "\nDecision: Stop for efficacy")
    }
    if(PP < brada_object$theta_U & PP > brada_object$theta_L){ # check for efficacy
      message("Predictive probability of trial success: ", round(PP,5), "\nEfficacy threshold: ", brada_object$theta_U, "\nDecision: Continue the trial")
    }
  }

  if(brada_object$method == "PPe"){
    if(PPe < brada_object$theta_L){ # check for futility
      message("Predictive probability of trial success:", round(PPe,5), "\nFutility threshold:v", brada_object$theta_L, "\nDecision: Stop for futility")
    }
    if(PPe > brada_object$theta_U){ # check for efficacy
      message("Predictive probability of trial success:", round(PPe,5), "\nEfficacy threshold:v", brada_object$theta_U, "\nDecision: Stop for efficacy")
    }
    if(PPe < brada_object$theta_U & PPe > brada_object$theta_L){ # check for efficacy
      message("Predictive probability of trial success:", round(PPe,5), "\nEfficacy threshold:v", brada_object$theta_U, "\nDecision: Continue the trial")
    }
  }
}

#' brada class
#'
#' Stores the results of a Full Bayesian Significance Test
#'
#' @slot data A named list for storing the user-accessible data of an fbst object
#' a0 A numeric shape1 parameter of the beta prior
#' b0 A numeric shape2 parameter of the beta prior
#' Nmax The numeric maximum trial size
#' batchsize The numeric batchsize for interim analyses
#' nInit The numeric minimum sample size at which first interim analysis is performed
#' p_true The numeric true response probability
#' p0 The numeric right boundary of null hypothesis
#' p1 The numeric left boundary of alternative hypothesis
#' theta_T The numeric threshold for PP and PPe designs
#' theta_L The numeric threshold for futility
#' theta_U The numeric threshold for efficacy
#' nsim The numeric number of Monte Carlo simulations 
#' seed The numeric seed for reproducibility
#' method The string for method, either "PP" or "PPe"
#' refFunc The string for reference function, either "flat" or "beta"
#' nu The numeric evidence threshold for Bayesian evidence interval
#' shape1 The numeric shape1 parameter for beta reference function
#' shape2 The numeric hape2 parameter for beta reference function
#' trials A Matrix containing the simulated trials
#' futility Percentage of the trials which stopped for futility
#' efficacy Percentage of the trials which stopped for efficacy
#' inconclusive Percentage of the trials which remained inconclusive at last interim analysis
#' expectedSampleSize Expected sample size of the design
#' @name brada-class
#' @rdname brada-class
#' @export
setClass("brada", representation(data="list"),
         prototype = NULL,
         validity = function(object) return(TRUE)
)

#' plot object of class brada
#' @usage \\method{plot}{brada}(x, ...)
#' @export
plot.brada <- function(x, trajectories = 100, ...){
  #par_temp = par() 
  partemp <- par(no.readonly = TRUE) # save graphics parameters to restore later
  on.exit(par(partemp))
  
  layout.matrix <- matrix(c(2, 1), nrow = 2, ncol = 1)
  layout(mat = layout.matrix,
         heights = c(1, 4), # Heights of the two rows
         widths = c(2,1)) # Heights of the two columns
  if(x$method=="PP"){
    plot(1, type="n", xlab="n", ylab=expression(PP), xlim=c(x$nInit, x$Nmax), ylim=c(0, 1),xaxt='n')
  }
  if(x$method=="PPe"){
    plot(1, type="n", xlab="n", ylab=expression(PP[e]), xlim=c(x$nInit, x$Nmax), ylim=c(0, 1),xaxt='n')
  }
  axis(1, at=seq(x$nInit,x$Nmax,by=x$batchsize), labels=seq(x$nInit,x$Nmax,by=x$batchsize))
  abline(h=x$theta_U,col="black",lwd=1.5,lty=1)
  abline(h=x$theta_L,col="black",lwd=1.5,lty=1)
  abline(v=x$Nmax,col="black",lty=3,lwd=2)
  abline(v=x$Nmax-x$batchsize,col="purple",lty=3,lwd=2) # last interim analysis point where PP / PPe value makes sense
  
  # plotting of trajectories (we replace the last PP / PPe value at Nmax with the second last at Nmax-batchsize; this is the last
  # time point where a PP / PPe value makes sense as a predictive probability)
  trajec = x$trials[,1:(ncol(x$trials)-1)]
  trajec = cbind(trajec,x$trials[,ncol(x$trials)-1])
  limit = trajectories
  if(limit > nrow(trajec)){
    limit = nrow(trajec)
  }
  for(t in 1:limit){
    if(trajec[t,ncol(trajec)] < x$theta_L){
      lines(seq(x$nInit,x$Nmax,by=x$batchsize),trajec[t,3:ncol(trajec)],ty="l",col=rgb(255, 0, 0, maxColorValue=255, alpha = 125))
    }
    if(trajec[t,ncol(trajec)] >= x$theta_U){
      lines(seq(x$nInit,x$Nmax,by=x$batchsize),trajec[t,3:ncol(trajec)],ty="l",col=rgb(100, 149, 237, maxColorValue=255, alpha = 125))
    }
    if(trajec[t,ncol(trajec)] < x$theta_U && trajec[t,ncol(trajec)] > x$theta_L){
      lines(seq(x$nInit,x$Nmax,by=x$batchsize),trajec[t,3:ncol(trajec)],ty="l",col=rgb(190, 190, 190, maxColorValue=255, alpha = 125))
    }
  }
  
  # compute Monte Carlo standard errors
  #valuesAtLastInterimAnalysis = x$trials[1:x$nsim,ncol(x$trials)-1]
  #futilityVector = as.numeric(valuesAtLastInterimAnalysis < x$theta_L)
  futilityVector = efficacyVector = inconclusiveVector = rep(0,x$nsim)
  futilityVector[which(x$trials[,1] == 1)] = 1
  efficacyVector[which(x$trials[,1] == 2)] = 1
  inconclusiveVector[which(x$trials[,1] == 0)] = 1
  #efficacyVector = as.numeric(valuesAtLastInterimAnalysis >= x$theta_U)
  #inconclusiveVector = as.numeric(valuesAtLastInterimAnalysis < x$theta_U & valuesAtLastInterimAnalysis > x$theta_L)
  sampleSizeVector = x$trials[1:x$nsim,2]
  # Bootstrap MCSEs
  bootstrap_se_futility = round(mean(replicate(10000, sd(sample(futilityVector, replace=TRUE))/sqrt(length(futilityVector)))),3)
  bootstrap_se_efficacy = round(mean(replicate(10000, sd(sample(efficacyVector, replace=TRUE))/sqrt(length(efficacyVector)))),3)
  bootstrap_se_inconclusive = round(mean(replicate(10000, sd(sample(inconclusiveVector, replace=TRUE))/sqrt(length(inconclusiveVector)))),3)
  bootstrap_se_samplesize = round(mean(replicate(10000, sd(sample(sampleSizeVector, replace=TRUE))/sqrt(length(sampleSizeVector)))),3)
  
  # add text
  mtext(paste0("Stopped for futility: ", round(x$futility*100,2), " % [ \u00B1", bootstrap_se_futility*100, "% ]"), side=1, line=2.5, at=9*x$Nmax/10)
  mtext(paste0("Stopped for efficacy: ", round(x$efficacy*100,2), " % [ \u00B1", bootstrap_se_efficacy*100, "% ]"), side=3, line=0.5, at=9*x$Nmax/10)
  mtext(paste0("Inconclusive: ", round(x$inconclusive*100,2), " % [ \u00B1", bootstrap_se_inconclusive*100, "% ]"), side=4, line=0.5, at=0.5)
  
  # add expected sample size + Monte Carlo SE to plot
  mtext(paste0("Expected sample size: ", round(x$expectedSampleSize,2), " [ \u00B1", bootstrap_se_samplesize, "]"), side=3, line=0.5, at=4*x$Nmax/10)
  
  # add boxplot of sample size
  par(mar = c(0.05, 3.5, 0, 2))
  boxplot(x$trials[,2], bty = "n", xlab="Sample size",
          col = "skyblue", frame = FALSE, horizontal = TRUE)
  # par(mar = par_temp$mar) # reset the margins to the old par margins
}

#' Calibrate an object of class brada
#' @usage \\method{calibrate}(brada_object, ...)
#' @export
calibrate <- function(brada_object, nsim = 100, cores = 2, seq, alpha=NULL, beta=NULL, calibration = "nu"){
  upperbound <- NULL
  rate = character(length(seq))
  if(calibration == "nu"){
    message(cli::symbol$info, " Starting to calibrate evidence threshold...")
    for(i in 1:length(seq)){
      res = brada(a0 = brada_object$a0, b0 = brada_object$b0,
                  Nmax = brada_object$Nmax,
                  batchsize = brada_object$batchsize,
                  nInit = brada_object$nInit,
                  p_true = brada_object$p_true,
                  p0 = brada_object$p0,
                  p1 = brada_object$p1,
                  theta_T = brada_object$theta_T,
                  theta_L = brada_object$theta_L,
                  theta_U = brada_object$theta_U,
                  nsim = nsim,
                  method = brada_object$method,
                  refFunc = brada_object$refFunc,
                  nu = seq[i],
                  cores = cores)
      
      # Monte Carlo Standard error
      #efficacyVector = rep(0,res$nsim)
      #efficacyVector[which(res$trials[,1] == 2)] = 1
      #bootstrap_se_efficacy = round(mean(replicate(10000, sd(sample(efficacyVector, replace=TRUE))/sqrt(length(efficacyVector)))),3)
      
      upperbound = res$efficacy #+ bootstrap_se_efficacy
      rate[i] = as.character(res$efficacy)
      message(cli::symbol$tick, " Iteration ", i, " of ", length(seq), " completed. False-positive rate: ", upperbound)
      if(upperbound <= alpha){
        message("\n", cli::symbol$info, " Evidence threshold is calibrated.\nFalse-positive rate is ", upperbound, " at last iteration.\nProceed with evidence threshold nu = ", seq[i], ".")
        if(upperbound == 0){
          message("\n Warning, false-positive rate has reduced to zero, you should probably increase the maximum trial size.")
        }
        break()
      }
      
    }
  }
  
  if(calibration == "theta_L"){
    message(cli::symbol$info, " Starting to calibrate futility threshold...")
    for(i in 1:length(seq)){
      res = brada(a0 = brada_object$a0, b0 = brada_object$b0,
                  Nmax = brada_object$Nmax,
                  batchsize = brada_object$batchsize,
                  nInit = brada_object$nInit,
                  p_true = brada_object$p_true,
                  p0 = brada_object$p0,
                  p1 = brada_object$p1,
                  theta_T = brada_object$theta_T,
                  theta_L = seq[i],
                  theta_U = brada_object$theta_U,
                  nsim = nsim,
                  method = brada_object$method,
                  refFunc = brada_object$refFunc,
                  nu = brada_object$nu,
                  cores = cores)
      
      # Monte Carlo Standard error
      #futilityVector = rep(0,res$nsim)
      #futilityVector[which(res$trials[,1] == 1)] = 1
      #bootstrap_se_futility = round(mean(replicate(10000, sd(sample(futilityVector, replace=TRUE))/sqrt(length(futilityVector)))),3)
      
      upperbound = res$futility #+ bootstrap_se_futility
      rate[i] = as.character(upperbound)
      message(cli::symbol$tick, " Iteration ", i, " of ", length(seq), " completed. False-negative rate: ", upperbound)
      if(upperbound <= beta){
        message(cli::symbol$info, " Futility threshold is calibrated.\nFalse-negative rate is ", upperbound, " at last iteration.\nProceed with futility threshold theta_L = ", seq[i], ".")
        if(upperbound == 0){
          message("\n Warning, false-negative rate has reduced to zero, you should probably increase the maximum trial size.")
        }
        break()
      }
    }
  }
  
  if(calibration == "nu"){
    if(upperbound >= alpha){
      message(cli::symbol$cross, " Calibration failed. Proceed with ", seq[length(seq)], " as the new minimum for the evidence threshold.")
    }
  }
  if(calibration == "theta_L"){
    if(upperbound >= beta){
      message(cli::symbol$cross," Calibration failed. Proceed with ", seq[length(seq)], " as the new maximum for the futility threshold.")
    }
  }
  upperbound = as.numeric(upperbound)
  upperbound
}

#' Print summary of an object of class brada
#' @usage \\method{summary}{brada}(x, ...)
#' @export
summary.brada <- function(object, ...){
  cat("BAYESIAN RESPONSE-ADAPTIVE DESIGN ANALYSIS:\n")
  cat("--------------------------------------\n")
  cat("Primary endpoint: binary \n")
  cat("Test of H_0: p <=",object$p0," against H_1: p >", object$p1, "\n")
  
  cat("Maximum sample size:", object$Nmax, "\n")
  cat("First interim analysis at:", object$nInit, "\n")
  cat("Interim analyses after each", object$batchsize, "observations \n")
  cat("Last interim analysis at:", object$Nmax-object$batchsize, "observations \n")
  cat("\n")
  cat("------------- RESULTS ---------------\n")
  cat("Expected number of patients:", object$expectedSampleSize, "\n")
  cat("Probability to stop for futility:", object$futility, "\n")
  cat("Probability to stop for efficacy:", object$efficacy, "\n")
  cat("Probability that PP / PPe is inconclusive at last interim analysis:", object$inconclusive, "\n")
}



#' Show method for an object of class brada
#' @usage \\method{show}{brada}(object)
#' @export
show.brada <- function(object){
  message("BAYESIAN RESPONSE-ADAPTIVE DESIGN ANALYSIS:\n")
  message("--------------------------------------\n")
  message("Primary endpoint: binary \n")
  message("Test of H_0: p <=",object$p0," against H_1: p >", object$p1, "\n")
  
  message("Maximum sample size:", object$Nmax, "\n")
  message("First interim analysis at:", object$nInit, "\n")
  message("Interim analyses after each", object$batchsize, "observations \n")
  message("Last interim analysis at:", object$Nmax-object$batchsize, "observations \n")
  message("\n")
  message("------------- RESULTS ---------------\n")
  message("Expected number of patients:", object$expectedSampleSize, "\n")
  message("Probability to stop for futility:", object$futility, "\n")
  message("Probability to stop for efficacy:", object$efficacy, "\n")
  message("Probability that PP / PPe is inconclusive at last interim analysis:", object$inconclusive, "\n")
}



#' Access results stored in the data slot of an object of class brada
#' @usage \\method{$}{brada}(x)
#' @export
setMethod('$', signature="brada",
          definition=function(x, name) {
            return_value = x@data[[name]]
            names(return_value) = name
            return(return_value)}
)




#' Get names of the data slot of an object of class brada
#' @usage \\method{names}{brada}(object)
#' @export
names.brada <- function(x){
  names(x@data)
}