##
##################################################################################
## version 1.0
##################################################################################
##
##################################################################################
## Required functions ##
##################################################################################
##
## optDesing for optimal sample size and position
##
runOpt <- function(design.matrix, coords, response=NULL,
                   initSampleSize = 10,
                   increment.factor = 5,
                   threshold = 0.05,
                   optType = "exact",
                   model = "univariate",
                   weight.response = NULL,
                   store.para = TRUE,
                   seed=1234){
  ##
  ## design.matrix = n x p in matrix/data format
  ## coords = n x 2 coordinates in matrix/data format
  ## response = n x q matrix, q = 1 if univariate
  ## initSampleSize = number of initial samples to select
  ## increment.factor = number of increments in each iteration
  ## threshold = threshold value for convergence of the rate of change
  ## optType = can take 2 arguments: (i) "exact" for exact algorithm and (ii) "montecarlo" for montecarlo algorithm
  ## model = can take "univariate" or "multivariate"
  ## weight.response = optional for multivariate model prediction, if null then use equal weights for multivariate response
  ## seed = a random seed number to reporduce the results
  ##
  start.time <- proc.time()[3]
  set.seed(seed)
  coords <- as.data.frame(coords)
  design.matrix <- as.data.frame(design.matrix)
  design.matrix <- data.frame(lapply(design.matrix, function(x) scale(x)))
  n <- nrow(coords)
  p <- ncol(design.matrix)
  if(nrow(design.matrix)<=initSampleSize){
    stop(paste0("value for 'initSampleSize' = ",initSampleSize," should be less than the sample size n = ",n))
  }
  if(!model%in%c("univariate","multivariate")){
    stop(" model argument can only take 'univariate' or 'multivariate'.\n")
  }
  if(!optType%in%c("exact","montecarlo")){
    stop(" optType argument can only take 'exact'or 'montecarlo'.\n")
  }
  if(initSampleSize<=p){
    print(paste0("initSampleSize <= p; i.e.,",initSampleSize," <= ",p))
    print(paste0("initSampleSize is considered as: ",p+1))
  }
  ## run a repeat loop
  sample.size <- initSampleSize
  vall <- 0.00000001 # initial criteria for threshold
  j <- 0
  sample.size0 <- sample.size
  newd <- run_fnc_opt(design.matrix=design.matrix,optSample=sample.size0,coords=coords,type=optType,seed=seed)
  ##
  if(unique(newd$flag)%in%1){
    repeat{
      j <- j+1
      sample.size <- sample.size0 + increment.factor
      sample.size0 <- sample.size
      newd <- run_fnc_opt(design.matrix=design.matrix,optSample=sample.size0,coords=coords,type=optType,seed=seed)
      if (unique(newd$flag)%in%0 | sample.size0 >= nrow(dat)){
        break
      }
    }
  }
  rm(j);
  ##
  ## output with response variable
  ##
  if(!is.null(response)){
    ##
    ## univariate model
    ##
    xvar <- names(design.matrix)
    dat <- design.matrix
    dat$id <- 1:nrow(dat)
    ##
    if(model%in%c("univariate")){
      if(!length(response)%in%n){
        stop(" you are running univariate model \n use vector input and/or check response variable \n ")
      }
      dat$response <- response
      yvar <- "response"
      q <- 1 # univariate case
      f <- as.formula(paste(paste("response")," ~ ", paste(xvar,collapse= "+")))
      m <- lm(f, data = dat[dat$id%in%newd$id,])
      pr <- predict(m, newd = dat)
    }
    ##
    if(model%in%c("multivariate")){
      if(!(is.matrix(response)|is.data.frame(response))){
        stop(" you are running multivariate model \n use matrix (n x q) input and/or check response variable \n ")
      }
      if(!nrow(response)%in%n){
        stop(" check the (n x q) matrix input for response. \n")
      }
      q <- ncol(response)
      response <- as.data.frame(response)
      names(response) <- paste0("response",1:q)
      yvar <- paste0("response",1:q)
      dat <- cbind(dat,response)
      #paste(paste0("response",1:q),collapse= ",")
      ff <- paste0("cbind(",paste(paste0("response",1:q),collapse= ","),")")
      f <- as.formula(paste(paste(ff)," ~ ", paste(xvar,collapse= "+")))
      m <- lm(f, data = dat[dat$id%in%newd$id,])
      pr <- predict(m, newd = dat)
    }
    ##
    if(model%in%c("univariate")){
      dd <- na.omit(cbind(dat[,yvar],pr))
      val <- abs(cor(dd)[1,2])^2
      val0 <- val
      val.para <- spT.validation(dd[,1],dd[,2])
    }
    ##
    if(model%in%c("multivariate")){
      val <- c()
      val.para <- list()
      for(k in 1:q){
        dd <- na.omit(cbind(dat[,yvar],pr))
        val[k] <- abs(cor(dd)[1,2])^2
        val.para[[k]] <- spT.validation(dd[,1],dd[,2])
      }
      if(!is.null(weight.response)){
        if(!length(weight.response)%in%q){
          stop(" check the weight.response argument.\n")
        }
        val0 <- weighted.mean(val,weight.response)
      }
      else{
        val0 <- mean(val,na.rm=TRUE)
      }
    }
    ##
    roc <- (val0-vall)/vall # rate of change
    j <- 0
    #print(paste0("Initialisation: ",j,"; Rate of Change: ",round(roc,4)))
    #print(c(j,vall,val0,roc))
    list.para <- list()
    repeat{
      j <- j+1
      sample.size <- sample.size0 + increment.factor
      if (roc <= threshold | sample.size >= nrow(dat)){
        break
      }
      vall <- val0
      sample.size0 <- sample.size
      newd0 <- run_fnc_opt(design.matrix=design.matrix,optSample=sample.size0,coords=coords,type=optType,seed=seed)
      ##
      k <- 0
      if(unique(newd0$flag)%in%1){
        repeat{
          k <- k+1
          sample.size <- sample.size0 + increment.factor
          sample.size0 <- sample.size
          newd0 <- run_fnc_opt(design.matrix=design.matrix,optSample=sample.size0,coords=coords,type=optType,seed=seed)
          cat(k)
          if (unique(newd0$flag)%in%0 | sample.size0 >= nrow(dat)){
            break
          }
        }
      }
      ##
      m <- lm(f, data = dat[dat$id%in%newd0$id,])
      pr <- predict(m, newdata = dat)
      ##
      if(model%in%c("univariate")){
        dd <- na.omit(cbind(dat[,yvar],pr))
        val <- abs(cor(dd)[1,2])^2
        val0 <- val
        val.para <- spT.validation(dd[,1],dd[,2])
      }
      ##
      if(model%in%c("multivariate")){
        val <- c()
        val.para <- list()
        for(k in 1:q){
          dd <- na.omit(cbind(dat[,yvar[k]],pr[,k]))
          val[k] <- abs(cor(dd)[1,2])^2
          val.para[[k]] <- spT.validation(dd[,1],dd[,2])
        }
        if(!is.null(weight.response)){
          if(!length(weight.response)%in%q){
            stop(" check the weight.response argument.\n")
          }
          val0 <- weighted.mean(val,weight.response)
        }
        else{
          val0 <- mean(val,na.rm=TRUE)
        }
      }
      ##
      roc <- (val0-vall)/vall # rate of change
      ##
      #print(paste0("Itr: ",j,"; Rate of Change: ",round(roc,4)))
      #print(paste0("Itr: ",j,"; Rate of Change: ",round(roc,4),"; Sample: ",nrow(newd0)))
      print(paste0("Sample Size: ",nrow(newd0),"; Rate of Change: ",round(roc,4)))
      #print(paste("replication:",j,"pre-r2:",vall,"curr-r2:",val0,"roc:",roc))
      list.para[[j]] <- list(roc=roc,pre_r2=vall,curr_r2=val0,other.para=val.para)
      ##
      if (roc > threshold){
        newd <- newd0
      }
      ##
    }
    ##
    if(model%in%c("univariate")){
      new.response <-response[newd$id]
    }
    if(model%in%c("multivariate")){
      new.response <-response[newd$id,]
    }
    ##
    end.time <- proc.time()[3]
    t <- end.time-start.time
    comp.time <- .fnc.time_(t)
    if(isTRUE(store.para)){
      out <- list(model=model,optType=optType,coords=coords,optCoords=newd,
                  optID=sort(newd$id),new.design.matrix=design.matrix[newd$id,],
                  new.response = new.response, comp.time=comp.time,para.list=list.para)
    }
    else{
      out <- list(model=model,optType=optType,coords=coords,optCoords=newd,
                  optID=sort(newd$id),new.design.matrix=design.matrix[newd$id,],
                  new.response = new.response, comp.time=comp.time)
    }
    class(out) <- "hcas"
    out
  }
  ##
  ## output without considering response variable
  ##
  else{
    end.time <- proc.time()[3]
    t <- end.time-start.time
    comp.time <- .fnc.time_(t)
    out <- list(model=model,optType=optType,coords=coords,optCoords=newd,
                optID=sort(newd$id),new.design.matrix=design.matrix[newd$id,],
                new.response = NULL, comp.time=comp.time)
    class(out) <- "hcas"
    out
  }
}
##
## optDesing for optimal pixel position using fixed sample size
##
run_fnc_opt <- function(design.matrix,optSample,coords,type,seed){
  ##
  ## design.matrix = n x p in matrix format
  ## optSample = number of samples to select
  ## seed = a random seed number to reporduce the results
  ## type = can take 2 arguments: (i) "exact" for exact algorithm and (ii) "montecarlo" for montecarlo algorithm
  ##
  ## criterion = always "I" to represent predictive performance, see details in "AlgDesign"
  criterion = "I"
  ##
  n <- nrow(coords)
  if(nrow(coords)!=nrow(design.matrix)){
    stop("coords and design.matrix should have same number of rows.\n")
  }
  ##
  set.seed(seed)
  if(type=="exact"){
    out <- try(optFederov(data = as.matrix(design.matrix), nTrials = optSample, criterion=criterion), TRUE)
  }
  else if(type=="montecarlo"){
    out <- try(optMonteCarlo(data = as.matrix(design.matrix), nTrials = optSample, criterion=criterion), TRUE)
  }
  else{
    stop("type argument can only take 'exact' and ''")
  }
  if(class(out) == "try-error"){
    set.seed(seed*round(runif(1,1,10)))
    if(type=="exact"){
      out <- try(optFederov(data = as.matrix(design.matrix), nTrials = optSample, criterion=criterion), TRUE)
    }
    else if(type=="montecarlo"){
      out <- try(optMonteCarlo(data = as.matrix(design.matrix), nTrials = optSample, criterion=criterion), TRUE)
    }
    else{
      stop("type argument can only take 'exact' and ''")
    }
    if(class(out) == "try-error"){
      newd <- as.data.frame(coords)
      newd$id <- 1:n
      names(newd) <- c("x","y","id")
      newd$flag <- 0
      if(nrow(newd)>optSample){
        set.seed(seed)
        newd <- newd[sample(1:nrow(newd),optSample),]
        newd$flag <- 1 # 1 refers to random selection
      }
    }
    else{
      newd <- as.data.frame(coords)
      newd$id <- 1:n
      names(newd) <- c("x","y","id")
      newd <- newd[out$rows,]
      newd <- data.frame(newd,proportion=out$design[,1])
      newd <- newd[order(newd$proportion,decreasing = TRUE),]
      newd$flag <- 0 # 0 refers to "optimal" selection
      newd$proportion <- NULL
    }
  }
  else{
    newd <- as.data.frame(coords)
    newd$id <- 1:n
    names(newd) <- c("x","y","id")
    newd <- newd[out$rows,]
    newd <- data.frame(newd,proportion=out$design[,1])
    newd <- newd[order(newd$proportion,decreasing = TRUE),]
    newd$flag <- 0 # 0 refers to "optimal" selection
    newd$proportion <- NULL
  }
  row.names(newd) <- NULL
  newd
}
##
## plot function
##
plot.hcas <- function(x, ...){
  plot(x$coords,pch="*",xlab="x",ylab="y")
  points(x$optCoords[,1:2], ...)
  title(sub=paste0("Total observation: ",nrow(x$coords),
                   ", Optimal sample: ",nrow(x$optCoords)))
}
##
## print function
##
print.hcas <-function(x, ...){
  cat("-------------------------"); cat('\n');
  cat("Model: "); cat(x$model); cat('\n');
  cat("Computation time: "); cat(x$comp.time); cat("\n")
  cat("-------------------------"); cat('\n');
}
##
## convert seconds into min. hour. and day
##
.fnc.time_<-function(t)
{
  #
  if(t < 60){
    t <- round(t,2)
    tt <- paste(t," - Sec.")
    cat(paste("##\n# Elapsed time:",t,"Sec.\n##\n"))
  }
  #
  if(t < (60*60) && t >= 60){
    t1 <- as.integer(t/60)
    t <- round(t-t1*60,2)
    tt <- paste(t1," - Mins.",t," - Sec.")
    cat(paste("##\n# Elapsed time:",t1,"Min.",t,"Sec.\n##\n"))
  }
  #
  if(t < (60*60*24) && t >= (60*60)){
    t2 <- as.integer(t/(60*60))
    t <- t-t2*60*60
    t1 <- as.integer(t/60)
    t <- round(t-t1*60,2)
    tt <- paste(t2," - Hour/s.",t1," - Mins.",t," - Sec.")
    cat(paste("##\n# Elapsed time:",t2,"Hour/s.",t1,"Min.",t,"Sec.\n##\n"))
  }
  #
  if(t >= (60*60*24)){
    t3 <- as.integer(t/(60*60*24))
    t <- t-t3*60*60*24
    t2 <- as.integer(t/(60*60))
    t <- t-t2*60*60
    t1 <- as.integer(t/60)
    t <- round(t-t1*60,2)
    tt <- paste(t3," - Day/s.",t2," - Hour/s.",t1," - Mins.",t," - Sec.")
    cat(paste("##\n# Elapsed time:",t3,"Day/s.",t2,"Hour/s.",t1,"Mins.",t,"Sec.\n##\n"))
  }
  #
  tt
}
