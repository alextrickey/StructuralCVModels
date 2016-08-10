# Filename: scvm_functions.r
# Purpose: 
#   Define function to be used to fit CV models

# List:
#   save_rgn: save a seed for random number generation
#   restore_rgn: restore a saved seed
#   make_data: make a sample data set which can be used for testing
#   getstart(MOD,S): Calculates Starting Values
#   SIGMAof(theta): Calculated Covariance of Free Parameters
#   vechSIGMA(theta): Returns half-vectorize of SIGMAof(theta)
#   calcDelta(theta,W): Calculates "Delta" for standard errors
#   getDelta(theta,W): Error Checking Wrapper for calcDelta
#   Criterion Functions: Return crit to be minimized given theta
#     fADFcv, fADFcvn, fGLScv, fMLcv
#   plzcon(theta,func): Optimizes Criterion Func to Estimate Pars
#   printRES: Prints results in form set up by plzcon and SEs

#######################################################################

#Saves current seed
#ref: http://stackoverflow.com/questions/13997444/
save_rng <- function(savefile=tempfile()) {
  if (exists(".Random.seed"))  {
    oldseed <- get(".Random.seed", .GlobalEnv)
  } else stop("don't know how to save before set.seed() or r*** call")
  oldRNGkind <- RNGkind()
  save("oldseed","oldRNGkind",file=savefile)
  invisible(savefile)
}

#Restores saved seed
#ref: http://stackoverflow.com/questions/13997444/
restore_rng <- function(savefile) {
  load(savefile)
  do.call("RNGkind",as.list(oldRNGkind))  ## must be first!
  assign(".Random.seed", oldseed, .GlobalEnv)
}


#Function to Generate Test Data with CV Structure
#   model_file_name 
#     - Name of the model file 
#   model_dir 
#     - Path to directory containing model file and seed file (if restoring)
#   save_data 
#     - Set to TRUE to save sample data to csv
#   use_random_seed 
#     - Set to TRUE to make or load seed file
#   existing_seed_file 
#     - Name of file containing seed to be restored
make_data <- function(model_file_name = 'MOD.dat',
                      model_dir = getwd(),
                      save_data = FALSE,
                      use_random_seed = FALSE,
                      existing_seed_file = NA
){
  
  #Set working directory
  setwd(model_dir)
  
  ########################
  ##Load and Prep Model ##
  ########################
  
  #Load MOD file (Must be copied to save folder)
  MOD <- read.delim(model_file_name)
  
  #Define Population Factor Structure (PSYID = 1)
  pF1  = 3 #Nodes per factor (pF1=pF2)
  pOBS = 6
  Fcov = 0.3 #Factor Covariance (yields cor(F1,F2)=0.3)
  OffBlocks = matrix(Fcov,pF1,pF1)
  DiagBlocks = matrix(1,pF1,pF1) + diag(pF1)
  PSY = rbind(cbind(DiagBlocks,OffBlocks),cbind(OffBlocks,DiagBlocks))
  
  #Define Population Mean and Covariance Matrices
  MU = rep(1,pOBS)
  SIGMA = diag(MU) %*% PSY %*% diag(MU)
  
  #Set sample size
  N =1000
  
  ###################################
  ## Count and Selection Variables ##
  ###################################
  
  #total number of variables 
  p = max(c(MOD$from,MOD$to))
  # ***Assumes no unused variables
  
  #number of observed variables
  pOBS = nrow(SIGMA) 
  # ***Assumes variables 1:pOBS are all observed and included in model
  
  #Number of unique covariances
  pSTAR = pOBS*(pOBS+1)/2
  print(pSTAR)
  
  #Select/Identify Paths and Vars
  sel.ip  <- (MOD$code == 1)#paths from ivs
  sel.dp  <- (MOD$code == 2)#paths between dvs
  sel.cov <- (MOD$code == 3)#vars and covs of ivs
  
  #Identify and count (strict) IVs 
  listIV = unique(MOD$from[sel.ip])
  pIV = length(listIV) #includes errors
  print(pIV)
  
  #Identify and count (all) DVs
  listDV = unique(MOD$to[MOD$code != 3])
  pDV = length(listDV)
  print(pDV)
  
  #Select/Identify Free Parameters
  sel.fixed <- (MOD$index == 0)
  sel.free <- !(sel.fixed)
  
  #number of parameters
  qTOT = nrow(MOD)
  qFREE = sum(sel.free)
  print(qFREE)
  # *** Assumes no constraints on parameters
  
  #df for chi-squared
  df = pSTAR - qFREE
  print(df)
  
  ##############
  ## Set Seed ##
  ##############
  
  ##Generate/Restore Seed for Random Number Generator
  if(use_random_seed) {
    #If no existing seed file create new
    if(is.na(existing_seed_file)){
      existing_seed_file = paste("RandomSeed_",
                       format(Sys.time(),"%Y%b%d_%H%M"),
                       ".Rdata",sep="")
      save_rng(savefile = existing_seed_file)
    }
    #Load Seed
    restore_rng(savefile = existing_seed_file)
  }
  
  #################
  ## Draw Sample ##
  #################
  Y=mvrnorm(N,MU,SIGMA)
  
  if(save_data == TRUE){
    write.csv(Y, paste('cv_sample_data.csv'))
  }
  
  #Return Sample Data Matrix
  Y
}


#######################################################################

#Function to Store Data (Must Match Col Names Below)
storeRES <- function(rMETH,filename="") {
  cat(c(N,REP,rMETH$METH,rMETH$conv,rMETH$oalg,rMETH$mini,rMETH$pars,rMETH$SEs),file=filename,append=T)
  cat("\n",file=filename,append=T)
}

#######################################################################

#Defines starting values for parameter estimation based on the specified
#model (MOD)
getstart <- function(MOD,S) {
  
  #Get selection Variables and Change Format
  MOD$value = as.character(MOD$value)
  sel.na <- (MOD$value == 'X')
  sel.ev <- (MOD$value == 'E')
  MOD$value = suppressWarnings(as.numeric(MOD$value))
  
  #number of observed variables
  pOBS = nrow(S) 
  # ***Assumes variables 1:pOBS are all observed and included in model
  
  #Select/Identify Paths and Vars
  sel.ip  <- (MOD$code == 1)#paths between dvs
  sel.dp  <- (MOD$code == 2)#paths from an iv
  sel.cov <- (MOD$code == 3)#vars and covs of ivs
  sel.var <- sel.cov & (MOD$to == MOD$from)#vars only
  
  #Start all NA path coefs at 1.0
  MOD$value[sel.ip] = 1.0
  MOD$value[sel.dp] = 1.0
  
  #Start all var/covs at 1.0 (update some below)
  MOD$value[sel.cov] = 1.0
  
  #Change NA covs to 0
  MOD$value[sel.cov & !sel.var] = 0.0
  
  #Update values for error vars
  listDV = unique(MOD$to[MOD$code != 3])
  sel.dvS = sort(listDV[listDV <= pOBS])
  MOD$value[sel.ev] = (0.9)*diag(S)[sel.dvS]
  # *** Error variances must be listed in order of Observed DVs!
  
  #Set start values for variance observed IVs
  listIV = unique(MOD$from[MOD$code != 3])
  sel.ivS = sort(listIV[listIV <= pOBS])
  sel.obsiv = sel.cov & sel.cov & (MOD$from <= pOBS)
  MOD$value[sel.obsiv] = diag(S)[sel.ivS]
  
  MOD
}

#######################################################################

#Model Covariance of Given free Parameters
SIGMAof <- function(theta) {
  #Combine New Free Values with Fixed Values
  pars = MOD$value
  pars[sel.free] = theta
  #Steps needed for iterations
  B[cbind(MOD$to[sel.dp],MOD$from[sel.dp])]=pars[sel.dp]
  Gam[cbind(MOD$to[sel.ip],ivfrom)] = pars[sel.ip]
  Gam[cbind(MOD$from[sel.ip],ivfrom)] = pars[sel.ip]
  Phi[Cindex] = pars[sel.cov]
  Phi = Phi + t(Phi) - diag(diag(Phi))
  #Get Sigma(pars)
  Binv = solve(diag(rep(1,p))-B)
  GBinvGam = G %*% Binv %*% Gam
  Sigma = GBinvGam %*% Phi %*% t(GBinvGam)
}

#######################################################################

#Vectorize Sigma
vechSIGMA <- function(theta) {vech(SIGMAof(theta))}

#Calculate standard errors
calcDelta <-function(theta,W) {
  #Tranform W for ML and GLS
  if(dim(W)[1] < pSTAR && dim(W)[1] > 1){
    Dp = D.matrix(pOBS)
    W = (t(Dp) %*% kronecker(W,W) %*% Dp)/2
  }
  #Get Derivatives
  Sdot = jacobian(vechSIGMA,t(theta))
  #Calculate Delta
  Delta = try(solve( t(Sdot) %*% W %*% Sdot),silent=TRUE)
}

#######################################################################

#Wrapper for calcDelta
#Checks for Convergence
getDelta <- function(rMETH,W) {
  #Check Convergence
  nlmconv = ((rMETH$oalg == 1) && (rMETH$conv <= 2)) #T if nlm converged
  optimconv = ((rMETH$oalg == 2) && (rMETH$conv == 0)) #T if optim converged
  #Get Delta if Convergence
  Delta = matrix(0,qFREE,qFREE)
  if(nlmconv || optimconv) {
    CheckDelta = calcDelta(rMETH$pars,W)
    if(typeof(CheckDelta) != "character") {
      Delta = CheckDelta
    }
  }
  Delta
}

#######################################################################

#Check for convergence and/or effectively singular covariance matrices
compsings <- function(rMETH) {
  #Check Convergence
  nlmconv = ((rMETH$oalg == 1) && (rMETH$conv <= 2)) #T if nlm converged
  optimconv = ((rMETH$oalg == 2) && (rMETH$conv == 0)) #T if optim converged
  SEtrouble = ( sum(!is.finite(rMETH$SEs) != 0) | sum( rMETH$SEs!=0) == 0 )
  if(nlmconv || optimconv) {
    if(SEtrouble == TRUE) {
      rMETH$conv = 6
    }
  }
  rMETH
}

#######################################################################

#Define Criterion Functions for various estimation procedures

#ADFcv criterion to be minimized
fADFcv <- function(theta) {
  #Get Sigma(pars)
  Sigma = SIGMAof(theta)
  #Minimize This
  diff = vech(PsyHat-Sigma)
  ADFval = t(diff) %*% SHPinv %*% (diff)
}

#ADFcvn criterion to be minimized
fADFcvn <- function(theta) {
  #Get Sigma(pars)
  Sigma = SIGMAof(theta)
  #Minimize This
  diff = vech(PsyHat-Sigma)
  ADFval = t(diff) %*% SHPNinv %*% (diff)
}

#GLScv criterion to be minimized 
fGLScv <- function(theta) {
  #Get Sigma(pars)
  Sigma = SIGMAof(theta)
  #Minimize This
  mid = (PsyHat-Sigma) %*% solve(PsyHat)
  GLSval = 0.5 * sum(diag( mid %*% mid )) 
}

#MLcv criterion to be minimized
fMLcv <- function(theta) {
  #Get Sigma(pars)
  Sigma = SIGMAof(theta)
  #Minimize This
  mid = (PsyHat-Sigma) %*% solve(Sigma)
  MLval = 0.5 * sum( diag( mid %*% mid ))
}

#######################################################################

#Apply minimization algorithms to criterion/cost function
plzcon <- function(theta,func) {
  result = new.env()
  #Try NLM: 
  resNLM = nlm(func,theta,iterlim = 1000)
  result$oalg = 1
  result$conv = resNLM$code
  result$mini = resNLM$minimum
  result$pars = resNLM$estimate
  #Check Code
  if(result$conv != 0) {
    resOPTIM = optim(theta,func,method='BFGS',control=list(maxit=1000))
    if(resOPTIM$value < result$mini) {
      result$oalg = 2
      result$conv = resOPTIM$convergence
      result$mini = resOPTIM$value
      result$pars = resOPTIM$par
    }
  }
  result$mini=N*(result$mini)
  result
}

#######################################################################

#Function to Print Results
printRES <- function(rMETH) {
  if(is.numeric(rMETH) == F) {
    #Convergence (0 good)
    print('Convergence (0 good)')
    print(rMETH$conv)
    #Optimization Method which produced minimum
    print('Optimization Method which produced minimum')
    print(rMETH$oalg)
    #Chi-squared Value
    print('Chi-squared Value and df')
    print(rMETH$mini)
    print(df)
    #Parameter Estimates
    print('Parameter Estimates (including fixed)')
    parlist = MOD$value
    parlist[sel.free] <- rMETH$pars
    print(parlist)
    #SEs of ADF Estimates
    print('SEs of Estimates')
    print(rMETH$SEs)
  }else {
    print("Still need to code printer for vectorized output.")
  }
}
