# Simulation to Test Standard Errors of Parameter Estimates
#   -11.24.14 Updated Computational Singularity Error Checking
#   -11.25.14 Change GLS W matrix in Delta formula (affects SEs only)
#   -11.25.14 Removed N/(N-pSTAR-1) correction from ADF methods

# Requires loading/entering/computing
#   MOD 
#     -the model specification matrix 
#     -based on RAM specification in "sem" package
#     -create mod from path diagram (see below)
#   S
#     -the covariance matrix of the observed variables
#     -variable order should be consistent with MOD

######################
## Simulation Setup ##
######################

#Save Files in This Folder:
savefolder = "~/Dropbox/Dissertation/SimResults/F2_normal/F2p6/" #Linux
#savefolder = "c:/Users/Trickey/Dropbox/Dissertation/SimResults/F2_normal/F2p6/" #Linux

#Save a new seed for this sim? 
getnewseed = 'yes' #'yes' or 'no' (if no specify seed file on next line)
#seedfile = paste(savefolder,"SeedNov04.Rdata",sep="")
#folder must contain MOD file, seed file, and SimCV program, 


#Load MOD file (Must be copied to save folder)
MOD <- read.delim(paste(savefolder,"MOD.dat",sep=""))

#Define Population Matrices (PSYID = 1?)
pF1  = 3 #pF1=pF2
pOBS = 6
Fcov = 0.3 #Factor Covariance (yields cor(F1,F2)=0.2)
OffBlocks = matrix(Fcov,pF1,pF1)
DiagBlocks = matrix(1,pF1,pF1) + diag(pF1)
PSY = rbind(cbind(DiagBlocks,OffBlocks),cbind(OffBlocks,DiagBlocks))


MU = rep(1,pOBS)
#MU = seq(1,3,length.out=pOBS)
SIGMA = diag(MU) %*% PSY %*% diag(MU)

#sample sizes
Nvals=c(100,300,500,1000,10000,100000)

#number of replications
NREP = 500

###########################
##Load Required Packages ##
###########################

require('MASS') #Multivariate Normal
require('matrixcalc') #Special Matrices
#require('bindata') #will need eventually
require('numDeriv') #Needed for SEs

#####################
##Define Functions ##
#####################
#####################
# List:
#   save_rgn, restore_rgn: get and restore current seed
#   getstart(MOD,S): Calculates Starting Values
#   SIGMAof(theta): Calculated Covariance of Free Parameters
#   vechSIGMA(theta): Returns half-vectorize of SIGMAof(theta)
#   calcDelta(theta,W): Calculates "Delta" for standard errors
#   getDelta(theta,W): Error Checking Wrapper for calcDelta
#   Criterion Functions: Return crit to be minimized given theta
#     fADFcv, fADFcvn, fGLScv, fMLcv
#   plzcon(theta,func): Optimizes Criterion Func to Estimate Pars
#   printRES: Prints results in form set up by plzcon and SEs
#   storeRES: Prints results to file?

#functions for handling random seeds
  #saves current seed
    save_rng <- function(savefile=tempfile()) {
  if (exists(".Random.seed"))  {
    oldseed <- get(".Random.seed", .GlobalEnv)
  } else stop("don't know how to save before set.seed() or r*** call")
  oldRNGkind <- RNGkind()
  save("oldseed","oldRNGkind",file=savefile)
  invisible(savefile)
}
  #restores saved seed
    restore_rng <- function(savefile) {
  load(savefile)
  do.call("RNGkind",as.list(oldRNGkind))  ## must be first!
  assign(".Random.seed", oldseed, .GlobalEnv)
}

#Computes Starting Values for MOD
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

#Check for Computational Singularities
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

#Optimization Methods
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


###################################
## Count and Selection Variables ##
###################################
###################################

#total number of variables 
p = max(c(MOD$from,MOD$to))
print(p)
# ***Assumes no unused variables

#number of observed variables
pOBS = nrow(SIGMA) 
print(pOBS)
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

###########################
## Prepare to Start Loop ##
###########################

##Seed File Name (Store Random Seed)
if(getnewseed == 'yes'){  
  seedfile = paste(savefolder,"Seed",format(Sys.time(), "%b%d.%H%M"),".Rdata",sep="")
}

#Function to Store Data (Must Math Col Names Below)
storeRES <- function(rMETH,filename="") {
  cat(c(N,REP,rMETH$METH,rMETH$conv,rMETH$oalg,rMETH$mini,rMETH$pars,rMETH$SEs),file=filename,append=T)
  cat("\n",file=filename,append=T)
}

#Set Up File and Column Names (Must Match Result List Above)
ofile = paste(savefolder,"Res",format(Sys.time(), "%b%d"),".dat",sep="")
cat("N","REP","METH","conv","oalg","mini",paste("par",1:qFREE,sep=""),paste("SE",1:qFREE,sep=""),file=ofile,append=T)
cat("\n",file=ofile,append=T)

#################
## Start Loops ##
#################

#Store Random Seed
if(getnewseed == 'yes'){
  save_rng(savefile = seedfile)
}

#Status:
cat("Starting Iterations: ", date())

#Loop Over Sample Sizes
for(N in Nvals){
cat("\n")
cat("Beginning",NREP, "replications for N =", N, "\n")
cat("Replication: ")

#Restore Seed
restore_rng(seedfile)
  
for(REP in 1:NREP) {
cat(REP," ")

#################
## Draw Sample ##
#################
Y=mvrnorm(N,MU,SIGMA)

############################
## Some Useful Statistics ##
############################
############################

#The mean
M = colMeans(Y)
#Warning: This doesn't return a transposable vector! 
DmInv = diag(1/M)
M

#Sample Covariance Matrix
S = var(Y)

#Sample CV Matrix
PsyHat = DmInv %*% S %*% DmInv

##########################
## Compute Start Values ##
##########################
MOD = getstart(MOD,S)

#############################################
## Define Matrices for Bentler-Weeks Model ##
#############################################
#############################################

#Get Matrix G to Select Observed Variables
  G = matrix(0,pOBS,p)
  G[cbind(1:pOBS,1:pOBS)] <- 1
  # ***Assumes variables 1:pOBS are observed and included

#Get Matrix B to relate DVs with DVs (contains "beta and known 0s")
  B = matrix(0,p,p)
  B[cbind(MOD$to[sel.dp],MOD$from[sel.dp])]=MOD$value[sel.dp]

#Define Reindex Function
  #function to recode IVs for Gam, Phi
  reindex <- function(x) {
    #initialize
    xnew = x
    #get ordered list of unique xvals
    sux = sort(unique(x))
    #loop through x
    for(i in 1:length(x)) {
      #loop over unique values
      for(k in 1:length(sux)) {
        if(x[i] == sux[k])
          xnew[i] = k
      }
    }
    xnew
  }

#Get Matrix Gamma to relate IVs with DVs
  #Initialize
    Gam = matrix(0,p,pIV)
  #List IVs with New Unique Indices
    ivfrom = reindex(MOD$from[sel.ip])
  #Calculate Gamma
    Gam[cbind(MOD$to[sel.ip],ivfrom)] = MOD$value[sel.ip]
    Gam[cbind(MOD$from[sel.ip],ivfrom)] = MOD$value[sel.ip]
  
#Get CovMat Phi containing covs of IVs
  #Initialize
    Phi = diag(1,pIV,pIV)
  #Find Covariances
    these = c(unique(MOD$from[sel.ip]),MOD$to[sel.cov],MOD$from[sel.cov])
  #List Covariances with New Unique Indices
    Cindex = reindex(these)
    Cindex = Cindex[-(1:pIV)]
    Cindex = matrix(Cindex,ncol=2,byrow=F)
  #Put Covariances into Phi Matrix
    Phi[Cindex] = MOD$value[sel.cov]
    Phi = Phi + t(Phi) - diag(diag(Phi))

###########################################
## Get Sampling Variances for ADF Method ##
###########################################

#################
## SigmaHatPsy ##
#################
#################

#Get SigmaPsyHat for CV-SEM Method
  #Identity
    Ip = diag(pOBS)
  #Kronecker Product of CV and Ip
    PsyHatI = kronecker(PsyHat,Ip)
  #Duplication and Elimination Matrices
    Dp = D.matrix(pOBS)
    Hp = solve(t(Dp) %*% Dp) %*% t(Dp)
  #Np (Funky Product of Dup Matrices)
    Np = Dp %*% solve( t(Dp) %*% Dp ) %*% t(Dp)
  #Lp
    Lp=matrix(0,pOBS^2,pOBS)
    for(i in 1:pOBS) {
      Lp = Lp + kronecker(Ip[,i],Ip[,i]) %*% t(Ip[,i])
    }
  #Transform Means Vector
    kronDmInv = kronecker(DmInv,DmInv)
  #Estimate Omegas [Prospectus (3.11) and (3.12)]
    #Get Sums of Products for OmegaHat12 and OmegaHat22
    sumthing12 = matrix(0,pOBS,pOBS^2)
    sumthing22 = matrix(0,pOBS^2,pOBS^2)
      for(i in 1:N) {
        #Precursors for Omega12
        ui = kronecker(Y[i,]-M,Y[i,] - M)
        prodthing12 = (Y[i,] - M) %*% t(ui)
        sumthing12 = sumthing12 + prodthing12
        #Precursors for Omega22
        prodthing22 = ui %*% t(ui)
        sumthing22 = sumthing22 + prodthing22
      }
    #Get OmegaHat12
      OmegaHat12 = N/((N-1)*(N-2))* sumthing12
      OmegaHat12str =DmInv %*% OmegaHat12 %*% kronDmInv
    #Get OmegaHat22
      OmegaHat22 = (sumthing22 - (N-2)*vec(S)%*%t(vec(S)) )/(N-pSTAR-1)
      OmegaHat22str = kronDmInv %*% OmegaHat22 %*% kronDmInv
  #Variance of CV Sampling Distribution [Prospectus (3.15)]
    SigmaHatPsy = (2*Np) %*% PsyHatI %*% Lp %*% PsyHat %*% t(Lp) %*% 
      PsyHatI %*% (2*Np) +
      -t(OmegaHat12str) %*% t(Lp) %*% PsyHatI %*% (2*Np) +
      -(2*Np) %*% PsyHatI %*% Lp %*% OmegaHat12str + OmegaHat22str
    SHPinv = solve(Hp %*% SigmaHatPsy %*% t(Hp))
  #Variance of CV Sampling Distribution w/Normality [Prospectus (3.16)]
    SigmaHatPsyN = 2*Np %*% PsyHatI %*% Lp %*% PsyHat %*% t(Lp) %*% 
      PsyHatI %*% (2*Np) +
      2*Np %*% kronecker(PsyHat,PsyHat)
    SHPNinv = solve(Hp %*% SigmaHatPsyN %*% t(Hp))


#########################################
## Fit CV Models and Examine Estimates ##
#########################################

#########################
## ADF - CV (Arbitrary)##
#########################
#########################
METH = 4
#Fit Model Using ADFcv procedure
  theta = MOD$value[sel.free]
  rADFcv = plzcon(theta,fADFcv)
  rADFcv$METH = METH
#SEs of ADFcv Estimates
  DeltaADFcv = getDelta(rADFcv,SHPinv)
  rADFcv$SEs = sqrt(diag(DeltaADFcv)/N)
  rADFcv=compsings(rADFcv)
#Results
  storeRES(rADFcv,filename=ofile)

#####################
## ADF - CV Normal ##
#####################
#####################
METH = 5
#Fit Model Using ADFcvn procedure
  theta = MOD$value[sel.free]
  rADFcvn = plzcon(theta,fADFcvn)
  rADFcvn$METH = METH
#SEs of ADF Estimates
  DeltaADFcvn = getDelta(rADFcvn,SHPNinv)
  rADFcvn$SEs = sqrt(diag(DeltaADFcvn)/N) 
  rADFcvn=compsings(rADFcvn)
#Results
  storeRES(rADFcvn,filename=ofile)

##############
## GLS - CV ##
##############
##############
METH = 6
#Fit Model Using GLS procedure
  theta = MOD$value[sel.free]
  rGLScv = plzcon(theta,fGLScv)
  rGLScv$METH = METH
#SEs of GLS Estimates
  DeltaGLScv = getDelta(rGLScv,solve(PsyHat))
  rGLScv$SEs = sqrt(diag(DeltaGLScv)/N)
  rGLScv=compsings(rGLScv)
#Results
  storeRES(rGLScv,filename=ofile)

#############
## ML - CV ##
#############
#############
METH = 7
#Fit Model Using ML procedure
  theta = MOD$value[sel.free]
  rMLcv = plzcon(theta,fMLcv)
  rMLcv$METH = METH
#SEs of ML Estimates
  DeltaMLcv = getDelta(rMLcv,solve(SIGMAof(rMLcv$pars)))
  rMLcv$SEs = sqrt(diag(DeltaMLcv)/N)
  rMLcv=compsings(rMLcv)
#Result
  storeRES(rMLcv,filename=ofile)

###############
## End Loops ##
###############

} #End Loop Over Replications
} #End Loop Over Sample Sizes

cat("End time: ", date())
