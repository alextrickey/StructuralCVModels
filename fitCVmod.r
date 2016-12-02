# Analyze Data

# Requires loading/entering/computing
#   MOD 
#     -the model specification matrix 
#     -based on RAM specification in "sem" package
#     -create mod from path diagram (see below)
#   Data
#     -the covariance matrix of the observed variables
#     -variable order should be consistent with MOD

###########################
##Load Required Packages ##
###########################

require('MASS') #Multivariate Normal Distributions
require('matrixcalc') #Special Matrices
require('numDeriv') #Needed to estimate SEs
source('scvm_functions.r') #Custom functions to fit CV models


#######################
## LOAD DATA & MODEL ##
#######################

#Enter data file name here:
Y = read.table('cv_sample_data.dat',header = TRUE)
head(Y)

#Specify MOD file here: 
MOD <- read.delim("MOD.dat") ## Read in model here


###################################
## Count and Selection Variables ##
###################################

#sample size
N = nrow(Y)

#total number of variables 
p = max(c(MOD$from,MOD$to))
cat(paste('Total Number of Variables: ',p))
# ***Assumes no unused variables

#number of observed variables
pOBS = ncol(Y) 
cat(paste('Number of Observed Variables: ', pOBS))
# ***Assumes variables 1:pOBS are all observed and included in model

#Number of unique covariances
pSTAR = pOBS*(pOBS+1)/2
cat(paste('Number of Unique Covariance Terms: ', pSTAR))

#Select/Identify Paths and Vars
sel.ip  <- (MOD$code == 1)#paths from ivs
sel.dp  <- (MOD$code == 2)#paths between dvs
sel.cov <- (MOD$code == 3)#vars and covs of ivs

#Identify and count (strict) IVs 
listIV = unique(MOD$from[sel.ip])
pIV = length(listIV) #includes errors
cat(paste('Number of Strictly Independent Variables: ', pIV))

#Identify and count (all) DVs
listDV = unique(MOD$to[MOD$code != 3])
pDV = length(listDV)
cat(paste('Number of Dependent Variables: ', pDV))

#Select/Identify Free Parameters
sel.fixed <- (MOD$index == 0)
sel.free <- !(sel.fixed)

#number of parameters
qTOT = nrow(MOD)
qFREE = sum(sel.free)
cat(paste('Number of Free Parameter: ', qFREE))
# *** Assumes no constraints on parameters

#df for chi-squared
df = pSTAR - qFREE
cat(paste('Number of Degrees of Freedom: ', df))

############################
## Some Useful Statistics ##
############################

#The mean
M = colMeans(Y)
#Warning: This doesn't return a transposable vector! 
DmInv = diag(1/M)
cat('Sample Means:')
print(M)

#Sample Covariance Matrix
S = var(Y)
cat('Sample Covariance Matrix:')
print(S)

#Sample CV Matrix
PsyHat = DmInv %*% S %*% DmInv
cat('Sample Coefficient of Variation (CV) Matrix:')
print(PsyHat)

##########################
## Compute Start Values ##
##########################
MOD = getstart(MOD,S)

#############################################
## Define Matrices for Bentler-Weeks Model ##
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

METH = 'ADF Estimation with Arbitrary Distribution'
#Fit Model Using ADFcv procedure
  theta = MOD$value[sel.free]
  rADFcv = plzcon(theta,fADFcv)
  rADFcv$METH = METH
#SEs of ADFcv Estimates
  DeltaADFcv = getDelta(rADFcv,SHPinv)
  rADFcv$SEs = sqrt(diag(DeltaADFcv)/N)
  rADFcv=compsings(rADFcv)
#Results
  printRES(rADFcv)
  cat('Sample CV Matrix:')
  print(PsyHat)
  cat('Reproduced CV Mat:')
  print(SIGMAof(rADFcv$pars))
  cat('Covariance Residuals:')
  sdres_ADF = (vech(PsyHat) - vech(SIGMAof(rADFcv$pars)))/sd(vech(PsyHat))
  print(PsyHat - SIGMAof(rADFcv$pars))
  cat('Covariance Residuals as Percent Difference:')
  print((PsyHat - SIGMAof(rADFcv$pars))*100/PsyHat)


#####################
## ADF - CV Normal ##
#####################

METH = 'ADF Estimation with Normal Distribution'
#Fit Model Using ADFcvn procedure
  theta = MOD$value[sel.free]
  rADFcvn = plzcon(theta,fADFcvn)
  rADFcvn$METH = METH
#SEs of ADF Estimates
  DeltaADFcvn = getDelta(rADFcvn,SHPNinv)
  rADFcvn$SEs = sqrt(diag(DeltaADFcvn)/N) 
  rADFcvn=compsings(rADFcvn)
#Results
  printRES(rADFcvn)
  cat('Sample CV Matrix:')
  print(PsyHat)
  cat('Reproduced CV Mat:')
  print(SIGMAof(rADFcvn$pars))
  cat('Covariance Residuals:')
  sdres_ADFN = (vech(PsyHat) - vech(SIGMAof(rADFcvn$pars)))/sd(vech(PsyHat))
  print(PsyHat - SIGMAof(rADFcvn$pars))
  cat('Covariance Residuals as Percent Difference:')
  print((PsyHat - SIGMAof(rADFcvn$pars))*100/PsyHat)

##############
## GLS - CV ##
##############

METH = 'GLS Estimation'
#Fit Model Using GLS procedure
  theta = MOD$value[sel.free]
  rGLScv = plzcon(theta,fGLScv)
  rGLScv$METH = METH
#SEs of GLS Estimates
  DeltaGLScv = getDelta(rGLScv,solve(PsyHat))
  rGLScv$SEs = sqrt(diag(DeltaGLScv)/N)
  rGLScv=compsings(rGLScv)
#Results
  printRES(rGLScv)
  cat('Sample CV Mat:')
  print(PsyHat)
  cat('Reproduced CV Mat:')
  print(SIGMAof(rGLScv$pars))
  cat('Covariance Residuals:')
  sdres_GLS = (vech(PsyHat) - vech(SIGMAof(rGLScv$pars)))/sd(vech(PsyHat))
  print(PsyHat - SIGMAof(rGLScv$pars))
  cat('Covariance Residuals as Percent Difference:')
  print((PsyHat - SIGMAof(rGLScv$pars))*100/PsyHat)

#############
## ML - CV ##
#############

METH = 'Maximum Likelihood Estimation'
#Fit Model Using ML procedure
  theta = MOD$value[sel.free]
  rMLcv = plzcon(theta,fMLcv)
  rMLcv$METH = METH
#SEs of ML Estimates
  DeltaMLcv = getDelta(rMLcv,solve(SIGMAof(rMLcv$pars)))
  rMLcv$SEs = sqrt(diag(DeltaMLcv)/N)
  rMLcv=compsings(rMLcv)
#Result
  printRES(rMLcv)
  cat('Sample CV Mat:')
  print(PsyHat)
  cat('Reproduced CV Mat:')
  print(SIGMAof(rMLcv$pars))
  cat('Covariance Residuals:')
  sdres_ML = (vech(PsyHat) - vech(SIGMAof(rMLcv$pars)))/sd(vech(PsyHat))
  print(PsyHat - SIGMAof(rMLcv$pars))
  cat('Covariance Residuals as Percent Difference:')
  print((PsyHat - SIGMAof(rMLcv$pars))*100/PsyHat)

#Plot Histograms of Standardized Residuals
par(mfrow = c(2,2),
    mar=c(3, 4, 1.5, 1.5)+0.1,
    cex=0.9)
xtitle = "Standardized Residual"
hist(sdres_ADF,
     ylab=NULL,xlab=NULL,main=NULL,
     ylim=c(0,8),xlim=c(-0.5,0.5),
     col="gray")
mtext("AGLS",3,line=-1.5,cex=0.9)
mtext(xtitle,1,line=2,cex=0.9)
mtext("Frequency",2,line=2,cex=0.9)

hist(sdres_ADFN,
     ylab=NULL,xlab=NULL,main=NULL,
     ylim=c(0,8),xlim=c(-0.5,0.5),
     col="gray")
mtext("NGLS",3,line=-1.5,cex=0.9)
mtext(xtitle,1,line=2,cex=0.9)
mtext("Frequency",2,line=2,cex=0.9)

hist(sdres_GLS,
     ylab=NULL,xlab=NULL,main=NULL,
     ylim=c(0,8),xlim=c(-0.5,0.5),
     col="gray")
mtext("MGLS",3,line=-1.5,cex=0.9)
mtext(xtitle,1,line=2,cex=0.9)
mtext("Frequency",2,line=2,cex=0.9)

hist(sdres_ML,
     ylab=NULL,xlab=NULL,main=NULL,
     ylim=c(0,8),xlim=c(-0.5,0.5),
     col="gray")
mtext("MRLS",3,line=-1.5,cex=0.9)
mtext(xtitle,1,line=2,cex=0.9)
mtext("Frequency",2,line=2,cex=0.9)

Parameter_Table = rbind(rADFcv$pars, rADFcv$SEs,
                  rADFcvn$pars,rADFcvn$SEs,
                  rGLScv$pars, rGLScv$SEs,
                  rMLcv$pars,  rMLcv$SEs)
View(Parameter_Table)

CV_Matrix_Table = cbind(vech(PsyHat), 
                  vech(SIGMAof(rADFcv$pars)), 
                  vech(SIGMAof(rADFcvn$pars)), 
                  vech(SIGMAof(rGLScv$pars)),
                  vech(SIGMAof(rMLcv$pars)))
View(CV_Matrix_Table)
