# Generate Sample Data for Test Runs
#   MOD 
#     -the model specification matrix 
#     -based on RAM specification in "sem" package
#     -create mod from path diagram (see below)
#   S
#     -the covariance matrix of the observed variables
#     -variable order should be consistent with MOD

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
SIGMA = diag(MU) %*% PSY %*% diag(MU)


###########################
##Load Required Packages ##
###########################

require('MASS') #Multivariate Normal
require('matrixcalc') #Special Matrices
#require('bindata') #will need eventually
require('numDeriv') #Needed for SEs

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

##Seed File Name (Store Random Seed)
if(getnewseed == 'yes'){  
  seedfile = paste(savefolder,"Seed",format(Sys.time(), "%b%d.%H%M"),".Rdata",sep="")
}

#Function to Store Data (Must Math Col Names Below)
storeRES <- function(rMETH,filename="") {
  cat(c(N,REP,rMETH$METH,rMETH$conv,rMETH$oalg,rMETH$mini,rMETH$pars,rMETH$SEs),file=filename,append=T)
  cat("\n",file=filename,append=T)
}

#Store Random Seed
if(getnewseed == 'yes'){
  save_rng(savefile = seedfile)
}


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
