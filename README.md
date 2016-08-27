## StructuralCVModels
Structural Models of Coefficient of Variation Matrices



#### Dependencies

The following R libraries are required: 
'MASS' #Multivariate Normal
'matrixcalc' #Special Matrices
'numDeriv' #Needed for SEs
** Currently uses R-studio View() [remove this]

#### Defining a Model (MOD) file


#### Fitting a Model with fitCVmod


#### Using Test Model and Sample Data
Generate Sample Data for Test Runs
MOD 
- the model specification matrix 
- based on RAM specification in "sem" package
- create mod from path diagram (see below)
S
- the covariance matrix of the observed variables
- variable order should be consistent with MOD

Y_try = make_data(use_random_seed=TRUE,
				  existing_seed_file = "RandomSeed_2016Aug09_2301.Rdata",
				  save_data = TRUE)
				  
#### Reading Results 


#### For More Info
add ref