## StructuralCVModels
Structural Models of Coefficient of Variation Matrices



#### Dependencies

The following R libraries are required: 
'MASS' #Multivariate Normal
'matrixcalc' #Special Matrices
'numDeriv' #Needed for SEs

#### Defining a Model (MOD) file

The model specification file should contain the necessary information to 
correctly specify a structural equation model and should represent all of the 
dependent variables in the associated data file. The file should be tab 
delimited, consisting of 5 columns (code, from, to, index, and value) with one
row for each path, variance, and covariance in the model (including error 
variances). 

The five columns should be specified as follows: 
* code
  + 1 if the row specifies a path from an independent 
  variable to a dependent variable
  + 2 if the row specifies a path from a dependent variable to another 
  dependent variable
  + 3 if the row specifies a variance, covariance, or error variance
* from
  + This is the number of the starting node in a given path. Note all observed
  variables must be numbered starting with 1 in the order that they appear in 
  the data file. 
  + For variances the start and end values should be the same. 
  + For covariances, the order of the start/end values is arbitrary. 
* to 
  + This is the number of the ending node in a given path. Note all observed
  variables must be numbered starting with 1 in the order that they appear in 
  the data file.
  + For variances the start and end values should be the same. 
  + For covariances, the order of the start/end values is arbitrary. 
* index
  + In order for the a structural equation model to be identified some variance
  or path coefficients will need to be fixed in order to ensure the model can be 
  estimated. The index value should be 0 if a particular row represents a fixed 
  coefficient, otherwise it should specify an number (from 1 to the number of 
  values to be estimated).
* value
  + The value column should contain the value for any fixed variables or any 
  desired coefficient starting values for to-be-estimated parameters. If 
  reasonable
  starting values are not known, enter E for unknown error variances and X for 
  any other unknow coefficient (including path, variance, and covariance 
  coefficients). 

See the MOD.dat file for an example of a model spefication file for a two-factor
model with 3 observed variables per factor. The observed variables are numbered 
1 through 6 and the factors are numbered 7 and 8. The remaining terms are errors
for the dependent variables. There is one row in the file for each path, variance
and covariance in the model. 

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


#### For More Information
The reference below contains relevant theoretical background, simulations and an application of the modeling techniques/algorithm provided here. 

Trickey, K. A. (2015). Structural Models of Coefficient of Variation Matrices. UCLA: Psychology 0780. Retrieved from: http://escholarship.org/uc/item/0zj3b6qk

