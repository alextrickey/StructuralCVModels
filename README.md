# StructuralCVModels
Structural Models of Coefficient of Variation Matrices

### Purpose
This repository contains the programs necessary to fit a structural equation 
model of coefficient of variation (CV) matrices. For more information, see Trickey 
(2015) referenced below. 

### Dependencies

The following R libraries are required: 
* `MASS` 
* `matrixcalc`
* `numDeriv` 

### Data File
A tab delimited data file should also be provided. See 'cv_sample_data.dat'
for an example. The data file should contain the values of the observed 
variables to be referenced in the model file. Do not include variables that 
will not be referenced in the model and preserve the variable ordering in the
model specification (see comments below). 

### Defining a Model (MOD) File
The model specification file should contain the necessary information to 
correctly specify a structural equation model and should represent all of the 
dependent variables in the associated data file. The file should be tab 
delimited, consisting of 5 columns (code, from, to, index, and value) with one
row for each path, variance, and covariance in the model (including error 
variances). 

The five columns should be specified as follows: 
* `code`
  + 1 if the row specifies a path from an independent 
  variable to a dependent variable
  + 2 if the row specifies a path from a dependent variable to another 
  dependent variable
  + 3 if the row specifies a variance, covariance, or error variance
* `from`
  + This is the number of the starting node in a given path. Note all observed
  variables must be **numbered starting with 1 in the order that they appear in 
  the data file**. 
  + For variances the start and end values should be the same. 
  + For covariances, the order of the start/end values is arbitrary. 
* `to` 
  + This is the number of the ending node in a given path. Note all observed
  variables must be numbered starting with 1 in the order that they appear in 
  the data file.
  + For variances the start and end values should be the same. 
  + For covariances, the order of the start/end values is arbitrary. 
* `index`
  + In order for the a structural equation model to be identified some variance
  or path coefficients will need to be fixed in order to ensure the model can be 
  estimated. The index value should be 0 if a particular row represents a fixed 
  coefficient, otherwise it should specify an number (from 1 to the number of 
  values to be estimated).
* `value`
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

### Fitting a Model with fitCVmod
Once you have created data and model specification files as described above. 
Place the files in the directory containing 'fitCVmod.r' and 'scvm_functions.r'.
Then modify the file names appropriately in fitCVmod.r (under the 'LOAD DATA 
& MODEL' heading). Then in the r console, run 'fitCVmod.r'. 
				  
### Results / Output
The program analyzes the input data using four difference estimation methods 
(see Trickey (2015) for details). The output for each estimation method 
includes the parameter estimates and standard errors (SEs) for any free 
variables (index numbers correspond with ordering specified in model file). The
output also displays the model chi-squared statistic for the test of model fit. 
(Recall in structural equation modeling low chi-squared values indicate better
fit.) The sample CV matrix, estimated CV matrix, and corresponding residuals are 
also provided. 

### For More Information
The reference below contains relevant theoretical background, performance tests 
and an application of the modeling techniques/algorithm provided here. 

### References
Trickey, K. A. (2015). Structural Models of Coefficient of Variation Matrices. 
UCLA: Psychology 0780. Retrieved from: http://escholarship.org/uc/item/0zj3b6qk
