# Animal Model

# Clean the working environment
rm(list = ls())  # Removes all objects from the workspace to ensure a clean environment

# Load required packages
library("pedigreemm")  # Provides functions for pedigree-based mixed models
library("tidyverse")   # A collection of data manipulation and visualization packages
# install.packages("pedigreemm")  # Uncomment to install if not already installed
# install.packages("tidyverse")  # Uncomment to install if not already installed

# Example 4.3 -------------------------------------------------------------

# Prepare dataset
wwg = c(4.5, 2.9, 3.9, 3.5, 5.0)  # Weaning weight gain values
sex = c("Male", "Female", "Female", "Male", "Male")  # Sex of individuals
an = seq(4, 8)  # Animal ID numbers
sn = c(1, 3, 1, 4, 3)  # Sire ID (father)
dn = c(NA, 2, 2, 5, 6)  # Dam ID (mother), NA if unknown

data = data.frame(an, sn, dn, sex, wwg)  # Create a data frame

data$sex = as.factor(sex)  # Convert sex to a factor variable
data$sex = relevel(factor(sex), ref = "Male")  # Set "Male" as the reference level for the factor
data$an = factor(x = data$an, levels = seq(1, 8))  # Convert animal ID to a factor with explicit levels

# Prepare pedigree data
a = seq(1, 6)  # Define pedigree ID numbers
s = c(NA, NA, NA, 1, 3, 1)  # Define sires (fathers)
d = c(NA, NA, NA, NA, 2, 2)  # Define dams (mothers)
ped = data.frame(a, s, d)  # Create pedigree data frame

# Create pedigree object
pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

# Compute the inverse of the additive relationship matrix (A-matrix)
Ainv = getAInv(ped = pedX)

# Define variance components
varA = 20  # Additive genetic variance
varE = 40  # Residual variance

# Compute variance ratio
alpha = varE / varA  # Ratio of environmental to genetic variance

# Setting up the incidence matrices for the Mixed Model Equations (MME)
X = model.matrix(wwg ~ -1 + sex, data = data)  # Design matrix for fixed effects (sex)

# Initialize the incidence matrix for random genetic effects (W)
W = matrix(0, nrow = nrow(data), ncol = nrow(Ainv))  # Matrix with zeros
rownames(W) = as.character(data$an)  # Assign row names as animal IDs
colnames(W) = as.character(pedX@label)  # Assign column names as pedigree IDs

# Populate the W matrix with identity values for direct genetic effects
for (i in 1:nrow(data)) {
  a_i = as.character(data[i, "an"])
  if (a_i %in% colnames(W)) {
    W[i, a_i] = 1  # Direct genetic effect assigned
  } else {
    s_i = as.character(data[i, "sn"])
    d_i = as.character(data[i, "dn"])
    W[i, s_i] = 0.5  # Half contribution from sire
    W[i, d_i] = 0.5  # Half contribution from dam
  }
}

# Initialize the residual covariance matrix (R)
R = matrix(0, nrow = nrow(data), ncol = nrow(data))
rownames(R) = as.character(data$an)
colnames(R) = as.character(data$an)

# Populate the residual covariance matrix
for (i in 1:nrow(data)) {
  a_i = as.character(data[i, "an"])
  if (a_i %in% unique(c(data$sn, data$dn))) {
    R[i, a_i] = varE  # Assign residual variance for offspring
  } else {
    R[i, a_i] = varE + varA * 0.5  # Adjusted variance for direct genetic effect
  }
}

# Compute the inverse of the residual covariance matrix
Ri = solve(R)

# Compute components of the Mixed Model Equations (MME)
XpRiX = crossprod(X, Ri) %*% X  # X'R⁻¹X
XpRiW = crossprod(X, Ri) %*% W  # X'R⁻¹W
WpRiX = crossprod(W, Ri) %*% X  # W'R⁻¹X
WpRiW = crossprod(W, Ri) %*% W  # W'R⁻¹W

# Extract response variable (weaning weight gain)
y = na.omit(data$wwg) 
XpRiy = crossprod(X, Ri) %*% y  # X'R⁻¹y
WpRiy = crossprod(W, Ri) %*% y  # W'R⁻¹y

# Construct the Left-Hand Side (LHS) matrix of MME
LHS = rbind(
  cbind(XpRiX, XpRiW), 
  cbind(WpRiX, WpRiW + Ainv * (1/varA))  # Incorporate A-matrix inverse
)

# Construct the Right-Hand Side (RHS) vector of MME
RHS = rbind(
  XpRiy, 
  WpRiy
)

# Solve for the fixed and random effects
solutions = solve(LHS, RHS)

# Display rounded solutions
round(solutions, 3)

