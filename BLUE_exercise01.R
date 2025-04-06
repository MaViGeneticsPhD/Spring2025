library(MASS)  # Load the MASS package for matrix operations

# Define the design matrix (X)
# This matrix represents different levels of categorical variables
X = matrix(data=c(1,0,0,  # First group
                  1,0,0,  # First group
                  1,0,0,  # First group
                  1,1,0,  # Second group
                  1,1,0,  # Second group
                  1,1,0,  # Second group
                  1,0,1,  # Third group
                  1,0,1,  # Third group
                  1,0,1), # Third group
           byrow = TRUE,  # Fill matrix by rows
           nrow=9,        # 9 observations
           ncol=3)        # 3 predictors (including intercept)

# Compute transpose of X (Xt)
Xt = t(X)  # Transpose of X

# Compute X'X (XtX), which is a key component in regression
XtX = Xt %*% X  # Matrix multiplication

# Display X'X matrix
XtX

# Define response variable matrix (Y)
# These are the observed values for the dependent variable
Y = matrix(data=c(28,30,32,  # First group
                  35,36,34,  # Second group
                  40,42,39), # Third group
           byrow = TRUE,     # Fill matrix by rows
           nrow=9,          # 9 observations
           ncol=1)          # 1 dependent variable

# Compute X'Y (XtY), another key regression component
XtY = Xt %*% Y  # Matrix multiplication

# Display X'Y matrix
XtY

# Compute the inverse of X'X (invXtX)
invXtX = solve(XtX)  # Solve for the inverse of XtX

# Display the inverse of X'X
invXtX

# Compute the regression coefficient estimates (B)
B = invXtX %*% XtY  # Solve the normal equations for B

# Display estimated regression coefficients
B
