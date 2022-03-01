rm(list = ls()) # Cleaning working environment

# Source useful functions
library(tidyverse) # For data manipulation
source("utils_formatting_functions.R") # Function for formatting

# Quanities

p <- 10; k <- 2; n_iterations <- 100


# Example when the ouput is a matrix (sigma2 for instance) ----------------

# Matrix where each column is a vector of estimate

example_sigma_2 = array(rnorm(p * n_iterations)^2, # Dummy values for example
                        dim = c(p, n_iterations)) # Dimensions

# Computing posterior mean

mean_sigma_2 <- rowMeans(example_sigma_2)
mean_sigma_2

# Plots

# First, we format to a data.frame (df)

sigma_2_df <- format_matrix(matrix_ = example_sigma_2, # Matrix to format
                            param_name = "sigma", # Name of the parameter
                            suffix_ = "^2") # Additional suffix (here the ^2, only useful here)
sigma_2_df # Look!
# We have three columns, 
#   -  the name of the parameter (p possible values here)
#   - The number of the iteration
#   - the estimate

# Representing the chain

ggplot(sigma_2_df) + # data plotted
  aes(x = iteration, y = Estimate) + # Which column on which axis
  facet_wrap(~Parameter, # facet_wrap gives one graph per group (here, Parameter)
             labeller = label_parsed, # Labels of graph parse to math expression 
             nrow = 2) + # 2 rows in the overall graph
  geom_line() # How do we represent? With a line

# Representing the densities

ggplot(sigma_2_df) + # data plotted
  aes(x = Estimate) + # Which column on which axis
  facet_wrap(~Parameter, # facet_wrap gives one graph per group (here, Parameter)
             labeller = label_parsed, # Labels of graph parse to math expression 
             scales = "free_y", # The y scale is not the same on each graph
             nrow = 2) + # 2 rows in the overall graph
  geom_density() # How do we represent? With a density




# Example when the ouput is an (lambda for instance) ----------------

example_lambda <- array(rnorm(p * k * n_iterations),
                       dim = c(p, k, n_iterations))

# Posterior mean, it is a matrix

apply(example_lambda, 
      MARGIN = c(1, 2), # The margins that are fixed (we fix the first two dimensions)
      mean) # Mean over the third dimension

lambda_df <- format_array(array_ = example_lambda, # Matrix to format
                            param_name = "Lambda") 
lambda_df  # Look!
# We have three columns, 
#   -  the name of the parameter (p times k possible values here)
#   - The number of the iteration
#   - the estimate

# Representing the chain

ggplot(lambda_df) + # data plotted
  aes(x = iteration, y = Estimate) + # Which column on which axis
  facet_wrap(~Parameter, # facet_wrap gives one graph per group (here, Parameter)
             labeller = label_parsed, # Labels of graph parse to math expression 
             nrow = 10) + # 2 rows in the overall graph
  geom_line() # How do we represent? With a line

# Representing the densities

ggplot(lambda_df) + # data plotted
  aes(x = Estimate) + # Which column on which axis
  facet_wrap(~Parameter, # facet_wrap gives one graph per group (here, Parameter)
             labeller = label_parsed, # Labels of graph parse to math expression 
             scales = "free_y", # The y scale is not the same on each graph
             nrow = 10) + # 2 rows in the overall graph
  geom_density() # How do we represent? With a density