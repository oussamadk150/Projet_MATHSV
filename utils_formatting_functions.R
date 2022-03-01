# Function to format an array of estimated parameters to a data.frame
# easily handle by ggplot2

format_array <- function(array_, param_name, row_indexes = NULL){
  if(is.null(row_indexes)){
    row_indexes <- 1:(dim(array_)[1])
  }
  else{
    row_indexes <- sort(row_indexes)
    array_ <- array_[row_indexes,,]
  }
  dim1 <- length(row_indexes)
  dim2 <- dim(array_)[2]
  n_iteration <- dim(array_)[3]
  names_levels <- paste0(param_name, "[",
                         rep(row_indexes, each = dim2),
                         "-",
                         rep(1:dim2, dim1),
                         "]")
  map_dfr(1:n_iteration,
          function(l){
            array_[,, l] %>% 
              as.data.frame() %>% 
              setNames(paste0("col", 1:dim2)) %>% 
              mutate(iteration = l,
                     row = row_indexes) %>% 
              pivot_longer(-c("iteration", "row"), 
                           names_to = "column",
                           names_prefix = "col",
                           values_to = "Estimate") %>% 
              unite(Parameter, row, column, sep = '-') %>% 
              mutate(Parameter = paste0(param_name, "[", Parameter, "]"),
                     Parameter = factor(Parameter,
                                        levels = names_levels)) %>% 
              arrange(Parameter)
          })
}

# Function to format a matrix of estimated parameters to a data.frame
# easily handle by ggplot2
format_matrix <- function(matrix_, # Matrix of values d x n_iterations
                          param_name, # String giving parameter name
                          suffix_ = NULL){ # suffix of the name typically "^2" when
  # the name is sigma
  dim1 <- dim(matrix_)[1] # DImension of the parameter
  n_iteration <- dim(matrix_)[2]
  names_levels <- paste0(param_name, "[", 1:dim1, "]", suffix_)
  matrix_ %>% 
    as.data.frame() %>% 
    setNames(paste0("it", 1:n_iteration)) %>%
    mutate(Parameter = names_levels) %>%  
    pivot_longer(-Parameter,
                 names_to = "iteration",
                 names_prefix = "it", 
                 names_transform = list(iteration = as.numeric),
                 values_to = "Estimate") %>% 
    mutate(Parameter = factor(Parameter, level = names_levels)) %>% 
    arrange(Parameter)
}