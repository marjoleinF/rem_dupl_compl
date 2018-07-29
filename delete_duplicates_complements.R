##########
##
## Internal function to remove duplicate and complement rules:
##
delete_duplicates_complements <- function(
  rules, data, removecomplements = TRUE, removeduplicates = TRUE,
  return.dupl.compl = FALSE, sparse = FALSE) {
  ## Generate rule variables:
  if(sparse){
    # See https://stackoverflow.com/a/8844057/5861244. 
    #
    # if all rules where binary then we could use the `lsparseMatrix-classes`
    # this would though require that we either check that they are all binary 
    # or we pass an argument which states that this is the case
    expr <- paste0("cbind_sparse_vec(", paste0(
      'as(as.numeric(', rules, '), "sparseVector")', collapse = ", "), ")")
    rulevars <- eval(parse(text = expr), data)
    
  } else {
    expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
    rulevars <- eval(expr, data)
    
  }
  
  colnames(rulevars) <- names(rules) <- paste0("rule", 1:length(rules))
  
  ## Remove duplicate rules
  if (removeduplicates) {
    # Remove rules with identical support
    duplicates <- duplicated(rulevars, MARGIN = 2)
    duplicates.removed <- rules[duplicates]
    rulevars <- rulevars[, !duplicates, drop = FALSE]
    rules <- rules[!duplicates]
    
  } else 
    duplicates.removed <- NULL
  
  ## Remove complement rules
  if (removecomplements) {
    if (sparse) {
      is_all_binary <- all(rulevars@x %in% 0:1)
      if (!is_all_binary)
        stop("method not implemented for non-binary rules")
      
      # find columns which length of non-zero entries are equal to the number 
      # of rows. The latter are complements if the union of the row indices are
      # equal to the number of rows
      js <- rep(1:ncol(rulevars), diff(rulevars@p))
      is <- rulevars@i + 1
      n <- rulevars@Dim[1]
      p <- rulevars@Dim[2]
      
      Jfac <- factor(js, levels = 1:p)
      row_indices <- split(is, Jfac)
      lengths <- table(Jfac)
      complements <- logical(p)
      for (i in 1:(p - 1)) {
        if (complements[i])
          next
        
        is_potential <- which(lengths[i] + lengths[(i + 1):p] == n) + i
        is_potential <- is_potential[!complements[is_potential]]
        
        if(length(is_potential) == 0)
          next
        
        union_len <- sapply(
          mapply(union, x = row_indices[i], 
                 y = row_indices[is_potential], SIMPLIFY = FALSE), 
          length)
        is_compl <- which(union_len == n) 
        if(length(is_compl) > 0)
          complements[is_potential[is_compl]] <- TRUE
      }
      
    } else {
      if (is.logical(rulevars)){
        tmp1 <- rulevars
        tmp2 <- !rulevars
        
      } else {
        is_binary <- which(apply(rulevars, 2, function(x) all(x %in% 0:1)))
        tmp1 <- rulevars
        
        # turn into -1 and 1 so tmp1 == -tmp2 for the binary columns
        tmp1[, is_binary] <- 
          ifelse(tmp1[, is_binary] == 0, -1, tmp1[, is_binary])
        
        tmp2 <- -tmp1
        
      }
      
      # want to remove columns which have flipped sign or have FALSE instead of 
      # TRUE. These columns will have equal variance
      vars <- apply(tmp1, 2, var)
      vars_distinct <- lapply(
        unique(vars), function(x){ 
          idx <- which(is_almost_eq(x, vars))
          list(var = x, n = length(idx), idx = idx)
        })
      
      complements <- logical(ncol(tmp1))
      for (va in vars_distinct) {
        if(va$n < 2L)
          next
        
        idx <- va$idx
        idx <- setdiff(idx, which(complements))
        
        if(length(idx) < 2)
          next
        
        n_idx <- length(idx)
        for(j in 1:(n_idx - 1)){
          if(complements[idx[j]])
            next
          
          this_val <- tmp1[, idx[j]]
          is_compl <- 
            which(apply(
              tmp2[, idx[(j + 1):n_idx], drop = FALSE], 2, 
              function(x) isTRUE(all.equal.numeric(this_val, x)))) + j
          
          if(length(is_compl) > 0)
            complements[idx[is_compl]] <- TRUE
        }
      }
      
    }
    
    complements <- which(complements)
    complements.removed <- rules[complements]
    if (length(complements) > 0)
      rules <- rules[-complements]
    
  } else 
    complements.removed <- NULL
  
  ## Return results:
  if (return.dupl.compl) {
    return(list(rules = rules, 
                duplicates.removed = duplicates.removed,
                complements.removed = complements.removed))
  } else {
    return(rules)
  }
  
}

# see https://stackoverflow.com/a/51457395/5861244
duplicated.dgCMatrix <- function (dgCMat, MARGIN) {
  MARGIN <- as.integer(MARGIN)
  n <- nrow(dgCMat)
  p <- ncol(dgCMat)
  J <- rep(1:p, diff(dgCMat@p))
  I <- dgCMat@i + 1
  x <- dgCMat@x
  if (MARGIN == 1L) {
    ## check duplicated rows
    names(x) <- J
    RowLst <- split(x, I)
    is_empty <- setdiff(1:n, I)
    result <- duplicated.default(RowLst)
  } else if (MARGIN == 2L) {
    ## check duplicated columns
    names(x) <- I
    ColLst <- split(x, J)
    is_empty <- setdiff(1:p, J)
    result <- duplicated.default(ColLst)
  } else {
    warning("invalid MARGIN; return NULL")
    result <- NULL
  }
  
  if (any(is_empty)){
    out <- logical(if (MARGIN == 1L) n else p)
    out[-is_empty] <- result
    if (length(is_empty) > 1)
      out[is_empty[-1]] <- TRUE
    result <- out
  }
  
  result
}

# see https://stackoverflow.com/a/30089750/5861244
cbind_sparse_vec <- function (...) {
  args <- list(...)
  stopifnot(all(sapply(args, inherits, what = "dsparseVector")))
  
  un_lengths <- unique(sapply(args,length))
  stopifnot(length(un_lengths) == 1)
  
  return(sparseMatrix( 
    x = unlist(lapply(args, slot, "x")), 
    i = unlist(lapply(args, slot ,"i")), 
    p = c(0, cumsum(sapply(args, function(x) { length(x@x) } ))),
    dims=c(un_lengths, length(args))))
}


##########
##
## Internal function to remove duplicate and complement rules:
##
delete_duplicates_complements_old <- function(rules, data, removecomplements = TRUE, 
                                          removeduplicates = TRUE, 
                                          return.dupl.compl = FALSE) {
  
  
  ## Generate rule variables:
  rulevars <- sapply(rules, function(x) {with(data, eval(parse(text = x)))})
  colnames(rulevars) <- names(rules) <- paste0("rule", 1:length(rules))
  
  ## Remove duplicate rules:
  if (removeduplicates) {
    # Remove rules with identical support:
    duplicates <- duplicated(rulevars, MARGIN = 2)
    duplicates.removed <- rules[duplicates]
    rulevars <- rulevars[, !duplicates, drop = FALSE]
    rules <- rules[!duplicates]
  } else {
    duplicates.removed <- NULL
  }
  
  
  ## Remove complement rules:
  if (removecomplements) {
    # Get variance of each rule:
    vars <- apply(rulevars, 2, var_bin)
    vars_distinct <- sapply(
      unique(vars), function(x) c(x, sum(is_almost_eq(x, vars))))
    complements <- vector(mode = "logical", length(vars))
    for (i in 1:ncol(vars_distinct)) {
      if (vars_distinct[2, i] < 2)
        next
      
      indices <- which(is_almost_eq(vars_distinct[1, i], vars))
      
      for (j in 2:length(indices)) {
        indices_prev <- indices[1:(j - 1)] 
        complements[indices_prev] <- 
          complements[indices_prev] | apply(
            rulevars[, indices_prev, drop = FALSE] != rulevars[, indices[j]], 2, all)
      }
    }  
    complements <- which(complements)
    complements.removed <- rules[complements]
    if (length(complements) > 0) {
      rulevars <- rulevars[, -complements, drop = FALSE]
      rules <- rules[-complements]
    }
  } else {
    complements.removed <- NULL
  }
  
  ## Return results:
  if (return.dupl.compl) {
    return(list(rules = rules, 
                duplicates.removed = duplicates.removed,
                complements.removed = complements.removed))
  } else {
    return(rules)
  }
  
}


##########
##
## Internal function to get variance of binary variable / rule (faster than function sd):
##
var_bin <- function(x) {
  p <- mean(x)
  p*(1L-p)
}


##########
##
## Internal function to check for near equality:
##
is_almost_eq <- function(x, y, tolerance = sqrt(.Machine$double.eps)) {
  stopifnot(is.numeric(x), length(x) == 1L)
  x_abs <- abs(x)
  xy <- if (x_abs > tolerance) {abs(x - y) / x_abs} else {abs(x - y)}
  xy <= tolerance
}

