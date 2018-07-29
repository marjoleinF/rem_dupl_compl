source("delete_duplicates_complements.R")
library("Matrix")

## Set up large matrix with random 0-1 coded variables:
n <- 10000 # nobs
p <- 10000 # nvars
set.seed(1)
data <- data.frame(replicate(p/2, sample(c(TRUE, FALSE), size = n, 
                                         replace = TRUE)))
data <- cbind(data, 
              !data[ , sample(1:(p/2), size = p/4)], # adds complements
              data[ , sample(1:(p/2), size = p/4)]) # adds duplicates
names(data) <- paste0("x", 1:p)
rules <- paste0(paste0("x", 1:p, " == TRUE"))


## Count time for evaluating rule as in old delete_duplicates_complements():
system.time(
  rulevars <- sapply(rules, function(x) {with(data, eval(parse(text = x)))}))
#user  system elapsed 
#41.12    1.01   60.97
#40.06    0.86   50.89 
rm(rulevars)


## Time for evaluating rules as in non-sparse new delete_duplicates_complements():
system.time({expr <- parse(text = paste0("cbind(", 
                                         paste0(rules, collapse = ", "), ")"))
  rulevars <- eval(expr, data)})
#user  system elapsed
#1.78    0.00    1.86 
#1.69    0.00    1.86 
rm(rulevars)


## Time for evaluating rules as in sparse new delete_duplicates_complements():
system.time({expr <- paste0("cbind_sparse_vec(", paste0(
  'as(as.numeric(', rules, '), "sparseVector")', collapse = ", "), ")")
rulevars <- eval(parse(text = expr), data)})
#user  system elapsed
#23.76    2.24   69.83 
#23.65    1.67   51.35 
rm(rulevars)

## Conclusion: old approach very slow, new non-sparse approach very fast, 
##  new sparse approach not so fast,



## compare speed and output of old and new delete_duplicates_complements():

## Old function:
system.time({
tmp_old1 <- delete_duplicates_complements_old(rules = rules, data = data, 
                                  return.dupl.compl = FALSE)})
#user  system elapsed
#203.88    2.74  247.02 
length(tmp_old1)
system.time({
tmp_old2 <- delete_duplicates_complements_old(rules = rules, data = data, 
                                  return.dupl.compl = TRUE)})
#user  system elapsed
#204.88    2.29  244.97
length(tmp_old2$rules)
table(tmp_old1 == tmp_old2$rules)



## New function non-sparse:
system.time({
tmp_new_std1 <- delete_duplicates_complements(rules = rules, data = data, 
                              return.dupl.compl = FALSE, 
                              sparse = FALSE)})
#user  system elapsed
#282.52    1.37  315.73 
system.time({
tmp_new_std2 <- delete_duplicates_complements(rules = rules, data = data, 
                              return.dupl.compl = TRUE, 
                              sparse = FALSE)})
#user  system elapsed
# 279.44    3.52  339.52 
## Check whether results match up:
table(tmp_new_std1 == tmp_new_std2$rules)
## Compare with old function:
table(tmp_new_std1 == tmp_old1) # Only one match out of 5000
table(tmp_new_std1 %in% tmp_old1) # 2500 out of 5000 match
table(tmp_new_std2$rules == tmp_old2$rules) # Only one match out of 5000
table(tmp_new_std2$rules %in% tmp_old2$rules) # 2500 out of 500 match
table(tmp_new_std2$duplicates.removed == tmp_old2$duplicates.removed) # all match
table(tmp_new_std2$complements.removed %in% tmp_old2$complements.removed) # none match



## New function sparse:
system.time({
tmp_new_sparse1 <- delete_duplicates_complements(rules = rules, data = data, 
                              return.dupl.compl = FALSE, 
                              sparse = TRUE)})
#user  system elapsed
#159.94    8.69  343.30
system.time({
tmp_new_sparse2 <- delete_duplicates_complements(rules = rules, data = data, 
                              return.dupl.compl = TRUE, 
                              sparse = TRUE)})
#user  system elapsed
#160.50   10.09  341.38 
table(tmp_new_sparse1 == tmp_new_sparse2$rules)  # all match
table(tmp_new_std1 == tmp_new_sparse1) # all match
table(tmp_new_std2$rules == tmp_new_sparse2$rules) # all match
table(tmp_new_sparse2$rules %in% tmp_old2$rules) # 2500 out of 500 match
table(tmp_new_sparse2$duplicates.removed == tmp_old2$duplicates.removed) # all match
table(tmp_new_sparse2$complements.removed %in% tmp_old2$complements.removed) # none match



## Check differences in rules saved:
length(rules)
hist(as.numeric(gsub("rule", "", names(tmp_old1), fixed = TRUE)), 
     main = "old")
hist(as.numeric(gsub("rule", "", names(tmp_new_std1), fixed = TRUE)), 
     main = "new non-sparse")
hist(as.numeric(gsub("rule", "", names(tmp_new_sparse1), fixed = TRUE)), 
     main = "new sparse")


## Conclusion: New delete_duplicates_complements() is slower than old one. 
##  Performance of new non-sparse and sparse delete_duplicates_complements()
##  is similar. 
##
## Old delete_duplicates_complements() keeps 'later' rules, 
##  new delete_duplicates_complements() keeps 'earlier' rules.
## Fastest solution may be to take old delete_duplicates_complements() and use 
##  new approach for evaluating rules.