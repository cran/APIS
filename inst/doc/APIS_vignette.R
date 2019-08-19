## ---- message=F, warning=F-----------------------------------------------
library(APIS)

## ------------------------------------------------------------------------
data("APIS_offspring")
data("APIS_sire")
data("APIS_dam")

head(APIS_offspring[,1:10])
rownames(APIS_offspring[1:6,])

## ------------------------------------------------------------------------
head(APIS_offspring[,1:10])
head(APIS_sire[,1:10])
head(APIS_dam[,1:10])
error <- 0.05 #I accept 5% of errors in the results

## ---- eval = F-----------------------------------------------------------
#  result <- APIS(off.genotype = APIS_offspring,
#                 sire.genotype = APIS_sire,
#                 dam.genotype = APIS_dam,
#                 error = error)

## ---- eval = F-----------------------------------------------------------
#  new.result <- personalThreshold(APIS.result = result,
#                                  method = 'Pmendel',
#                                  threshold = 0.7)

## ---- eval = F-----------------------------------------------------------
#  new.result <- personalThreshold(APIS.result = result,
#                                  method = 'exclusion',
#                                  threshold = 1)

## ----full_data, echo = F, message = F, warning = F, fig.width=8, fig.height=4----
result <- APIS(off.genotype = APIS_offspring,
               sire.genotype = APIS_sire,
               dam.genotype = APIS_dam,
               error = 0.05,
               verbose = FALSE)

## ---- eval = T, fig.width=8, fig.height=4--------------------------------
new.result <- personalThreshold(APIS.result = result,
                                method = 'exclusion',
                                threshold = 1,
                                verbose = FALSE)

## ----degraded_data, echo = F, message = F, warning = F, results = "hide", fig.width=8, fig.height=4----
result <- APIS(off.genotype = APIS_offspring[, 1:35],
               sire.genotype = APIS_sire[sample(c(1:nrow(APIS_sire)), 10), 1:35],
               dam.genotype = APIS_dam[, 1:35],
               error = 0.05,
               verbose = FALSE)

## ----info----------------------------------------------------------------
print(sessionInfo(), locale=FALSE)

