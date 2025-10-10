pkgname <- "excode"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "excode-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('excode')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("excodeFamily")
### * excodeFamily

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: excodeFamily
### Title: Create a family of probability distributions for excess count
###   detection.
### Aliases: excodeFamily

### ** Examples


excode_family_pois <- excodeFamily("Poisson")
excode_family_pois




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("excodeFamily", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("excodeFormula")
### * excodeFormula

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: excodeFormula
### Title: Create a formula of the model for excess count detection
### Aliases: excodeFormula

### ** Examples


excode_formula_har <- excodeFormula("Harmonic")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("excodeFormula", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("excodeModel")
### * excodeModel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: excodeModel
### Title: Create a model for excess count detection
### Aliases: excodeModel

### ** Examples

# Initialisation of a mean model without timetrend with Poisson emission

excode_formula_mean <- excodeFormula("Mean", timeTrend = FALSE)
excode_family_pois <- excodeFamily("Poisson")
excodeModel(excode_family_pois, excode_formula_mean)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("excodeModel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("init_excode")
### * init_excode

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: init_excode
### Title: Initialize a multi-state EXCODE model from a surveillance time
###   series
### Aliases: init_excode

### ** Examples

## Not run: 
##D # surv_ts must provide the columns referenced by the EXCODE formula.
##D mod <- init_excode(
##D   surv_ts   = my_surveillance_df,
##D   timepoint = nrow(my_surveillance_df),
##D   time_units_back = 156,           # ~3 years of weekly data
##D   distribution = "nbinom",
##D   states = 3,
##D   time_trend = "Linear",
##D   periodic_model = "Harmonic",
##D   period_length = 52,
##D   intercept = TRUE
##D )
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("init_excode", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot-excodeModel-ANY-method")
### * plot-excodeModel-ANY-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot,excodeModel,ANY-method
### Title: Summary of an excodeModel.
### Aliases: plot,excodeModel,ANY-method

### ** Examples


# TODO




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot-excodeModel-ANY-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("run_excode")
### * run_excode

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: run_excode
### Title: Detect Excess Counts in Epidemiological Time Series
### Aliases: run_excode

### ** Examples

# Create a Poisson harmonic model
excode_family_pois <- excodeFamily("Poisson")
excode_formula_har <- excodeFormula("Harmonic")
excode_har_pois <- excodeModel(excode_family_pois, excode_formula_har)

# Example: data.frame as input
data(shadar_df)
result_shadar_har <- run_excode(excode_har_pois, shadar_df, 209:295)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("run_excode", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary-excodeModel-method")
### * summary-excodeModel-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary,excodeModel-method
### Title: Summary of an excodeModel.
### Aliases: summary,excodeModel-method

### ** Examples


# Looking at summary of the results using a harmonic Poisson model on the shadar_df
## Not run: 
##D #' excode_family_pois <- excodeFamily("Poisson")
##D excode_formula_har <- excodeFormula("Harmonic")
##D excode_har_pois <- excodeModel(excode_family_pois, excode_formula_har)
##D # perform excess count detection for time points 209:295
##D result_shadar_har <- run_excode(shadar_df, excode_har_pois, 209:295)
##D # obtain the summary of the results for the time points 209:295
##D summary(result_shadar_har)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary-excodeModel-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
