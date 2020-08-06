



# Created 22.07.2020


# Main analysis, sensitivity 3 and 4 and post-hoc



##########################################################################################


# libraries
library(INLA)
library(dplyr)
library(sf)
library(maptools)
library(parallel)


# Read main covariates file
findata <- readRDS("data/200722_FindataGit")
colnames(findata)


#------------------------------------------------------------------------------------------------
# Explanation of dataframe

# lsoa11cd:               The LSOA code as of 2011
# Temperature:            Mean temperature during 2014-2018 [oC]
# RelativeHumidity:       Mean relative during 2014-2018 [%]
# diabetes:               Diabetes mellitus prevalence as of 2018-2019
# hypertension:           Hypertension prevalence as of 2018-2019
# obesity:                Obesity prevalence as of 2018-2019
# smoking:                Smoking prevalence as of 2018-2019
# copd:                   Chronic obstructive pulmonary disease prevalence as of 2018-2019
# high_risk_occ:          Proportion of worker in high risk COVID19 exposure occupations as of 2011
# log.pop:                Logged population per LSOA as of 2018
# id:                     An id for each LSOA
# IMD:                    Index of multiple deprivation in quintiles as of 2011
# TotalICUBeds            Number of intensive care unit beds per population as of February 2020
# days.diff               Number of days since first reported COVID19 case per LSOA (up to June 30, 2020)
# NumberCases             Number of tested positive COVID19 cases per LTLA (up to June 30, 2020)
# no2.weighted            Population weighted mean during 2014-2018 NO2 exposure [??g/m^3]
# pm25.weighted           Population weighted mean during 2014-2018 PM2.5 exposure [??g/m^3]
#------------------------------------------------------------------------------------------------


# select pollutant

pol <- "no2"
pol <- "pm25"

if(pol == "no2"){
  colnames(findata)[colnames(findata) %in% "no2.weighted"]  <- "pol"
}

if(pol == "pm25"){
  colnames(findata)[colnames(findata) %in% "pm25.weighted"]  <- "pol"
}

# Compute quintiles
findata$pol.quintiles <- cut(findata$pol, breaks = quantile(findata$pol, probs = 
                                                              seq(from = 0, to = 1, by = 0.2), na.rm = TRUE),
                             labels = paste0("Q", 1:5), include.lowest = TRUE)


# Scale the covariates
cov2scale <- c("days.diff", "NumberCases", "TotalICUBeds", "Temperature", "RelativeHumidity",
               "log.pop", "high_risk_occ", "smoking", "obesity", "copd", "diabetes", "hypertension")

findata <- as(findata, "Spatial")
findata@data[,cov2scale] <- apply(findata@data[,cov2scale], 2, scale)
findata <- st_as_sf(findata)

# load the downscaling samples
dat.sample <- readRDS("data/SamplesDownscaleCOVID")


# Define formulas for the analysis

form.anal <- list(
  
  #------- Main Analysis
  # Model 1
  cmb1 = deaths ~ 1 + pol, 
  
  # Model 2
  cmb2 = deaths ~ 1 + pol +  
    f(id, model='bym2', graph='W.adj', scale.model = TRUE, 
      constr = TRUE, hyper = hyper.bym),
  
  # Model 3
  cmb3 = deaths ~ 1 + pol +  
    factor(IMD) + days.diff + NumberCases + 
    TotalICUBeds + Temperature + RelativeHumidity + 
    log.pop + factor(urbanicity) + high_risk_occ + smoking + 
    obesity,
  
  # Model 4
  cmb4 = deaths ~ 1 + pol +  
    # confounders
    factor(IMD) + days.diff + NumberCases + 
    TotalICUBeds + Temperature + RelativeHumidity + 
    log.pop + factor(urbanicity) + high_risk_occ + smoking + 
    obesity + 
    # LF
    f(id, model='bym2', graph='W.adj', scale.model = TRUE, 
      constr = TRUE, hyper = hyper.bym), 
  
  #-------Sensitivity 3: 
  # No spread
  cmb5 = deaths ~ 1 + pol +  
    # confounders
    factor(IMD) + 
    TotalICUBeds + Temperature + RelativeHumidity + 
    log.pop + factor(urbanicity) + high_risk_occ + smoking + 
    obesity + 
    # LF
    f(id, model='bym2', graph='W.adj', scale.model = TRUE, 
      constr = TRUE, hyper = hyper.bym), 
  
  # Only spread
  cmb6 = deaths ~ 1 + pol +  
    # confounders
    days.diff + NumberCases +
    # LF
    f(id, model='bym2', graph='W.adj', scale.model = TRUE, 
      constr = TRUE, hyper = hyper.bym), 
    
  #-------Sensitivity 4, quintiles: 
  # Model 1
  cmb7 = deaths ~ 1 + factor(pol.quintiles), 
    
  # Model 2
  cmb8 = deaths ~ 1 + factor(pol.quintiles) +  
    # LF
    f(id, model='bym2', graph='W.adj', scale.model = TRUE, 
      constr = TRUE, hyper = hyper.bym), 
    
  # Model 3
  cmb9 = deaths ~ 1 + factor(pol.quintiles) +  
      # confounders
      factor(IMD) + days.diff + NumberCases + 
      TotalICUBeds + Temperature + RelativeHumidity + 
      log.pop + factor(urbanicity) + high_risk_occ + smoking + 
      obesity, 
    
  # Model 4
  cmb10 = deaths ~ 1 + factor(pol.quintiles) +  
      # confounders
      factor(IMD) + days.diff + NumberCases + 
      TotalICUBeds + Temperature + RelativeHumidity + 
      log.pop + factor(urbanicity) + high_risk_occ + smoking + 
      obesity +
    # LF
    f(id, model='bym2', graph='W.adj', scale.model = TRUE, 
      constr = TRUE, hyper = hyper.bym), 
  
  #-------Post-hoc comorbidities 
  # Model 4
  cmb11 = deaths ~ 1 + 
    pol +  
    # confounders
    factor(IMD) + days.diff + NumberCases + 
    TotalICUBeds + Temperature + RelativeHumidity +
    log.pop + factor(urbanicity) + high_risk_occ + smoking + 
    obesity + copd_nopol + hyper_nopol + diab_nopol +
    # LF
    f(id, model='bym2', graph='W.adj', scale.model = TRUE, 
      constr = TRUE, hyper = hyper.bym)
   #-------
)



# define the function to run INLA

in.mod.Poisson <- function(X, dat.inla){
  
  Y <- inla(formula = as.formula(X),
            data=dat.inla,
            family="poisson",
            E = expected,
            verbose = FALSE, 
            
            # priors fixed effects
            control.fixed = list(
              
              mean = list(pol = 0, 
                          `factor(pol.quintiles)2` = 0, `factor(pol.quintiles)3` = 0, 
                          `factor(pol.quintiles)4` = 0, `factor(pol.quintiles)5` = 0,
                          days.diff = 0, NumberCases = 0, TotalICUBeds = 0, 
                          Temperature = 0, RelativeHumidity = 0, log.pop = 0, 
                          diabetes = 0, hypertension =0, copd = 0, 
                          obesity = 0, smoking = 0, 
                          `factor(IMD)2` = 0, `factor(IMD)3` = 0, `factor(IMD)4` = 0, 
                          `factor(IMD)5` = 0, `factor(urbanicity)urban` = 0
              ), 
              
              prec = list(pol = 0.01, 
                          `factor(pol.quintiles)2` = 0.01, `factor(pol.quintiles)3` = 0.01, 
                          `factor(pol.quintiles)4` = 0.01, `factor(pol.quintiles)5` = 0.01,
                          days.diff = 0.01, NumberCases = 0.01, TotalICUBeds = 0.01, 
                          Temperature = 0.01, RelativeHumidity = 0.01, log.pop = 0.01,
                          diabetes = 0.01, hypertension =0.01, copd = 0.01,
                          obesity = 0.01, smoking = 0.01, 
                          `factor(IMD)2` = 0.01, `factor(IMD)3` = 0.01, `factor(IMD)4` = 0.01, 
                          `factor(IMD)5` = 0.01, `factor(urbanicity)urban` = 0.01
              ),
              
              prec.intercept = 1, mean.intercept = 0),
            
            control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE), 
            control.predictor=list(link = 1),
            control.inla=list(strategy="simplified.laplace", int.strategy="eb"),
            control.mode=list(restart = T, theta = thet), 
            num.threads = 1
  )
  
  Y <- inla.rerun(Y)
  return(Y)
}



# neighbor matrix
W.nb <- poly2nb(findata)
nb2INLA("W.adj", W.nb) 

# bym prior
hyper.bym <- list(theta1 = list('PCprior', c(1,0.01)), theta2 = list('PCprior', c(0.5, 0.5)))
hyper.rw<- list(theta = list(prior="pc.prec", param=c(1,0.01), prec=list(initial=log(10^-2), fixed=TRUE)))

# k is the an indicator to run on the different samples
k <- as.list(1:5)


t_0 <- Sys.time()

# The following loop is computationally expansive. If RAM explodes, you can decrease the number of cores for parallelisation.
# In addition, the parallel procedure here works for Windows machines, but can easily extended to linux or macOS using mclapply.

for(i in 1:8){
  
  print(paste("Fitting model ", i, " of ", 8))
  
  if(i %in% c(2,4:8)){
    thet  <- c(2.136760, 4.782639) # here I give good starting values for the hyperparameters of the latent field to decrease computation time.
                                   # These values can be obtained by the inla.object once you run it: inla.object$mode$theta
                                          
  }else{
    thet <- NULL
  }
  
  form <- form.anal[[i]]
  
  inla.parallel <- function(dat.inla){
    return(in.mod.Poisson(X = form, dat.inla = dat.inla))
  } 
  
  # parallel function
  par.fun <- function(k){
    dat.K <- left_join(findata, dat.sample[[k]], by = c("lsoa11cd" = "LSOA"))
    return(inla.parallel(dat.inla = dat.K))
  }
  
  
  # Set up parallel environment
  ncores <- 10
  cl_inla <- makeCluster(ncores, methods=FALSE)
  
  # extract packages on parallel environment 
  clusterEvalQ(cl_inla, {
    library(INLA)
    library(dplyr)
    library(sf)
  })
  
  # extract R objects on parallel environment
  clusterExport(cl_inla, c("inla.parallel", "form", "hyper.bym", "hyper.rw","in.mod.Poisson", "par.fun",
                           "findata", "W.nb", "dat.sample", "k", "thet"))
  
  # run the the function in parallel
  outpar <- parLapply(cl = cl_inla, k, par.fun)
  
  
  # close parallel environment
  stopCluster(cl_inla)
  
  # Propagate the uncertainty of the sampling
  fin_mod <- inla.merge(outpar)
  
  # and store
  saveRDS(fin_mod, file = paste0("data/BMA_CMB_", i, "_", pol))
  
}

t_1 <- Sys.time()
t_1 - t_0




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
