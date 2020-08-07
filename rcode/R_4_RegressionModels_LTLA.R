
# Created 22.07.2020


# Main analysis for data aggregated at Lower Tier Local Authority level



##########################################################################################


# libraries
library(INLA)
library(dplyr)
library(sf)
library(maptools)
library(parallel)
library(spdep)

# Read main covariates file
findata <- readRDS("data/FindataGit_LTLA")
colnames(findata)


#------------------------------------------------------------------------------------------------
# Explanation of dataframe

# LTLA:                   The LTLA code as of 2019
# Temperature:            Population weighted mean temperature during 2014-2018 [oC]
# RelativeHumidity:       Population weighted mean relative humidity during 2014-2018 [%]
# obesity:                Population weighted obesity prevalence as of 2018-2019
# smoking:                Population weighted smoking prevalence as of 2018-2019
# high_risk_occ:          Proportion of worker in high risk COVID19 exposure occupations as of 2011
# log.pop:                Logged population per LTLA as of 2018
# id:                     An id for each LTLA
# IMD:                    Index of multiple deprivation in quintiles as of 2011
# TotalICUBeds            Number of intensive care unit beds per population as of February 2020
# days.diff               Number of days since first reported COVID19 case per LTLA (up to June 30, 2020)
# NumberCases             Number of tested positive COVID19 cases per LTLA (up to June 30, 2020)
# no2.weighted            Population weighted mean during 2014-2018 NO2 exposure [??g/m^3]
# pm25.weighted           Population weighted mean during 2014-2018 PM2.5 exposure [??g/m^3]
#------------------------------------------------------------------------------------------------


# compute expected


pop_df = readRDS("data/pop_df_GIT")

pop_df$LSOA = NULL
pop_df$ID = NULL
pop_df_LTLA = pop_df %>% group_by(Age, Sex, Ethnicity, LTLA) %>% summarise(PopEst2018 = sum(PopEst2018), Deaths = mean(Deaths))

pop_df_LTLA = pop_df_LTLA %>% group_by(Age, Ethnicity, Sex) %>% 
  mutate(PopTot   = sum(PopEst2018),
         RatesTot = sum(Deaths)/PopTot,
         Eijk     = RatesTot* PopEst2018) 

pop_df_LTLA$RatesTot[is.na(pop_df_LTLA$RatesTot)] <- 0
pop_df_LTLA$Eijk[is.na(pop_df_LTLA$Eijk)] <- 0
expected_df = pop_df_LTLA %>% group_by(LTLA) %>% summarise(expected = sum(Eijk), deaths = sum(Deaths))


# select pollutant

pol <- "no2"
pol <- "pm25"

if(pol == "no2"){
  colnames(findata)[colnames(findata) %in% "no2.weighted"]  <- "pol"
}

if(pol == "pm25"){
  colnames(findata)[colnames(findata) %in% "pm25.weighted"]  <- "pol"
}



# Scale the covariates
cov2scale <- c("days.diff", "NumberCases", "TotalICUBeds", "Temperature", "RelativeHumidity",
               "log.pop", "high_risk_occ", "smoking", "obesity")

findata <- as(findata, "Spatial")
findata@data[,cov2scale] <- apply(findata@data[,cov2scale], 2, scale)
findata <- st_as_sf(findata)

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
      constr = TRUE, hyper = hyper.bym)
  
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
            control.mode=list(restart = T), 
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

dat.inla = left_join(findata, expected_df, by = "LTLA")

t_0 <- Sys.time()

# initialize return object
res.list = list()

for(i in 1:4){
 
  print(paste("Fitting model ", i, " of ", 4))
  
  form <- form.anal[[i]]
  
  res.list[[i]] = in.mod.Poisson(X = form, dat.inla = dat.inla)
 
  
}

t_1 <- Sys.time()
t_1 - t_0




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
