


# Last updated 23.07.2020


# Weights to downscale the COVID19 to LSOAs



############################################################################################################


library(parallel)
library(dplyr)

#library(data.table)
library(pbapply)

# import population data 
dat <- readRDS("data/pop_df_GIT")
dat = as_tibble(dat)

dat = dat %>% group_by(LTLA, Age, Ethnicity, Sex) %>% mutate(PopTot2018 = sum(PopEst2018), 
                                                             pi = PopEst2018/PopTot2018  )
dat$pi[is.na(dat$pi)] = 0
dat$Deaths[is.na(dat$Deaths)] = 0


# weigth deaths by population weights
dat$deathsd <- dat$Deaths*dat$pi


# check |sex| * |ethnicity| * |ageclasses| * |LTLA| = dat$ID -> check that all combinations of age, 
# sex ethnicity and LTLA are considered in the weigths from the deaths file
max(dat$ID) == 2*5*20*317

IDs <- unique(dat$ID)
#dat <- as.data.table(dat)
dat$IDsample <- 1:nrow(dat)



# the following parallelization is written for Windows and can be simplified for Linux/mac machines
# use 10 cores for the selection step 
n_cores <- 10
cl<-makeCluster(n_cores)

# load objects in the newly defined cluster
clusterExport(cl, "dat")
clusterExport(cl, "IDs")
clusterEvalQ(cl, library(dplyr))


## Selection Step: 
# create a list of dataframes, one for each sex-age-ethnicity-ltla combination,
# to be indipendently sampled from

dd = parLapply(cl, IDs, function(id) dat %>% dplyr::filter(ID == id) %>% dplyr::select(deathsd,pi, IDsample))

# for linux/mac 
# dd = parallel::mclapply(IDs, function(id) dat[ID == id, c('deathsd', 'pi', 'IDsample')], mc.cores = n_cores)

stopCluster(cl)


dn = dat %>% group_by(ID) %>% summarise( nlist = n(), ndeaths = mean(Deaths) )


set.seed(11)

n_samples = 10 #number of samples   
list.loop = list()

## Sampling Step:
# for each sex-age-ethnicity-ltla combination sample ndeaths LSOAs, where 
# ndeaths is the number of deaths in the selected sex-age-ethnicity-ltla group


for(j in 1:n_samples){
  
  print(paste("Extracting sample ", j ," of ", n_samples))
  list_sample = list()
  foo.list = list()
  for(i in 1:nrow(dn)){
    dat.loop <- dd[[i]]
    ndeaths <- dn[i,]$ndeaths
    nlist   <- dn[i,]$nlist
    n <- round(sum(dat.loop$deathsd))

    
    if(n != 0){
      # sample positions of the LSOAs        
      list_sample[[i]] = sample(x = 1:nrow(dat.loop), 
                                size = n, replace = TRUE, 
                                prob = dat.loop$pi)
      

    }else{
      
    # if there are no deaths in the selected sex-age-ethnicity-ltla group return NA
      list_sample[[i]] <- NA
    }
    
    foo.list[[i]] = dd[[i]]$IDsample[list_sample[[i]]]
    
  }
  
  list.loop[[j]] = foo.list
}





# replace samples with only NAs with one NA only to reduce storing space
for(i in 1:n_samples){
  for(j in 1:length(list.loop[[i]])){
    if(mean(is.na(list.loop[[i]][[j]])==1)) list.loop[[i]][[j]] = NA
  }
}


# Now calculate the age*sex*ethnicity specific expected based on the different samples
library(data.table)


dat <- as.data.table(dat)
dat$IDsample <- 1:nrow(dat)

# clean file and remove NAs
list.int <- lapply(list.loop, function(X){
  Y <- do.call(c, X)
  Y <- Y[complete.cases(Y)]
  return(Y)
})

# aggregate per ID
dat.s <- lapply(list.int, function(X){
  
  dat.sample <- data.frame(ID = X, cases = 1)
  dat.sample <- aggregate(dat.sample$cases, by = list(ID = dat.sample$ID), sum)
  colnames(dat.sample)[2] <- "n.deaths"
  return(dat.sample)
})

# and now calculate the expected, needs ~1/2hour
dat.expected <- pblapply(dat.s, function(X){
  
  dat.tmp <- left_join(dat, X, by = c("IDsample" = "ID"))
  dat.tmp$n.deaths[is.na(dat.tmp$n.deaths)] <- 0
  
  # get the death rate by age sex and ethinicity
  dat.marg <- aggregate(cbind(dat.tmp$PopEst2018, dat.tmp$n.deaths), 
                        by = list(Ethnicity = dat.tmp$Ethnicity, Sex = dat.tmp$Sex, Age = dat.tmp$Age), 
                        sum)
  
  colnames(dat.marg)[c(4:5)] <- c("Pop", "n.deaths.sum")
  dat.marg$death.rate4exp <- dat.marg$n.deaths/dat.marg$Pop
  
  dat.tmp <- left_join(dat.tmp, dat.marg, by = c("Ethnicity", "Sex", "Age"))
  
  dat.tmp$expected <- dat.tmp$PopEst2018*dat.tmp$death.rate4exp
  
  dat.fin <- aggregate(cbind(dat.tmp$n.deaths, dat.tmp$expected), by = list(LSOA = dat.tmp$LSOA), sum)
  colnames(dat.fin)[c(2,3)] <- c("deaths", "expected") 
  
  return(dat.fin)
})


# check
sum(sapply(dat.expected, function(X) sum(X$deaths) == sum(X$expected))) == n_samples
# correct

# store it
saveRDS(dat.expected, file = "data/SamplesDownscaleCOVID")



############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################


