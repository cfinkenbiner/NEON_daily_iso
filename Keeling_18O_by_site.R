# Keeling_CO2_batch.R
# rich fiorella 191008
# updated for batch processing on 191014.

rm(list=ls())
# make a quick keeling plot based on a data file.

library(tidyverse)
library(rhdf5)
library(lubridate)
library(gridExtra)
library(lutz)
library(sp)
library(ggrepel)

# get a list of available files.
flist <- list.files("C:/Users/libon/Box/neon_extrac_data/CalibratedWaterFiles/210405_Wiso",pattern='.h5',full.names=TRUE,recursive=TRUE)


sites.tmp <- strsplit(flist,split=".",fixed=TRUE)   


layer.isotope <- function(hdf5_file, sitename, which_part = "H2") {
  ### hdf5_file     : HDF5 file's absolute path
  ### sitename      : The name of the required neon sites, correspondent to "codes[i]" 
  ### which_part    : Either "H2" or "O18"
    iso_con <- "rtioMoleWetH2o"
  
    if (which_part == "H2") {
        
         iso_data <- "dlta2HH2o"
         
         } else if (which_part == "O18"){
           
           iso_data <- "dlta18OH2o"
      
        } else {
        
            stop("your are requiring water isotope but the 'which_part' is 
             not correctly specified. either which_part = H2 or 
             which_part = O18, or just keep the default which is  which_part = 'H2'!")
        }
      
 
    print(paste0("We are extracting isotope: ", which_part))
    print(paste0("The isotope data field being extracted is: ",iso_data))
    print(paste0("The isotope concerntration being extracted is: ",iso_con))
    
    unique_names <- unique(h5ls(hdf5_file)$name) ### the unique names of the HDF5 groups 
    layer_names <- vector()
    
    for (m in unique_names) {
      start_str <- startsWith(m, "000_0")
      end_str <- endsWith(m, "0_09m")
      if (start_str && end_str) {
        layer_names <- c(layer_names, m)
      }
    }
    data_flux_height <- list() ###
    print(paste0("The isotope height measurements are: ",layer_names))
    
    for (k in 1:length(layer_names)){
      ###read the isotope and concentration at each layer
      all_layer <- h5read(hdf5_file, paste0('/', sitename,'/dp01/data/isoH2o/', layer_names[k]))
      ###The isotope dataset
      data_iso <- all_layer[[iso_data]][,c("mean_cal","timeBgn")]
      ###The Dry mole H2O dataset
      data_rtio <- all_layer[[iso_con]][,c("mean","timeBgn")] ## no way to generate a ¡®calibrated¡¯ version of rtioMoleWetH2o 
                                                              ## because the reference measurements don¡¯t 
                                                              ## have a known value of rtioMoleWetH2o
      
      names(data_iso) <- c(paste0(which_part, "_isoH2o"),"timeBgn")
      names(data_rtio) <- c(paste0(which_part, "_",iso_con),"timeBgn")
      
      data.df <- merge(data_iso,data_rtio, by=c("timeBgn"))
      data.df$timeBgn <- as.POSIXct(data.df$timeBgn,format="%Y-%m-%dT%H:%M:%OSZ",tz="UTC")
      data.df <-  data.df %>%
        mutate(dom=date(timeBgn))
      
      data_flux_height[[k]] <- data.df
      names(data_flux_height)[k] <- paste0(which_part,"(",layer_names[k],")")
      
    }
    H5close()
    return(data_flux_height)
    
}

site.keeling <- list()


# loop through sfiles and load the calibration datasets 
for (j in 1:length(flist)) {
  site <-  strsplit(flist[j],split=".",fixed=TRUE)[[1]][3]   
  print(paste0("*******************Star working on ", site, " *****************"))
  tower_data <- layer.isotope(flist[j], site, which_part = "O18")
  binded_data <- do.call(rbind, tower_data) 
	  # fit the model using all heights
  keelings <- binded_data %>% group_by(dom) %>% na.omit() %>% 
	    do(model = lm(I(O18_isoH2o*O18_rtioMoleWetH2o) ~ O18_rtioMoleWetH2o, data=.))
  print("*******************Finished the Keeling *****************")
  
  tempdf <- data.frame(date = rep(NA, nrow(keelings)), 
                  flux =rep(NA,nrow(keelings)), 
                  rsq =rep(NA,nrow(keelings)),
                  nptsReg = rep(NA,nrow(keelings)),
                  slopeStd= rep(NA,nrow(keelings)),
                  site =rep(site,nrow(keelings)))
  
  for (rs in 1:nrow(keelings)) {
    tempdf[rs, 1] <-  toString(keelings[rs,1]$dom)
    tempdf[rs, 2] <- coef(keelings[rs,2]$model[[1]])[2]
    tempdf[rs, 3] <- summary(keelings[rs,2]$model[[1]])$r.squared
    tempdf[rs, 4] <- length(keelings[rs,2]$model[[1]]$residuals)
    summaryTable <- coef(summary(keelings[rs,2]$model[[1]]))
    if (nrow(summaryTable) == 2) {
      tempdf[rs, 5] <- summaryTable[2,2] 
    }
        
   
    
  }
  site.keeling[[j]] <- tempdf
  print("******************done**********************************")
}


output.all.height <- do.call(rbind, site.keeling)
write.csv(x = output.all.height, file = "C:/Users/libon/Box/neon_extrac_data/Scientific Data Code/keeling isotopes/et_O18_iso.csv")










