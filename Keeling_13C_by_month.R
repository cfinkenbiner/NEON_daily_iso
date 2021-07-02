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
flist <- list.files("H:/NEON_cal_carbon_bymonth_Bowling",pattern='.calibrated.h5',full.names=TRUE,recursive=TRUE)

# open csv of NEON site metadata.
neon_sites <- read_csv("C:/Users/libon/Box/190115-field-sites.csv")

#------------------------------------------------------------------
# extract list of stations that script will run through.
sites.tmp <- strsplit(flist,split=".",fixed=TRUE)   # split file names to extract part of file name w/ site ID
codes <- unique(sapply(sites.tmp,'[[',3))           # extract and unlist unique site IDs
names <- sapply(codes,
    function(x){neon_sites$`Site Name`[neon_sites$`Site ID` == x]}) # match site names to list of unique site IDs.
domain.number <- sapply(codes,
    function(x){neon_sites$`Domain Number`[neon_sites$`Site ID` == x]}) # get domain numbers
domain.name <- sapply(codes,
    function(x){neon_sites$`Domain Name`[neon_sites$`Site ID` == x]}) # get domain numbers
domain.points <- sapply(codes,
		function(x){as.numeric(strsplit(neon_sites$`Lat./Long.`[neon_sites$`Site ID` == x],",")[[1]])})

layer.isotope <- function(hdf5_file, sitename, which_part = "C13") {
  ### hdf5_file     : HDF5 file's absolute path
  ### sitename      : The name of the required neon sites, correspondent to "codes[i]" 
  ### which_part    :
    iso_con <- "rtioMoleDryCo2"
  
    if (which_part == "C13") {
        
         iso_data <- "dlta13CCo2"
      
        } else {
        
            stop("your are requiring water isotope but the 'which_part' is 
             not correctly specified keep the default if you are doing C13 keelings")
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
      all_layer <- h5read(hdf5_file, paste0('/', sitename,'/dp01/data/isoCo2/', layer_names[k]))
      ###The isotope dataset
      data_iso <- all_layer[[iso_data]][,c("mean_cal","timeBgn")] ##make sure to use mean_cal
     
      ###The Dry mole H2O dataset
      data_rtio <- all_layer[[iso_con]][,c("mean_cal","timeBgn")] ##make sure to use mean_cal

      names(data_iso) <- c(paste0(which_part, "_isoCo2"),"timeBgn")
      names(data_rtio) <- c(paste0(which_part, "_",iso_con),"timeBgn")
    
      data.df <- merge(data_iso,data_rtio,by=c("timeBgn"))
    
      data.df$timeBgn <- as.POSIXct(data.df$timeBgn,format="%Y-%m-%dT%H:%M:%OSZ",tz="UTC")
      data.df <-  data.df %>%
        mutate(dom=date(timeBgn))
      
      data_flux_height[[k]] <- data.df
      names(data_flux_height)[k] <- paste0(which_part,"(",layer_names[k],")")
      
    }
    H5close()
    return(data_flux_height)
    
}


# get tz by site.
site.tz <- sapply(codes,
                  function(x){tz_lookup_coords(as.vector(domain.points[1,x]),as.vector(domain.points[2,x]))})



  iso.all <- list() ##a list stores all iso datasets
  
  output.all.height <- list() ##a list stores all keeling isotopes
  
  fitted.all.height <- list() 
  

for (i in 1:length(codes)) {
  print("--------------------------------------------")
  print(paste0("Working on site: ",codes[i]))
  slist <- list.files("H:/NEON_cal_carbon_bymonth_Bowling", pattern=codes[i],full.names=TRUE,recursive=TRUE)

  # get just the h5 files.
  slist <- slist[!grepl(".gz",slist)]
  print("Let's look at some files of this site!")
  print(head(slist))
  print("--------------------------------------------")
  
  iso.all[[i]] <- list()
  
  output.all.height[[i]] <- list()
  

  fitted.all.height[[i]] <- list()

  # loop through sfiles and load the calibration datasets 
  for (j in 1:length(slist)) {        
	  ##################################This is the keeling for the flux using all heights
	  tower_data <- layer.isotope(slist[j], codes[i], which_part = "C13")	  
	  # the binded dataset with all tower heights
	  binded_data <- do.call(rbind, tower_data) 
	  iso.all[[i]][[j]] <- tower_data
	  # fit the model using all heights
	  fitted.all.height[[i]][[j]] <- binded_data %>% group_by(dom) %>% na.omit() %>% 
	    do(model = lm(I(C13_isoCo2*C13_rtioMoleDryCo2) ~ C13_rtioMoleDryCo2, data=.))
	  print("-----Keeling plot done------------------------")
	  cat("\n")
	  H5close()

  }
  # squash by site.
  fitted.all.height[[i]]  <- do.call(rbind,fitted.all.height[[i]])

  slopes.all.height <- vector()
  rsq.all.height <- vector()
  nptsReg.all.height <- vector()
  slopes.std <- vector()
  
  for (z in 1:nrow(fitted.all.height[[i]])) {
    slopes.all.height[z] <- coef(fitted.all.height[[i]]$model[[z]])[2]
    rsq.all.height[z]    <- summary(fitted.all.height[[i]]$model[[z]])$r.squared
    summaryTable <- coef(summary(fitted.all.height[[i]]$model[[z]]))
    if (nrow(summaryTable) == 1) {
      slopes.std[z] <- NA
      } else{
      slopes.std[z] <- coef(summary(fitted.all.height[[i]]$model[[z]]))[2,2] 
      }
    
    nptsReg.all.height[z]  <- length(fitted.all.height[[i]]$model[[z]]$residuals) 
    }
 
  output.all.height[[i]] <- data.frame(date =fitted.all.height[[i]]$dom,
                                      flux=slopes.all.height,
                                       rsq=rsq.all.height,
                                       nptsReg=nptsReg.all.height,
                                       slopeStd = slopes.std,
                                       site=rep(codes[i],nrow(fitted.all.height[[i]]))  
                                       )
                                  
}

# #squash all across sites.
output.all.height <- do.call(rbind, output.all.height)
write.csv(x = output.all.height, 
          file =  "C:/Users/libon/Box/neon_extrac_data/Scientific Data Code/keeling isotopes/et_C13_iso.csv")
