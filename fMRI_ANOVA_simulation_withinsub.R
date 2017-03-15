library(dplyr)
library(reshape2)
library(neuRosim)
library(oro.nifti)
library(parallel)


# Make stimulus timing -------------------------------------------------
  ## 30 events per condition
  ## 216 total dynamic scans == 432 seconds (7.2 minutes)
  ## 2 runs

  
  # Run 1
    onset_cond1_level1_run1 <- c(0,24,38,42,48,50,64,78,82,98,126,146,154,166,174,200,214,240,242,276,282,290,294,296,298,308,324,338,346,352) + 10  # add 10 seconds to the begining of the timeseries
    onset_cond2_level1_run1 <- c(4,20,22,30,34,52,68,72,74,88,90,100,108,136,142,144,160,176,178,188,194,210,220,234,266,288,310,322,328,340) + 10
    onset_cond1_level2_run1 <- c(18,28,54,60,106,112,120,130,134,168,172,180,186,206,208,224,244,248,258,280,284,314,316,334,342,344,348,358,366,376) + 10
    onset_cond2_level2_run1 <- c(12,32,36,62,76,80,86,94,104,114,124,132,158,184,192,196,222,226,228,256,264,268,304,318,332,354,360,368,372,374) + 10
  
  # Run 2
    onset_cond1_level1_run2 <- c(0,2,6,26,32,38,40,62,78,82,92,96,120,144,146,164,172,200,210,214,244,274,286,288,298,308,320,322,344,368) + 10
    onset_cond2_level1_run2 <- c(4,12,14,16,20,36,48,56,74,110,134,140,174,176,186,192,194,196,206,228,238,254,256,262,280,290,300,302,316,354) + 10
    onset_cond1_level2_run2 <- c(30,52,76,86,90,98,108,114,118,128,132,138,150,152,160,166,212,240,276,292,312,314,332,350,356,362,364,366,370,376) + 10
    onset_cond2_level2_run2 <- c(10,28,46,50,58,60,66,80,94,102,158,184,188,202,208,222,224,236,248,250,260,264,278,294,306,328,336,340,374,378) + 10
      

# Make subject names ---------------------------------------------------  
  subject_list <- list("sub001","sub002","sub003","sub004","sub005","sub006","sub007","sub008","sub009","sub010",
                       "sub011","sub012","sub013","sub014","sub015","sub016","sub017","sub018","sub019","sub020",
                       "sub021","sub022","sub023","sub024","sub025","sub026","sub027","sub028","sub029","sub030",
                       "sub031","sub032","sub033","sub034","sub035","sub036","sub037","sub038","sub039","sub040",
                       "sub041","sub042","sub043","sub044","sub045","sub046","sub047","sub048","sub049","sub050",
                       "sub051","sub052","sub053","sub054","sub055","sub056","sub057","sub058","sub059","sub060",
                       "sub061","sub062","sub063","sub064","sub065","sub066","sub067","sub068","sub069","sub070",
                       "sub071","sub072","sub073","sub074","sub075","sub076","sub077","sub078","sub079","sub080",
                       "sub081","sub082","sub083","sub084","sub085","sub086","sub087","sub088","sub089","sub090",
                       "sub091","sub092","sub093","sub094","sub095","sub096","sub097","sub098","sub099","sub100")



# Define simulation function----------------------------------------------------------
sim_func <- function(list){   
    
    print(list)
    
   # Simulated fMRI data generation Run 1 
    TR <- 2
    nscan <- 216
    total <- TR * nscan  # total time in seconds
    os1 <- onset_cond1_level1_run1
    os2 <- onset_cond2_level1_run1
    os3 <- onset_cond1_level2_run1
    os4 <- onset_cond2_level2_run1
    dur <- list(2, 2, 2, 2)
    os <- list(os1, os2,os3, os4)
    w <- c(0.3, 0.3, 0.01, 0.09, 0.1, 0.2)
    
    #ROIs
    roi1 <- c(16,30,37)
    roi2 <- c(44,30,37)
    roi3 <- c(17,49,23)
    roi4 <- c(42,49,23)
    roi5 <- c(28,52,32)

    regions <- simprepSpatial(regions = 5, 
                              coord = list(roi1,roi2,roi3,roi4,roi5), 
                              radius = c(3,3,3,3,3), 
                              form = "sphere")
    
    # Simulate response variability  
    effect  <- list(rnorm(1,1,1),rnorm(1,1,1),rnorm(1,2,1),rnorm(1,1,1))
    effect1 <- effect
    effect2 <- effect
    effect3 <- effect
    effect4 <- effect
    effect5 <- effect
    
    effects <- list(effect1,effect2,effect3,effect4,effect5)
    
    # Apply the design temporally
    design1 <- simprepTemporal(regions = 5, totaltime = total, onsets = list(os,os,os,os,os), 
                              durations = list(dur,dur,dur,dur,dur), effectsize = effects, TR = TR, hrf = "double-gamma")
    
    # Nifti dimentions: 60,70,55
    data <- simVOLfmri(dim = c(60, 70, 55), base = 0, design = design1,
                       image = regions, SNR = .50, noise = "mixture", type = "rician",
                       weights = w, verbose = FALSE)
    
    # Adds a base level to images for use in FSL
    data <- data + 200
    
    dir_path <- paste("/lab/neurodata/anova_sims/raw/snr_50_within/",list, sep="")
    dir.create(dir_path)
    
    path <- paste("/lab/neurodata/anova_sims/raw/snr_50_within/",list,"/",list,"_run1", sep="")
    writeNIfTI(data, path)
    
    
  # Simulated fMRI data generation Run 2 
    TR <- 2
    nscan <- 216
    total <- TR * nscan  # total time in seconds
    os1 <- onset_cond1_level1_run2
    os2 <- onset_cond2_level1_run2
    os3 <- onset_cond1_level2_run2
    os4 <- onset_cond2_level2_run2
    dur <- list(2, 2, 2, 2)
    os <- list(os1, os2,os3, os4)
    w <- c(0.3, 0.3, 0.01, 0.09, 0.1, 0.2)
    
    regions <- simprepSpatial(regions = 5, 
                              coord = list(roi1,roi2,roi3,roi4,roi5), 
                              radius = c(3,3,3,3,3), 
                              form = "sphere")
  
    # Assign effect with the same values as run 1
    effect1 <- effect
    effect2 <- effect
    effect3 <- effect
    effect4 <- effect
    effect5 <- effect
    
    effects <- list(effect1,effect2,effect3,effect4,effect5)
    
    # Apply the design temporally
    design2 <- simprepTemporal(regions = 5, totaltime = total, onsets = list(os,os,os,os,os), 
                               durations = list(dur,dur,dur,dur,dur), effectsize = effects, TR = TR, hrf = "double-gamma")
    
    # Nifti dimentions: 60,70,55
    data <- simVOLfmri(dim = c(60, 70, 55), base = 0, design = design2,
                       image = regions, SNR = .50, noise = "mixture", type = "rician",
                       weights = w, verbose = FALSE)
    
    # Adds a base level to images for use in FSL
    data <- data + 200
  
    path <- paste("/lab/neurodata/anova_sims/raw/snr_50_within/",list,"/", list,"_run2", sep="")
    writeNIfTI(data, path)
    
  }  

# Run simulation and generate data files ---------------------------------
  set.seed(10) # for reproducability
  mclapply(subject_list,sim_func, mc.cores = 30)

    
