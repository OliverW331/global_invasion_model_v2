rm(list=ls());gc()
graphics.off()


library(abind) # For abind
library(pROC)
library(data.table)

###Switches###
sw.processData = FALSE
sw.model = FALSE
sw.maxit = 5
sw.separateTaxa = FALSE
sw.predict = TRUE
sw.varswitch = "model_data/regionalized_model_data/Var-variable_switches.csv"

###Starting parameter###
#other starting parameters will be added through the variable switch
pars.start = c(
  alpha = 0.00001
)

###Output###
timestamp <- format(Sys.time(), "%Y-%m-%d-%H:%M:%S")
fdir <- tempfile(pattern = paste0("1-", timestamp, "-"), tmpdir = "output")
fdir <- 'result/all_sp_time_v4_0.757'

###Functions###
source("pp_funcs_v4.R")

###Data pre-processing###
if (sw.processData){
  source("../Scripts/1a_newData.R")
}else{
  ###import pre-processed data###
  #A: 22,15,4
  soc_eco = readRDS("model_data/regionalized_model_data/A-socioEcoDat-array.rds") 
  
  #C: 15,16082
  first_sight = readRDS("model_data/regionalized_model_data/C-firstSightings-matrix.rds") 
  
  #D: 15,15,1
  pairwise_data = readRDS("model_data/regionalized_model_data/D-pairwiseData-array.rds")
  
  #B (outdated): 6026,4
  # species_traits = readRDS("model_data/B-speciesSighitngsTraits-taxa-df.rds")
  
  #Tr:15,15,22
  trade = readRDS("model_data/Tr-baci1995to2018-array.rds") 
  
  ###Other Parameters###
  #number of countries
  nc = nrow(first_sight)
  
  #number of species
  ns = ncol(first_sight)
  
  #number of years
  last_yr = max(first_sight, na.rm=T)
  first_yr = min(first_sight, na.rm=T)
  tot_yrs = last_yr - first_yr + 1
  
  #Classify the species into Taxonomy
  # sp.taxa = lapply(unique(species_traits$taxa), FUN = function(x){which(species_traits$taxa == x)})
  # names(sp.taxa) = unique(species_traits$taxa)
  
  #Construct a Time frame: 25,6202,24
  time_array = array(NA, dim=c(nc,ns,tot_yrs))
  for(i in 1:nc){
    for(j in 1:ns){
      # The last column is year, so make it 0
      time_array[i,j,] = 1:tot_yrs
    }
  }
  
  #Generate variable needs to be included through varCSV (Variable switch)
  varCSV = read.csv(sw.varswitch, stringsAsFactor = FALSE)
  sw.vars = lapply(unique(varCSV$type), FUN = function(x){
    temp.out = varCSV$variable[varCSV$type==x & varCSV$include == TRUE]
    if(length(temp.out) > 0) return(temp.out)
    return(NA)
  })
  names(sw.vars) = unique(varCSV$type)
  
  #Adding starting parameters from the varCSV(Variable Switch)
  if(length(which(varCSV$include))){
    pars.start = c(pars.start,
                   varCSV$starting.value[varCSV$include]
    )
    names(pars.start)[2:length(pars.start)] = varCSV$variable[varCSV$include]
  }
  
  # #if we have species variables, then exclude alpha
  # if(!all(is.na(sw.vars$sp.vars))){
  #   pars.start = pars.start[names(pars.start) %in% c(unlist(sw.vars), "pp.y0", "B")]
  # }
  # #Exclude B and pp.y0 since we don't need them if there are no PP variables
  # if(all(is.na(sw.vars$pp.vars))){
  #   pars.start = pars.start[-which(names(pars.start) %in% c("pp.y0", "B"))]
  # }
  
}

### Calculate the propagule pressure###
t1 = Sys.time()
pps = mk_ppstruct(sw.vars, first_sight, soc_eco, pairwise_data, trade)
t2 = Sys.time()
cat("Generating mk_ppstruct\n")
print(t2-t1)


### Model Fitting###
if (sw.model){
  # Create the directory if it does not exist
  if (!dir.exists(fdir))dir.create(fdir)
  
  # Save metadata
  saveRDS(sw.vars, paste0(fdir, "/variable_included.rds"))
  saveRDS(pars.start, paste0(fdir, "/starting_params.rds"))
  
  #
  fname = paste0(fdir, "/opt_fit-", ifelse(sw.separateTaxa, "byTaxa-", "All-"), Sys.time())
  fname = gsub("\\s+.*", ".rds", fname)
  cat("Output model:\n")
  print(fname)
  
  pps = pps[names(pps) %in% c(na.omit(unlist(sw.vars)), "remove")]
  gc()
  
  # Default method is Nelder-Mead
  if(sw.separateTaxa){
    pps.tot = pps
    first_sight.tot = first_sight
    time_array.tot = time_array
    
    for(tt in names(sp.taxa)){
      print(tt)
      first_sight = first_sight.tot[,sp.taxa[[tt]]]
      time_array = time_array.tot[,sp.taxa[[tt]],]
      pps = lapply(pps.tot, FUN = function(x){
        # It will always be the last one!
        if(length(dim(x)) == 4){
          return(x[,,,sp.taxa[[tt]]])
        }else if(length(dim(x)) == 3){
          return(x[,,sp.taxa[[tt]]])
          stop("Incorrect dimensions")
        }
      })
      ns = length(sp.taxa[[tt]])
      
      opt.out = optim(par = pars.start, le, echo = TRUE, control = list(maxit = sw.maxit))
      # Save as we go, in case it crashes
      saveRDS(opt.out, gsub("\\.rds",paste0("-",tt,".rds"),fname))
    }
  }else{
    opt.out = optim(par = pars.start, le, echo = TRUE, control = list(maxit = sw.maxit))
    saveRDS(opt.out, fname)
  }
}





#### MODEL PREDICTION ####
if(sw.predict){
  first_sight.tot = first_sight
  pps.tot = pps
  time_array.tot = time_array
  sw.vars.tot = sw.vars
  
  lf = list.files(fdir, pattern = "opt_fit.*rds", full.names = TRUE)
  
  opt = lapply(lf, readRDS)
  names(opt) = lf
  opt.results = list()
  
  # Subset sw.vars to the 1st model
  sw.vars = lapply(sw.vars.tot, FUN = function(x){
    if(all(is.na(x))) return(NA)
    vv = x[x %in% names(opt[[1]]$par)]
    if(length(vv)==0){
      return(NA)
    }else{
      return(vv)
    }
  })
  
  pps = pps[names(pps) %in% c(na.omit(unlist(sw.vars)), "remove")]
  gc() # Otherwise we can have some issues
  
  for(i in 1:length(opt)){
    print(i)
    print(lf[i])
    temp.opt = opt[[i]]
    
    if(grepl("All", lf[i])){
      first_sight = first_sight.tot
      pps = pps.tot
      time_array = time_array.tot
    }else{
      tt = gsub("(.*[0-9]-)|(\\.rds)","",lf[i])
      first_sight = first_sight.tot[,sp.taxa[[tt]]]
      time_array = time_array.tot[,sp.taxa[[tt]],]
      pps = lapply(pps.tot, FUN = function(x){
        # It will always be the last one!
        if(length(dim(x)) == 4){
          return(x[,,,sp.taxa[[tt]]])
        }else if(length(dim(x)) == 3){
          return(x[,,sp.taxa[[tt]]])
          stop("Incorrect dimensions")
        }
      })
    }
    ns = dim(first_sight)[2]
    
    # Default is species-country (time = F), not species-country-time (time = T)
    pred = exp(predict.le(temp.opt$par, time = TRUE))
    
    ## Generate observed values using C
    obs = array(data = NA, dim(pred))
    # country j, species k
    for(j in 1:dim(first_sight)[1]){
      for(k in 1:dim(first_sight)[2]){
        if(!is.na(first_sight[j,k])){ #this can't be is.na - if uninvaded, still important to know it remains a zero.
          #### FILL IN ####
          obs[1:first_sight[j,k],j,k] = 0 # If C[j,k] is 1, then it will be replaced by 1 below
          obs[first_sight[j,k],j,k] = 1
        }
      }
    }
    
    ## Melt pred values onto observed values
    obs2 = reshape2::melt(obs)
    pred2 = reshape2::melt(pred)
    
    ## Note: pred accounts for first invasions so we don't have to for obs
    pvo = data.frame(year = pred2[,1], country = pred2[,2], species = pred2[,3], pred = pred2[,4], obs = obs2[,4])
    pvo$obs[which(is.na(pvo$obs))] = 0
    pvo$pred = ifelse(pvo$obs==0, 1-pvo$pred, pvo$pred)
    pvo = na.omit(pvo)
    
    # Print alpha generated from mean(obs)
    print(-log(1-mean(pvo$obs)))
    
    opt.results[[i]] = list(converged = ifelse(temp.opt$convergence==0, TRUE, FALSE),
                            ll = temp.opt$value, 
                            dev = dev.expl(pvo$pred, pvo$obs), 
                            auc = as.numeric(roc(pvo$obs~pvo$pred)$auc),
                            par = temp.opt$par)
  }
  
  names(opt.results) = 'result'
  saveRDS(opt.results, paste0(fdir, "/results.rds"))
  saveRDS(pvo, "result/all_sp_time_v3/pvo_df.rds")
  print(opt.results)
}

###Dataframe of the result###

if(1){
  test = glm(obs ~ year, pvo, family = "binomial")
  pvo$logit.pred = test$fitted.values
  
  (test$null.deviance-test$deviance)/test$null.deviance
  roc(pvo$obs ~ pvo$logit.pred)$auc
  
  # Plot curves
  pvo2 = pvo[!duplicated(pvo[,c("year", "obs")]),]
  pvo2 = pvo2[order(pvo2$year),]
  plot(obs~year, pvo2, pch = 20)
  lines(logit.pred~year,pvo2, lwd = 2, col = "blue")
  lines(pred~year,pvo2, lwd = 2, col = "red")
  
  # pvo[, 2] <- rownames(first_sight)[pvo[, 2]]
  saveRDS(list(
    plot_data = pvo2,
    dev = (test$null.deviance-test$deviance)/test$null.deviance,
    auc = roc(pvo$obs ~ pvo$logit.pred)$auc
  ),"result/all_sp_time_v3/logit_results.rds")
  
}



