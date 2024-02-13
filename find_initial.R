#this script is used to find the best initial value for optim
rm(list=ls());gc()
graphics.off()


library(abind) # For abind
library(pROC)
library(data.table)
library(dplyr)
library(combinat)
source("gim_background_pp_funcs_v2.R")
fst_rnd = FALSE
sec_rnd = TRUE
data_path = "model_data/2023/"
sw.varswitch = paste0(data_path, "Var-variable_switches.csv")
#every time we run the model, we need to change the initial values of each variables
#including, gdp, pop, dist, pp_bk, alpha

###import pre-processed data###
soc_eco = readRDS(paste0(data_path, "A-socioEcoDat-array.rds"))
 
first_sight = readRDS(paste0(data_path, "C_new-firstSightings-matrix.rds"))
#set.seed(123) # for reproducibility
#first_sight <- first_sight[, sample(1:ncol(first_sight), 3000)]
pairwise_data = readRDS(paste0(data_path, "D-pairwiseData-array.rds"))

#trade = readRDS("model_data/Tr-baci1995to2018-array.rds") 
#policy = readRDS("model_data/regionalized_model_data/policy.rds")

###Other Parameters###
#number of countries
nc = nrow(first_sight)

#number of species
ns = ncol(first_sight)

#number of years
last_yr = max(first_sight, na.rm=T)
first_yr = 1
tot_yrs = last_yr - first_yr + 1

#Classify the species into Taxonomy
# sp.taxa = lapply(unique(species_traits$taxa), FUN = function(x){which(species_traits$taxa == x)})
# names(sp.taxa) = unique(species_traits$taxa)

#Construct a Time frame: 25,6202,24
time_array = array(NA, dim=c(nc,ns,tot_yrs))
for(i in 1:nc){
  for(j in 1:ns){
    # The last column is year, so make it 0
    # if want to use policy, replace 1:tot_yrs with policy$NUM
    time_array[i,j,] = 1:tot_yrs
  }
}

###create the model as a function###
global_model = function(pps, pars.start, sw.model = TRUE, sw.predict = TRUE, sw.maxit = 2000){
  sw.model = sw.model
  sw.maxit = sw.maxit
  sw.predict = sw.predict
  
  ###Starting parameter###
  #other starting parameters will be added through the variable switch
  pars.start = pars.start
  
  ###Output###
  timestamp <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
  fdir <- tempfile(pattern = paste0("1-", timestamp, "-"), tmpdir = "output/choose_init/")
  
  if (sw.model){
    # Create the directory if it does not exist
    if (!dir.exists(fdir))dir.create(fdir)
    
    # Save metadata
    saveRDS(sw.vars, paste0(fdir, "/variable_included.rds"))
    saveRDS(pars.start, paste0(fdir, "/starting_params.rds"))
    
    #
    fname = paste0(fdir, "/opt_fit-All-", Sys.time())
    fname = gsub("\\s+.*", ".rds", fname)
    cat("Output model:\n")
    print(fname)
    
    pps = pps[names(pps) %in% c(na.omit(unlist(sw.vars)), "remove")]
    gc()
    
    # Default method is Nelder-Mead

    opt.out = optim(par = pars.start, le, echo = TRUE, control = list(maxit = sw.maxit))
    saveRDS(opt.out, fname)
  }
  
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
          if(!is.na(first_sight[j,k]) && first_sight[j,k] != 0){ #this can't be is.na - if uninvaded, still important to know it remains a zero.
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
      gim_aic = cal_aic(length(temp.opt$par),pvo$pred, pvo$obs)
      opt.results[[i]] = list(converged = ifelse(temp.opt$convergence==0, TRUE, FALSE),
                              ll = temp.opt$value, 
                              dev = dev.expl(pvo$pred, pvo$obs), 
                              auc = as.numeric(roc(pvo$obs~pvo$pred)$auc),
                              aic = gim_aic,
                              par = temp.opt$par)
    }
    
    names(opt.results) = 'result'
    saveRDS(opt.results, paste0(fdir, "/results.rds"))
    saveRDS(pvo, paste0(fdir, "/pvo_df.rds"))
    print(opt.results)
  }
}

### Find best init###
#Since we have the model as the function, we just need to change the input of pars.start when finding the initial values
#first is to consider which parameters are considered, this can be find in the variable switch file
#Generate variable needs to be included through varCSV (Variable switch)
varCSV = read.csv(paste0(data_path, "new_variable_switches.csv"),stringsAsFactor = FALSE)
a_bkpp = read.csv(paste0(data_path, "alpha_bkpp.csv"),stringsAsFactor = FALSE)
#var_cases = permn(varCSV$variable)
mk_sw.vars = function (varCSV){
  sw.vars = lapply(unique(varCSV$type), FUN = function(x){
    temp.out = varCSV$variable[varCSV$type==x & varCSV$include == TRUE]
    if(length(temp.out) > 0) return(temp.out)
    return(NA)
  })
  names(sw.vars) = unique(varCSV$type)
  return(sw.vars)
}

#how to prepare the pars.start
#needs to manipulate which varibles to include in the model
#iterate all combinations to find the lowest AIC

#create an algorithm to iterate all possible cases
#use permn to get all possible permutation of the variables

#alpha and pp_bk will always be included. It's just that the values will change
if (fst_rnd){
  alpha_init = a_bkpp$starting.value[which(a_bkpp$variable == "alpha")]
  pbk = a_bkpp$starting.value[which(a_bkpp$variable == "pp_bk")]
  pars.init = c(
    alpha = alpha_init,
    pp_bk = pbk
  )
  results <- data.frame(alpha=numeric(), auc = numeric(), dev = numeric(), ll = numeric(), convergence = logical())
  ### First step ###
  #to get it started, we test only one variable first in the model and start from here
  for ( i in c(1:nrow(varCSV))) {
    print("##########Preparing pars.start##########")
    varCSV$include[i] = TRUE
    sw.vars = mk_sw.vars(varCSV)
    pars.start = c(pars.init, varCSV$starting.value[varCSV$include])
    names(pars.start)[3:length(pars.start)] = varCSV$variable[varCSV$include]
    print(pars.start)
    print(sw.vars)
    print("##########Start Modeling##########")
    pps = mk_ppstruct(sw.vars, first_sight, soc_eco, pairwise_data, trade)
    model_res = global_model(pps, pars.start,TRUE,TRUE,2000)
    print("##########Preparing results##########")
    auc = model_res$result$auc
    aic = model_res$result$aic
    dev = model_res$result$dev
    ll = model_res$result$ll
    convergence = model_res$result$converged
    final_alpha = as.numeric(model_res$result$par["alpha"])
    final_pp_bk = as.numeric(model_res$result$par["pp_bk"])
    final_gdp_src = as.numeric(model_res$result$par["gdp.src"])
    final_gdp_dst = as.numeric(model_res$result$par["gdp.dst"])
    final_pop_src = as.numeric(model_res$result$par["population.src"])
    final_pop_dst = as.numeric(model_res$result$par["population.dst"])
    final_dis = as.numeric(model_res$result$par["phys_dist"])
    final_t = as.numeric(model_res$result$par["t"])
    # Store the results
    results <- rbind(results, data.frame(alpha=alpha_init, pp_bk = pbk, auc = auc, aic = aic, dev = dev, ll = ll, convergence = convergence, final_alpha = final_alpha, final_pp_bk=final_pp_bk, final_gdp_src = final_gdp_src, final_pop_src = final_pop_src, final_gdp_dst = final_gdp_dst, final_pop_dst = final_pop_dst, final_dis = final_dis, final_t = final_t))
    saveRDS(results, "result/choose_init/2023/gimbkpp20.rds")
    varCSV$include[i] = FALSE
    i = i+1
  }
}else{
  results = readRDS("result/choose_init/2023/gimbkpp20.rds")
}


#Now we have finished the first round of optim, using only one variable.
#to fit the rest of parameter, we can use var_cases to get the order of fitting
#1. determine which parameters have already been fitted.
#2. use the fitted value as the starting value
#3. determine the unfitted parameters and choose one to fit
#4. repeat the process. 

#let's do a second round
if (sec_rnd) {
  results_2nd <- data.frame(alpha=numeric(), auc = numeric(), dev = numeric(), ll = numeric(), convergence = logical())
  row_min_aic = which.min(results$aic)
  include = which(!is.na(results[row_min_aic, 10:15]))
  varCSV$include[include] = TRUE
  #now we need to decide the next parameter to include
  #gimbkpp for now, need to change to result
  varCSV$starting.value[include] = results[row_min_aic, include+9]
  alpha_init = results[row_min_aic,]$final_alpha
  pbk = results[row_min_aic,]$final_pp_bk
  pars.init = c(
    alpha = alpha_init,
    pp_bk = pbk
  )
  for (j in which(varCSV$include == FALSE)){
    varCSV$include[j] = TRUE
    print("##########Preparing pars.start##########")
    sw.vars = mk_sw.vars(varCSV)
    pars.start = c(pars.init, varCSV$starting.value[varCSV$include])
    names(pars.start)[3:length(pars.start)] = varCSV$variable[varCSV$include]
    print(pars.start)
    #print(sw.vars)
    print("##########Start Modeling##########")
    pps = mk_ppstruct(sw.vars, first_sight, soc_eco, pairwise_data, trade)
    model_res = global_model(pps, pars.start,TRUE,TRUE,2000)
    print("##########Preparing results##########")
    auc = model_res$result$auc
    aic = model_res$result$aic
    dev = model_res$result$dev
    ll = model_res$result$ll
    convergence = model_res$result$converged
    final_alpha = as.numeric(model_res$result$par["alpha"])
    final_pp_bk = as.numeric(model_res$result$par["pp_bk"])
    final_gdp_src = as.numeric(model_res$result$par["gdp.src"])
    final_gdp_dst = as.numeric(model_res$result$par["gdp.dst"])
    final_pop_src = as.numeric(model_res$result$par["population.src"])
    final_pop_dst = as.numeric(model_res$result$par["population.dst"])
    final_dis = as.numeric(model_res$result$par["phys_dist"])
    final_t = as.numeric(model_res$result$par["t"])
    # Store the results
    results_2nd <- rbind(results_2nd, data.frame(alpha=alpha_init, pp_bk = pbk, auc = auc, aic = aic, dev = dev, ll = ll, convergence = convergence, final_alpha = final_alpha, final_pp_bk=final_pp_bk, final_gdp_src = final_gdp_src, final_pop_src = final_pop_src, final_gdp_dst = final_gdp_dst, final_pop_dst = final_pop_dst, final_dis = final_dis, final_t = final_t))
    saveRDS(results_2nd, "result/choose_init/2023/gimbkpp20_rnd2.rds")
    varCSV$include[j] = FALSE
  }
  varCSV$include[include] = FALSE
  varCSV$starting.value[include] = 0
}


