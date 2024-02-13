rm(list=ls());gc()
graphics.off()

library(abind) # For abind
library(pROC)
library(data.table)

###Functions###
source("pp_funcs_v4.R")
sw.varswitch = "model_data/regionalized_model_data/Var-variable_switches.csv"
###Data pre-processing###
###import pre-processed data###
#A: 22,15,2
soc_eco = readRDS("model_data/regionalized_model_data/A-socioEcoDat-array.rds") 

#C: 15,16082
first_sight = readRDS("model_data/regionalized_model_data/C-firstSightings-matrix.rds") 
set.seed(123) # for reproducibility
first_sight <- first_sight[, sample(1:ncol(first_sight), 1500)]
#D: 15,15,1
pairwise_data = readRDS("model_data/regionalized_model_data/D-pairwiseData-array.rds")

#B (outdated): 6026,4
# species_traits = readRDS("model_data/B-speciesSighitngsTraits-taxa-df.rds")

#Tr:205,205,22 this is not used in this model; for now it only consider gdp, pop and dis
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


global_model = function(pars.start, sw.model = TRUE, sw.predict = TRUE, sw.maxit = 2000){
  ###Switches###
  sw.model = sw.model
  sw.maxit = sw.maxit
  sw.predict = sw.predict
  sw.varswitch = sw.varswitch
  
  ###Starting parameter###
  #other starting parameters will be added through the variable switch
  pars.start = pars.start
  
  ###Output###
  timestamp <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
  fdir <- tempfile(pattern = paste0("1-", timestamp, "-"), tmpdir = "output/grid_search/")
  #fdir <- "output/1-2023-09-28-10-41-23-a1fc510160be"
  

  

  
  
  ### Model Fitting###
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
    opt.out = optim(par = pars.start, le, echo = TRUE, ,control = list(maxit = sw.maxit))
    saveRDS(opt.out, fname)
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
    return(opt.results)
  }
  return(NA)
}
  

alpha_range <- c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001)
results <- data.frame(alpha=numeric(), auc = numeric(), dev = numeric(), ll = numeric(), convergence = logical())

for (alpha_init in alpha_range) {
  for (beta_init in beta_range) {
    # Use these initial values in the optim function
    print(paste0("alpha is ", alpha_init))
    print(paste0("beta is ", beta_init))
    pars.start = c(
      alpha = alpha_init
    )
    
    #Adding starting parameters from the varCSV(Variable Switch)
    if(length(which(varCSV$include))){
      pars.start = c(pars.start,
                     varCSV$starting.value[varCSV$include]
      )
      names(pars.start)[2:length(pars.start)] = varCSV$variable[varCSV$include]
    }
    
    #if we have species variables, then exclude alpha
    if(!all(is.na(sw.vars$sp.vars))){
      pars.start = pars.start[names(pars.start) %in% c(unlist(sw.vars), "pp.y0", "B")]
    }
    #Exclude B and pp.y0 since we don't need them if there are no PP variables
    if(all(is.na(sw.vars$pp.vars))){
      pars.start = pars.start[-which(names(pars.start) %in% c("pp.y0", "B"))]
    }
    
    
    
    ### Calculate the propagule pressure###
    t1 = Sys.time()
    pps = mk_ppstruct(sw.vars, first_sight, soc_eco, pairwise_data, trade)
    t2 = Sys.time()
    cat("Generating mk_ppstruct\n")
    print(t2-t1)
    model_res = global_model(pars.start,TRUE,TRUE,500)
    auc = model_res$result$auc
    dev = model_res$result$dev
    ll = model_res$result$ll
    convergence = model_res$result$converged
    final_alpha = as.numeric(model_res$result$par["alpha"])
    final_beta = as.numeric(model_res$result$par["B"])
    final_y0 = as.numeric(model_res$result$par["pp.y0"])
    final_gdp_src = as.numeric(model_res$result$par["gdp.src"])
    final_gdp_dst = as.numeric(model_res$result$par["gdp.dst"])
    final_pop_src = as.numeric(model_res$result$par["population.src"])
    final_pop_dst = as.numeric(model_res$result$par["population.dst"])
    final_dis = as.numeric(model_res$result$par["phys_dist"])
    final_t = as.numeric(model_res$result$par["t"])
    # Store the results
    results <- rbind(results, data.frame(alpha=alpha_init, beta=beta_init, auc = auc, dev = dev, ll = ll, convergence = convergence, final_alpha = final_alpha, final_beta = final_beta, final_y0 = final_y0, final_gdp_src = final_gdp_src, final_gdp_dst = final_gdp_dst, final_pop_src = final_pop_src, final_pop_dst = final_pop_dst, final_dis = final_dis, final_t = final_t))
    saveRDS(results, "result/grid_search_res/grid_search_1500_sp_500it_nm_v3.rds")
  }
}





