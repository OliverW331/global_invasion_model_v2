#This version of pp_funcs change the formula of probability. Time and N are not independent.

# Deviance explained, using mean observations as the null predictions
#'@pred = vector predicted values
#'@obs = vector observed values
dev.expl = function(pred, obs, null = NA){
  if(all(is.na(null))){
    null = mean(obs)
  }
  ll = sum(log(1-abs(obs - pred)))
  nl = sum(log(1-abs(obs - null)))
  return((nl-ll)/nl)
}

### mk_ppstruct to generate propagule pressure ####
#elements of PP remain constant, except for parameters. So we can percalcuate everything?
#for each variable determining PP, we want to keep track of the relevant elements of data structures so we can multiply them along with par later
#one type of array is by yrs, src, dst, spscies
## Now modified to only consider variables in sw.var$pp.vars

fill_pp_dim <- function(x, dim, reord, pp){
  #get a template the right dimensions
  pp1 = pp 
  #copies the changing dimensions a number of times
  #e.g.: yrs, dst, copies, src # times
  tmp = array(as.vector(x), dim = dim)
  #reorders dimensions to match pp (yrs, src, dst, species)
  tmp = aperm(tmp, reord)
  #fills in pp, and replicates tmp for any higher dimensions (e.g., species)
  pp1[] = tmp 
  return(pp1)
}

process_fst_sight <- function(ns, first_sight, pp){
  for (s in 1:ns){
    #x is the countries
    x = first_sight[,s]
    mn = min(x, na.rm = T)
    
    #iterate through each destination country
    for (i in 1:length(x)){
      #countries have to be in consistent order, get sources, and each yr contribute to dst i
      #if not invaded, all invaded are sources
      #src is the list of index of source countries
      if (is.na(x[i])){
        src = which(!is.na(x))
        mx = tot_yrs
      }else{
        src = which(x<x[i])
        mx = x[i]
      }
      
      #don't include th first invasion in Europe
      if (length(src)>0){
        #determine which src/yr combinations can provide propagules
        #for each year, need to sum across all source countries contribute to a given destination
        #for each source, put 0 before inv, up to the date of invasion: generates a yr by src matrix
        l1 <- lapply(x[src], function(x2){c(rep(NA, x2-mn), rep(1, mx-x2+1))})
        #need to have the years, associated with location
        pp[mn:mx, src, i, s] = unlist(l1)
        
      }else{
        #this country/species has no sources (it is the first invasion in Europe)
        pp[mn:mx, src, i, s] = NA
      }
    }
  }
  return(pp)
}

process_soc_eco <- function(sw.vars, soc_eco, pps, pp){
  soc_eco_var = dimnames(soc_eco)[[3]]
  soc_eco_var = soc_eco_var[soc_eco_var %in% gsub("(\\.src)|(\\.dst)","",sw.vars$pp.vars)]
  if (length(soc_eco_var) > 0){
    for (var in soc_eco_var){
      #src
      pps[[paste0(var, ".src")]] = pp * as.vector(soc_eco[,,var])
      #dst
      pps[[paste0(var, ".dst")]] = fill_pp_dim(soc_eco[,,var], c(tot_yrs, nc, nc), c(1,3,2), pp)
    }
  }
  return(pps)
}

process_pairwise <- function(sw.vars, pairwise_data, pps, pp){
  pairwise_var = dimnames(pairwise_data)[[3]]
  pairwise_var = pairwise_var[pairwise_var %in% sw.vars$pp.vars]
  if(length(pairwise_var) > 0){
    pps[pairwise_var] = lapply(pairwise_var, FUN = function(x){
      return(fill_pp_dim(pairwise_data[,,x],c(nc,nc,tot_yrs),c(3,1,2), pp))
    })
  }
  return(pps)
}

process_trade <- function(sw.vars, trade, pps, pp){
  if("bilateral_trade" %in% sw.vars$pp.vars){
    pps[['bilateral_trade']]=fill_pp_dim(trade, c(nc,nc,tot_yrs), c(3,1,2), pp)
  }
  return(pps)
}

add_remove <- function(pps, pp){
  pps[['remove']]=apply(pp,c(1,3,4),function(x){
    if(all(is.na(x)))
      return(NA)
    return(1)
  })
  return(pps)
}

mk_ppstruct = function(sw.vars, first_sight, soc_eco, pairwise_data, trade){
  #mk_ppstruct:
  # used in generation of pps
  # destination/yr predictors. Need to keep order yrs, src, dst, sp, so need aperm
  # the vector, the dimensions up to the one being filled, the reordered dimensions needed to match pp
  
  ### the basic parameters for dimension
  nc = nrow(first_sight)
  
  #number of species
  ns = ncol(first_sight)
  
  #number of years
  last_yr = max(first_sight, na.rm=T)
  first_yr = min(first_sight, na.rm=T)
  tot_yrs = last_yr - first_yr + 1
  
  #Question:
  #given that structure, how would one multiply by distance?
  #probably an each function?
  
  pp = array(NA, dim = c(tot_yrs, nc, nc, ns))
  
  #we want to keep tract of each relevant variable in l
  #we need to also keep track of how they add or multiply together
  #the dimensions of the matrix in pps are all (yr, src, dst, sp)
  
  pps = list()
  
  #process the matrix of first sight, separating by species, and taking the years each country was invaded
  #only consider the time periods that are not invaded yet
  #for a country that is not invaded: 
  # other countries not invaded throughout time = NA; 
  # all invaded countries are sources starting from the time they are invaded = 1; other time period = NA;
  #for a country that is invaded:
  # other countries not invaded throughout time = NA;
  # the countries that are invaded before it are sources from the time they are invaded to the time first sight showed up in this country = 1;
  # other period = NA;
  
  pp = process_fst_sight(ns, first_sight, pp)
  
  # Generate all variables in the socio_eco data, for src and dst, and multiply it with pp;
  # pps contains the result of multiplications with the names of corresponding variables;
  # e.g.: "gdp.src", "population.dst";
  # 
  pps = process_soc_eco(sw.vars, soc_eco, pps, pp)
  
  # Generate all variables in the pairwise data
  # e.g.: "phys_dist"
  pps = process_pairwise(sw.vars, pairwise_data, pps, pp)
  
  # Generate all variables in the trade data
  # e.g.: "bilateral_trade"
  pps = process_trade(sw.vars, trade, pps, pp)
  
  #remove: shows whether there were invasions happened in this destination country in each year (1: Yes; NA: No)
  pps = add_remove(pps, pp)
  
  gc()
  
  return(pps)
}

#Calculate prop pressure for each species, location, and each year
#go through each species
#return a single value: summed propagule pressure for each year yr, summed across all sources
pp_func <- function(start_par){
  #iterate through each country, find previous invasions
  #each element (species) has array of country by year
  #   pps - from mk_ppstruct - dim:tot_yrs,nsrc,ndst,nspecies #is in this dimensionality 
  # The propagule pressure to each destination needs to be aggregated across all sources, so the new dimensionality will be yr, dst country, and species, from a 4 dimensional array
  # for a gravity formulation, values in pps would already be log transformed, so would get sum(X1^b1X2^b2): i.e., sum logs, take exponential
  # This doesn't include the intercept though!
  agg1 = lapply(sw.vars$pp.vars, FUN = function(vv){
    return(start_par[vv]*pps[[vv]])
  })
  ## Agg1 is a list of arrays that we can sum using Reduce
  # add up the values through the variables; the dimension of tmp is yr,nsrc,ndst,sp
  tmp = Reduce('+', agg1)
  # Also include intercept
  tmp = tmp + start_par["pp.y0"]
  
  # NOTE: RETURNS NA IF ANY VALUE IS ALSO NA
  # test=par["pp.y0"]+par["gdp.src"]*pps[["gdp.src"]]+par["gdp.dst"]*pps[["gdp.dst"]]+par["pop"]*pps[["pop"]]+par["dist"]*pps[["dist"]] #if this could be calculated for each src,dst pair, for each year, and then summed separately from all sources, this could be N.
  # convert back to N
  tmp=exp(tmp)
  #sum across sources (i.e., by yr, dest, and species))
  #NOTE: May also want a "c" coefficient - i.e., the shape parameter in a weibull. It would be applied to the entire sum in a given year
  return(apply(tmp,c(1,3,4),sum,na.rm=T))
  #we are asking the question what is the functional relationship between prop pressure and these macroecological variables. 
}

#calculate the probability
get_p <- function(x){
  yoi = x[length(x)]
  
  # Remove the final column
  x = x[-length(x)]
  if(all(is.na(x))){ # If first invasion, don't included
    return(NA)
  }
  
  #it is uninvaded, return prob of remaining uninvaded the entire time
  if(is.na(yoi)){
    return(sum(x,na.rm=T))
  }
  #it has some non-NA values, so need to get rid of entry from C, held in the last position
  pi=log(1-exp(x[yoi]))
  pu=sum(x[-yoi],na.rm=T)
  
  #return back prob i in yr invasion, and prob uninvaded up til then.
  #need to figure out which ones were uninvaded
  return(pi+pu)
}

# Returns log probabilities array of country, species, year
# Extension of get_p, and identical to the verson in '11_estabFuncs.R' with minor changes
get_p_time = function(x){
  yoi = x[length(x)]
  
  # Remove the final column
  x = x[-length(x)]
  if(all(is.na(x))){ # If first invasion, don't included
    x[]=NA
    return(x)
  }
  
  # Never invaded, return pu (which is just probability of not invading the whole 60 years)
  if(is.na(yoi)) return(x)
  # If invaded, return pu from 1:(yoi-1) and pi at yoi, then NA for the rest
  # Therefore, we maximize (log) probability of not establishing up to yoi,
  # and the probability of establishing at yoi (because use multiply by -1)
  x[yoi] = log(1-exp(x[yoi]))
  
  # If yoi is max year, then we don't have to do this
  # This should already happen, but we'll do it again
  if(yoi != length(x)) x[(yoi+1):length(x)] = NA
  
  return(x)
}

predict.le <- function(start_par, time = TRUE, echo = FALSE){
  if(echo){
    print(Sys.time())
    print(start_par)
  }
  
  #pp_func returns propagule pressure array of dim [year,country,species]
  # Currently in order: tot_yrs, nc, ns
  if(all(is.na(sw.vars$pp.vars))){
    # If we have no propagule pressure variables, then we will exclude the pp_func step (since it is long)
    # Instead we just add an array of all 1s
    N = array(1, dim = dim(pps[["remove"]]))
    start_par["B"] = 0
  }else{
    N = pp_func(start_par)
  }
  # 'remove' here will turn all time-steps after invasion to NA
  # To deal with the invaded-at-t1 issue
  # Test cases (country 1): species 1 = never invaded, species 2 = invaded at t1, species 220 = invaded at t10
  # All 0s -> throws error
  N=N*pps[["remove"]]
  
  N = aperm(N, c(2,3,1)) #needs to be in order nc, ns, tot_yrs (plus one, for year of first record), to multiply by coefficient a
  
  #Make C + tmp in the appropriate format for N? N is row countries, col years, array species. C and tmp are row country by species.
  # If there are no variables provided, will instead fit using a fixed parameter (alpha)
  if(all(is.na(sw.vars$sp.vars))){
    a = start_par["alpha"] # single scalar value
  }else{
    a = as.vector(apply(B, 1, get_a, par=start_par)) #returns back a vector - order same as matrix nrow:country, ncol:species 
  }
  
  if(echo){
    print(paste("min(alpha) =", min(a)))
    print(paste("min(N) =", min(N, na.rm = T)))
  }
  
  # If EITHER a or PP are negative, then skip
  if(min(a) < 0 | min(N, na.rm = T) < 0){
    return(NULL)
  }
  
  # U gives log(prob of remaining uninvaded), then need to sum over time (for each species location combo).
  # Therefore, U must be negative
  #we add a year to tot_yrs, cause we need to pass additional info about whether it is invaded or not
  U=array(NA, dim=c(nc,ns,tot_yrs+1))
  #### A LOT OF 0s WHICH WILL THROW AN ERROR IN GET_P, PRODUCTING INFINITES
  if(start_par["B"] != 0){
    U[,,-(tot_yrs+1)] =-a*(N)^start_par["B"]
  }else{
    U[,,-(tot_yrs+1)] =-a*aperm(pps[["remove"]], c(2,3,1))
    
  }
  #if par["B"] == 0, even NAs in N will be treated as 1.
  
  # If we want time as a predictor
  if("t" %in% sw.vars$time.vars){
    # Add time as predictor to U
    # We subtract because technically the whole thing needs to be negative (see working doc)
    U[,,-(tot_yrs+1)] = U[,,-(tot_yrs+1)] * time_array ^ start_par["t"]
  }
  
  if(max(U, na.rm = T) > 0) return(NULL)
  
  #vector with country (row) by species (column).
  #for each U, need to extract the year of invasion, for each species and country. Info is in first_sight
  U[,,(tot_yrs+1)] = first_sight #C is country by species
  #U[,,-(tot_yrs+1)] = 1-exp(U[,,-(tot_yrs+1)])
  # returns back a matrix of log probabilities, with row:countries, col:species
  # return(apply(U,c(1,2),get_p))
  # returns array of log probabilities, country, species and year
  if(time){
    return(apply(U,c(1,2), get_p_time))
  }else{
    return(apply(U,c(1,2), get_p))
  }
}

le <- function(par, echo = FALSE, time = TRUE){
  if(echo) cat("FITTING\n")
  
  if(all(is.na(sw.vars$pp.vars))){
    par["B"] = 0
  }
  
  # B cannot be negative!
  if(par["B"] < 0){
    return(10^24)
  }
  
  # returns back a matrix of log probabilities, with row:countries, col:species
  P = predict.le(par, time = time, echo = echo)
  
  # Basically remove if outside of contraints
  if(all(is.null(P))){
    return(10^24)
  }
  
  #then we sum to get the joint likelihood
  # NA are excluded variables, which replace 0 values in earlier versions
  P = sum(P, na.rm = TRUE)
  
  if(echo){
    print(paste("-P =", -P))
    cat("DONE\n\n")
  }
  
  return(-P) #return neg, because optim minimizes (large negative P => low prob, so make it large positive)
}

