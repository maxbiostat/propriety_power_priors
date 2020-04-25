quantities <- function(x, alpha = .95){
  quants <- quantile(x, probs = c( (1-alpha)/2 , (1 + alpha)/2)) 
  out <- data.frame(
    lwr = quants[1],
    mean = mean(x, na.rm = TRUE),
    upr = quants[2]
  )
  return(out)
}
##
summarise_post <- function(obj){
  K <- length(obj)
  obj <- lapply(obj, as.matrix) ## taking care of one-dimensional "arrays"
  dims <- lapply(obj, ncol)
  par.names <- names(obj)
  all.names <-  vector(K, mode = "list")
  for(k in 1:K){
    all.names[[k]] <- paste(par.names[k], '[', 1:dims[[k]], ']', sep = "")
  }
  all.names.flat <- unlist(all.names)
  summys <- lapply(obj, function(par.dt) apply(par.dt, 2, quantities))
  res <- data.frame(
    parameter = all.names.flat,
    do.call(rbind, lapply(summys, function(y) do.call(rbind, y)))
  )
  row.names(res) <- NULL
  return(res)
}
##
get_deriv_and_mal <- function(a, stan.list, pars = NA, strict = FALSE){
  stan.list$a_0 <- a
  if(strict){
    power.prior <- sampling(compiled.model.prior, data = stan.list, iter = 5000,
                            control = list(adapt_delta = 0.99, max_treedepth = 15), refresh = 0)  
  }else{
    power.prior <- sampling(compiled.model.prior, data = stan.list, refresh = 0)
  }
  
  first_deriv <- mean(extract(power.prior, 'logL')$logL)
  second_deriv_part <- mean(extract(power.prior, 'logL_sq')$logL_sq)
  result <- list(
    chain = power.prior,
    lc = bridgesampling::bridge_sampler(power.prior, silent = TRUE)$logml,
    deriv_lc = first_deriv,
    second_deriv_lc = second_deriv_part - (first_deriv)^2
  )
  if(!is.na(pars[1])) result$summaries <- summarise_post(extract(power.prior, pars = pars))
  return(result)
}
## 
get_deriv_only <- function(a, stan.list, pars = NA, strict = FALSE){
  stan.list$a_0 <- a
  if(strict){
    power.prior <- sampling(compiled.model.prior, data = stan.list, iter = 5000,
                            control = list(adapt_delta = 0.99, max_treedepth = 15), refresh = 0)
  }else{
    power.prior <- sampling(compiled.model.prior, data = stan.list, refresh = 0)  
  }
  second_deriv_part <- mean(extract(power.prior, 'logL_sq')$logL_sq)
  first_deriv <- mean(extract(power.prior, 'logL')$logL)
  result <- list(
    chain = power.prior,
    deriv_lc = first_deriv,
    second_deriv_lc = second_deriv_part - (first_deriv^2)
  )
  if(!is.na(pars[1])) result$summaries <- summarise_post(extract(power.prior, pars = pars))
  return(result)
}
###
plug_gap <- function(x){
  sx <- sort(x)
  deltas <- diff(sx)
  pos <- which.max(deltas)
  new.x <- mean(sx[c(pos, pos + 1)])
  return(new.x)
}
###
build_grid <- function(eps = .01,  M = 10, J = 15, v1 = 5, v2 = 10, stan.list, pars, strict = FALSE){
  
  f0 <- get_deriv_and_mal(0.0, stan.list = stan.list, pars = pars, strict = strict) ## need to run for a_0 = 0 to get parameter summaries.
  
  ## Preliminary checking 
  fm <- get_deriv_and_mal(eps, stan.list = stan.list, pars = pars, strict = strict)
  fM <- get_deriv_and_mal(M, stan.list = stan.list, pars = pars, strict = strict)
  
  all.outs <- vector(J + 1, mode = "list")
  all.outs[[1]]$summaries <- f0$summaries
  all.outs[[2]]$summaries <- fm$summaries
  all.outs[[J + 1]]$summaries <- fM$summaries 

  
  same_sign <- sign(fm$deriv_lc) == sign(fM$deriv_lc)
  
  if(same_sign){ ## if they are the same sign, can construct a regularly spaced grid
    a0_grid <- seq(2*eps, M-eps, length.out = J - 2)
    all.outs[3:J] <- lapply(a0_grid, get_deriv_and_mal, stan.list = stan.list, pars = pars, strict = strict)
    
    inner.mals <- unlist(lapply(all.outs, function(x) x$lc ))
    inner.derivs <- unlist(lapply(all.outs, function(x) x$deriv_lc ))
    inner.2nd.derivs <- unlist(lapply(all.outs, function(x) x$second_deriv_lc ))
    
    res <- data.frame(
      a0 = c(0, eps, a0_grid, M),
      lc_a0 = c(f0$lc, ## should be zero. Maybe gives a gauge of how accurate estimation is.
                fm$lc, inner.mals, fM$lc),
      deriv_lc =  c(f0$deriv_lc, fm$deriv_lc, inner.derivs, fM$deriv_lc),
      second_deriv_lc = c(f0$second_deriv_lc, fm$second_deriv_lc, inner.2nd.derivs, fM$second_deriv_lc) 
    )
    
  }else{
    ## if signs change, then we need to be a bit more thoughtful
    ### Initialising
    a0s <- c(0, eps) 
    mals <- c(f0$lc, fm$lc)
    lderivs <- c(f0$deriv_lc, fm$deriv_lc)
    l2derivs <- c(f0$second_deriv_lc, fm$second_deriv_lc)
    
    ### Now start the "search"
    counter <- J - 2
    l <- eps
    u <- M
    z <- (u + l)/2
    d.old <- fm$deriv_lc
    delta <- 100 * eps
    while(counter > 0 && delta > v1*eps){
      lz <- get_deriv_and_mal(z, stan.list = stan.list, pars = pars, strict = strict)
      a0s <- c(a0s, z)
      mals <- c(mals, lz$lc)
      lderivs <- c(lderivs, lz$deriv_lc)
      l2derivs <- c(l2derivs, lz$second_deriv_lc)
      all.outs[[3 + (J-2)-counter]]$summaries <- lz$summaries
      
      same_sign_below <- sign(d.old) == sign(lz$deriv_lc)
      if(same_sign_below){
        l <- z
        u <- u
      }else{ ## then it must be the case that sign(du) == sign(lz$deriv_lc)
        l <- l
        u <- z
      }
      z.old <- z
      z <- (u + l)/2
      delta <- abs(z - z.old)
      counter <- counter - 1
    }
    if(counter > 1){## in case it stopped via getting close enough to the zero
      lz <- get_deriv_and_mal(z, stan.list = stan.list, pars = pars, strict = strict)
      a0s <- c(a0s, z)
      mals <- c(mals, lz$lc)
      lderivs <- c(lderivs, lz$deriv_lc)
      l2derivs <- c(l2derivs, lz$second_deriv_lc)
      all.outs[[3 + (J-2)-counter]]$summaries <- lz$summaries
      
      counter <- counter - 1
      newrange <- c(max(0, z - v2*eps), min(z + v2 * eps, M))
      new_a0s <- c( a0s[a0s > newrange[1] & a0s < newrange[2]], newrange)
      
      while(counter > 0){
        ## above 
        z <-  plug_gap(new_a0s)
        new_a0s <- c(new_a0s, z)
        lz <- get_deriv_and_mal(z, stan.list = stan.list, pars = pars, strict = strict)
        lderivs <- c(lderivs, lz$deriv_lc)
        l2derivs <- c(l2derivs, lz$second_deriv_lc)
        a0s <- c(a0s, z)
        mals <- c(mals, lz$lc)
        all.outs[[3 + (J-2)-counter]]$summaries <- lz$summaries
        
        counter <- counter - 1
      }
    }
    
    a0s <- c(a0s, M)
    mals <- c(mals, fM$lc)
    lderivs <- c(lderivs, fM$deriv_lc)
    l2derivs <- c(l2derivs, fM$second_deriv_lc)
    
    res <- data.frame(a0 = a0s, lc_a0 = mals, deriv_lc = lderivs, second_deriv_lc = l2derivs)
  }
  
  if(!is.na(pars[1])) parameter.summaries <- lapply(all.outs, function(x) x$summaries)
  
  if(!is.na(pars[1])){
    out <- list(
      result = res,
      summaries = parameter.summaries
    )
  }else{
    out <- list(
      result = res,
    )
  }
  return(out)
}
#
build_grid_derivOnly <- function(eps = .01,  M = 10, J = 15, v1 = 5, v2 = 10, stan.list, pars, strict = FALSE){
  
  f0 <- get_deriv_only(0.0, stan.list = stan.list, pars = pars, strict = strict) ## need to run for a_0 = 0 to get parameter summaries.
  
  ## Preliminary checking 
  fm <- get_deriv_only(eps, stan.list = stan.list, pars = pars, strict = strict)
  fM <- get_deriv_only(M, stan.list = stan.list, pars = pars, strict = strict)
  
  all.outs <- vector(J + 1, mode = "list")
  all.outs[1]$summaries <- f0$summaries
  all.outs[2]$summaries <- fm$summaries
  all.outs[J + 1]$summaries <- fM$summaries
  
  same_sign <- sign(fm$deriv_lc) == sign(fM$deriv_lc)
  
  if(same_sign){ ## if they are the same sign, can construct a regularly spaced grid
    a0_grid <- seq(2*eps, M-eps, length.out = J - 2)
    all.outs[3:J] <- lapply(a0_grid, get_deriv_only, stan.list = stan.list, pars = pars, strict = strict)
    
    inner.derivs <- unlist(lapply(all.outs, function(x) x$deriv_lc ))
    inner.2nd.derivs <- unlist(lapply(all.outs, function(x) x$second_deriv_lc ))
    
    res <- data.frame(
      a0 = c(0, eps, a0_grid, M),
      deriv_lc = c(f0$deriv_lc, fm$deriv_lc, inner.derivs, fM$deriv_lc),
      second_deriv_lc = c(f0$second_deriv_lc, fm$second_deriv_lc, inner.2nd.derivs, fM$second_deriv_lc)
    )
    
  }else{
    ## if signs change, then we need to be a bit more thoughtful
    ### Initialising
    
    a0s <- c(0, eps) 
    lderivs <- c(f0$deriv_lc, fm$deriv_lc)
    l2derivs <- c(f0$second_deriv_lc, fm$second_deriv_lc)
    
    ### Now start the "search"
    counter <- J - 2
    l <- eps
    u <- M
    z <- (u + l)/2
    d.old <- fm$deriv_lc
    delta <- 100 * eps
    while(counter > 0 && delta > v1*eps){
      lz <- get_deriv_only(z, stan.list = stan.list, pars = pars)
      a0s <- c(a0s, z)
      lderivs <- c(lderivs, lz$deriv_lc)
      l2derivs <- c(l2derivs, lz$second_deriv_lc)
      all.outs[[3 + (J-2)-counter]]$summaries <- lz$summaries
      
      same_sign_below <- sign(d.old) == sign(lz$deriv_lc)
      if(same_sign_below){
        l <- z
        u <- u
      }else{ ## then it must be the case that sign(du) == sign(lz$deriv_lc)
        l <- l
        u <- z
      }
      z.old <- z
      z <- (u + l)/2
      delta <- abs(z - z.old)
      counter <- counter - 1
    }
    if(counter > 1){## in case it stopped via getting close enough to the zero
      lz <- get_deriv_only(z, stan.list = stan.list, pars = pars, strict = strict)
      a0s <- c(a0s, z)
      lderivs <- c(lderivs, lz$deriv_lc)
      l2derivs <- c(l2derivs, lz$second_deriv_lc)
      all.outs[[3 + (J-2)-counter]]$summaries <- lz$summaries
      
      counter <- counter - 1
      newrange <- c(max(0, z - v2*eps), min(z + v2 * eps, M))
      new_a0s <- c( a0s[a0s > newrange[1] & a0s < newrange[2]], newrange)
      while(counter > 0){
        ## above 
        z <-  plug_gap(new_a0s)
        new_a0s <- c(new_a0s, z)
        lz <- get_deriv_only(z, stan.list = stan.list, pars = pars, strict = strict)
        lderivs <- c(lderivs, lz$deriv_lc)
        l2derivs <- c(l2derivs, lz$second_deriv_lc)
        all.outs[[3 + (J - 2) - counter]]$summaries <- lz$summaries
        a0s <- c(a0s, z)
        
        counter <- counter - 1
      }
    }
    
    a0s <- c(a0s, M)
    lderivs <- c(lderivs, fM$deriv_lc)
    l2derivs <- c(l2derivs, fM$second_deriv_lc)
    
    res <- data.frame(a0 = a0s, deriv_lc = lderivs, second_deriv_lc = l2derivs)
  }
  
  if(!is.na(pars[1])) summaries <- lapply(all.outs, function(x) x$summaries)
  
  if(!is.na(pars[1])){
    out <- list(
      result = res,
      summaries = summaries
    )
  }else{
    out <- list(
      result = res
    )
  }
  return(out)
  
}
#
create_lc_df <- function(a0_grid, stan.list, pars = NA, strict = FALSE){
  a0_grid <- unique(c(0, a0_grid)) # adds zero and avoids double computing if it's already there
  all.outs <- lapply(a0_grid, get_deriv_and_mal, stan.list = stan.list, pars = pars, strict = strict)
  mals <- unlist(lapply(all.outs, function(x) x$lc ))
  derivs <- unlist(lapply(all.outs, function(x) x$deriv_lc ))
  second.derivs <- unlist(lapply(all.outs, function(x) x$second_deriv_lc ))
  if(!is.na(pars[1])) summaries <- lapply(all.outs, function(x) x$summaries)
  
  res <- data.frame(
    a0 = a0_grid,
    lc_a0 =  mals,
    deriv_lc = derivs,
    second_deriv_lc = second.derivs
  )
  
  if(!is.na(pars[1])){
    out <- list(
      result = res,
      summaries = summaries
    )
  }else{
    out <- list(
      result = res
    )
  }
  return(out)
}