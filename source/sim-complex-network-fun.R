
generate.get.Jf_y <- function(ct.lst,
                           rct.lst,
                           p,
                           jsens,
                           envir = .GlobalEnv,
                           flickering = FALSE){

  if(flickering){
    rct.lst <- c(rct.lst,
                 "F->O",
                 "O->F")
  }

  ########### fx
  V <- get.V(ct.lst = ct.lst, rct.lst = rct.lst)
  Vth <- matrix(data = "0", nrow(V), nrow(V))
  rownames(Vth) <- colnames(Vth) <- ct.lst
  for (ct in ct.lst) {
    for (r in names(V[ct,V[ct,] != 0])) {
      ct.rgt <- as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))[1]
      V_ct_r <- eval(parse(text=paste("V['", ct, "','", r, "']", sep = "")))
      if(Vth[ct, ct.rgt] == "0"){
        Vth[ct, ct.rgt] <- paste(ifelse(V_ct_r == 1, "", "-"),"lambda['", r, "']", sep = "")
      }else{
        Vth[ct, ct.rgt] <- paste(Vth[ct, ct.rgt], ifelse(V_ct_r == 1, " + ", " -"), paste("lambda['", r, "']", sep = ""), sep = "")
      }
    }
  }

  bfx <- Vth[-ncol(Vth),ncol(Vth)]
  Vth <- Vth[-ncol(Vth),]

  Jfx.x <- Vth <- t(sapply(1:nrow(Vth), function(i){
    sapply(1:(ncol(Vth)-1), function(j){
      paste(ifelse(Vth[i,j] == "0", ifelse(Vth[i,ncol(Vth)] == "0", "0", ""), Vth[i,j]),
            ifelse(Vth[i,ncol(Vth)] == "0", "", Vth[i,ncol(Vth)]),
            sep = ifelse(Vth[i,ncol(Vth)] == "0", "", " -"))
    })
  }))

  ############ fp

  wPidx <- which(lower.tri(matrix(data = NA,
                                  nrow = length(ct.lst) - 1,
                                  ncol = length(ct.lst) - 1), diag = TRUE), arr.ind = TRUE)

  Jq.string <- matrix(data = "0",
                      nrow = nrow(wPidx),
                      ncol = length(ct.lst) - 1)
  rownames(Jq.string) <- paste0("p",1:nrow(wPidx))
  colnames(Jq.string) <- head(ct.lst,-1)



  for (r in 1:nrow(wPidx)) {
    i <- ct.lst[wPidx[r,1]]
    j <- ct.lst[wPidx[r,2]]

    wri <- sapply(rct.lst, function(rct){i %in% unlist(strsplit(x = rct,
                                                                split = "->",
                                                                fixed = TRUE))})

    wrj <- sapply(rct.lst, function(rct){j %in% unlist(strsplit(x = rct,
                                                                split = "->",
                                                                fixed = TRUE))})
    if(length(rct.lst[wri & wrj]) > 0){
      rgts <- sapply(rct.lst[wri & wrj], function(rct){unlist(strsplit(x = rct,
                                                                       split = "->",
                                                                       fixed = TRUE))[1]})
      for (c in 1:length(rgts)) {
        rgt <- rgts[c]
        if(rgt != tail(ct.lst,1)){
          if(Jq.string[r,rgt] == "0"){
            Jq.string[r,rgt] <- ifelse(i == j, paste0("lambda['", names(rgt),"']"), paste0("-lambda['", names(rgt),"']"))
          }else{
            Jq.string[r,rgt] <- ifelse(i == j,
                                       paste0(Jq.string[r,rgt], " + lambda['", names(rgt),"']"),
                                       paste0(Jq.string[r,rgt], " -lambda['", names(rgt),"']"))
          }
        }else{
          Jq.string[r,] <- sapply(Jq.string[r,], function(sr){
            if(sr == "0"){
              return(paste0("-lambda['", names(rgt),"']"))
            }else{
              return(paste0(sr, " -lambda['", names(rgt),"']"))
            }
          })
        }
      }

    }

  }

  Jfp.x <- Jq.string

  bfp.string <- matrix(data = "0",
                       nrow = nrow(wPidx),
                       ncol = 1)
  rownames(bfp.string) <- paste0("p",1:nrow(wPidx))
  colnames(bfp.string) <- tail(ct.lst, 1)

  for (r in 1:nrow(wPidx)) {
    i <- ct.lst[wPidx[r,1]]
    j <- ct.lst[wPidx[r,2]]

    wri <- sapply(rct.lst, function(rct){i %in% unlist(strsplit(x = rct,
                                                                split = "->",
                                                                fixed = TRUE))})

    wrj <- sapply(rct.lst, function(rct){j %in% unlist(strsplit(x = rct,
                                                                split = "->",
                                                                fixed = TRUE))})
    if(length(rct.lst[wri & wrj]) > 0){
      rgts <- sapply(rct.lst[wri & wrj], function(rct){unlist(strsplit(x = rct,
                                                                       split = "->",
                                                                       fixed = TRUE))[1]})
      for (c in 1:length(rgts)) {
        rgt <- rgts[c]
        if(rgt == tail(ct.lst,1)){
          bfp.string[r,1] <- paste0("lambda['", names(rgt),"']")
        }
      }
    }
  }

  bfp <- bfp.string

  ######################
  Vth.forP_1 <- matrix(data = "0",
                       nrow = nrow(Vth)*(nrow(Vth) - 1)/2 + nrow(Vth),
                       ncol = nrow(Vth)*(nrow(Vth) - 1)/2 + nrow(Vth))

  Vth.forP_1[1:nrow(Vth),
             1:nrow(Vth)] <- Vth

  wPidx <- which(lower.tri(matrix(data = NA,
                                  nrow = length(ct.lst) - 1,
                                  ncol = length(ct.lst) - 1), diag = TRUE), arr.ind = TRUE)

  idxs <- matrix(data = NA,
                 nrow = length(ct.lst) - 1,
                 ncol = length(ct.lst) - 1)
  idxs[rbind(wPidx,
             wPidx[,c(2,1)])] <- 1:(nrow(Vth)*(nrow(Vth) - 1)/2 + nrow(Vth))

  nblocks <- nrow(Vth) - 1
  bsizes <- (nrow(Vth) - 1):1
  r <- nrow(Vth) + 1
  for (i in 1:nblocks) {

    Vth.forP_1[r:(r+bsizes[i] - 1),idxs[,i + 1]] <- Vth[-(1:i),]

    r <- r + bsizes[i]
  }

  ######################
  Vth.forP_2 <- matrix(data = "0",
                       nrow = nrow(Vth)*(nrow(Vth) - 1)/2 + nrow(Vth),
                       ncol = nrow(Vth)*(nrow(Vth) - 1)/2 + nrow(Vth))

  nblocks <- nrow(Vth)
  bsizes <- (nrow(Vth)):1
  r <- 1
  for (i in 1:nblocks) {

    Vth.forP_2[cbind(rep(r:(r+bsizes[i] - 1), each = nblocks),
                     c(idxs[,i:nrow(idxs)]))] <- rep(Vth[i,], bsizes[i])

    r <- r + bsizes[i]
  }

  Vth.forP_1_plus_Vth.forP_2 <- Vth.forP_1

  for (i in 1:nrow(Vth.forP_1_plus_Vth.forP_2)) {
    for (j in 1:nrow(Vth.forP_1_plus_Vth.forP_2)) {
      if(Vth.forP_2[i,j] != "0"){
        if(Vth.forP_1_plus_Vth.forP_2[i,j] != "0"){
          Vth.forP_1_plus_Vth.forP_2[i,j] <- paste(Vth.forP_1_plus_Vth.forP_2[i,j], Vth.forP_2[i,j], sep = ifelse(substring(Vth.forP_2[i,j], 1, 1) == "-", " ", " + "))
        }else{
          Vth.forP_1_plus_Vth.forP_2[i,j] <- Vth.forP_2[i,j]
        }
      }
    }
  }

  Jfp.p <- Vth.forP_1_plus_Vth.forP_2

  Jfy.y <- rbind(cbind(Jfx.x, matrix("0",
                                     nrow = nrow(V) - 1,
                                     ncol = (nrow(V) - 1)*(nrow(V) - 2)/2 + (nrow(V) - 1))),
                 cbind(Jfp.x, Jfp.p))



  get.Jfy.y.string <- paste("get.Jfy.y <- function(lambda){
  matrix(data = c(",
                              paste(as.vector(unlist(lapply(split(as.vector(t(Jfy.y)), ceiling(seq_along(1:(nrow(Jfy.y)^2))/nrow(Jfy.y))), function(lst){
                                paste(lst, collapse = ",")
                              }))), collapse = ",\n"),
                              "),\n byrow = TRUE, nrow = ", nrow(Jfy.y), ", ncol = ", nrow(Jfy.y),
                              ")\n}", sep = "")

  eval(parse(text=get.Jfy.y.string), envir = envir)

  get.Jfy.y.forgrad.string <- paste("get.Jfy.y.forgrad <- function(",
                                    paste0(paste(paste0("p", 1:p), collapse = ",\n"), ",\n"),
                                    "t){\n",
                                    paste0("theta <- c(",
                                           paste(paste0("p", 1:p), collapse = ",\n"),
                                           ")\n"),
                                    "lambda <- get.lambda(theta = theta,\n V = get.voltage(t))\n",
                                    "c(",
                                    paste(as.vector(unlist(lapply(split(as.vector(t(Jfy.y)), ceiling(seq_along(1:(nrow(Jfy.y)^2))/nrow(Jfy.y))), function(lst){
                                      paste(lst, collapse = ",")
                                    }))), collapse = ",\n"),
                                    ")\n}", sep = "")

  eval(parse(text=get.Jfy.y.forgrad.string), envir = environment())

  get.Jfy.y.dj_c <- Deriv(f = get.Jfy.y.forgrad,
                      x = c(paste0("p", jsens)),
                      combine = "c",
                      cache.exp = FALSE)

  get.Jfy.y.dj_c_deparsed <- deparse(get.Jfy.y.dj_c)
  get.Jfy.y.dj_c_deparsed[1] <- paste0("get.Jfy.y.dj_c <- ", get.Jfy.y.dj_c_deparsed[1])
  eval(parse(text=paste0(get.Jfy.y.dj_c_deparsed, collapse = " \n")), envir = envir)

  eval(parse(text=paste0(paste0("get.Jfy.y.dj <- function(t, \n parms){\n matrix(data = get.Jfy.y.dj_c(",
                                paste(paste0("p",1:p, "=parms[", 1:p,"]"), collapse = ',\n'),
                                ",\n",
                                "t = t), nrow = ", nrow(Jfy.y)*length(jsens), ", ncol = ",  nrow(Jfy.y), ", byrow = TRUE)\n}"), collapse = " \n")),
       envir = envir)

  get.bfy.string <- paste("get.bfy <- function(lambda){\nc(",
                          paste(paste(bfx, collapse = ",\n"),
                                paste(bfp, collapse = ",\n"), sep = ",\n"),
                            ")\n}", sep = "")

  eval(parse(text=get.bfy.string), envir = envir)

  get.bfy.forgrad.string <- paste("get.bfy.forgrad <- function(",
                                  paste0(paste(paste0("p", 1:p), collapse = ",\n"), ",\n"),
                                  "t){\n",
                                  "lambda <- get.lambda(theta = c(",
                                  paste(paste0("p", 1:p), collapse = ",\n"),
                                  "),\n V = get.voltage(t))\n\n",
                                  "c(",
                                  paste(paste(bfx, collapse = ",\n"),
                                        paste(bfp, collapse = ",\n"), sep = ",\n"),
                                  ")\n}", sep = "")

  eval(parse(text=get.bfy.forgrad.string), envir = environment())

  get.bfy.dj_c <- Deriv(f = get.bfy.forgrad,
                        x = c(paste0("p", jsens)),
                        combine = "c",
                        cache.exp = TRUE)

  get.bfy.dj_c_deparsed <- deparse(get.bfy.dj_c)
  get.bfy.dj_c_deparsed[1] <- paste0("get.bfy.dj_c <- ", get.bfy.dj_c_deparsed[1])
  eval(parse(text=paste0(get.bfy.dj_c_deparsed, collapse = " \n")), envir = envir)

  eval(parse(text=paste0(paste0("get.bfy.dj <- function(t, \n parms){\n get.bfy.dj_c(",
                                paste(paste0("p",1:p, "=parms[", 1:p,"]"), collapse = ',\n'),
                                ",\n",
                                "t = t)\n}"), collapse = " \n")),
       envir = envir)
}


generate.get.x.steady_djs <- function(p,
                                      jsens,
                                      ct.lst,
                                      rct.lst,
                                      envir = .GlobalEnv,
                                      flickering = FALSE){
  
  generate.get.x.steady.forgrad(p = p,
                                ct.lst = ct.lst,
                                rct.lst = rct.lst,
                                envir = environment(),
                                flickering = flickering)
  
  get.x.steady_djs_c <- Deriv(f = get.x.steady.forgrad,
                              x = c(paste0("p", jsens)),
                              combine = "c",
                              cache.exp = FALSE)
  
  get.x.steady_djs_c_deparsed <- deparse(get.x.steady_djs_c)
  get.x.steady_djs_c_deparsed[1] <- paste0("get.x.steady_djs_c <- ", get.x.steady_djs_c_deparsed[1])
  eval(parse(text=paste0(get.x.steady_djs_c_deparsed, collapse = " \n")), envir = envir)
  
  eval(parse(text=paste0(paste0("get.x.steady_djs <- function(t,
                   parms){\nget.x.steady_djs_c(",
                                paste(paste0("p",1:p, "=parms[", 1:p,"]"), collapse = ',\n'),
                                ",\n",
                                "t = t)\n}"), collapse = " \n")),
       envir = envir)
}

generate.get.dx_dj_di <- function(p,
                                  n,
                                  jsens,
                                  ct.lst,
                                  rct.lst,
                                  envir = .GlobalEnv){
  
  generate.get.dx.forgrad(p = p,
                          n = n,
                          ct.lst = ct.lst,
                          rct.lst = rct.lst,
                          envir = environment())
  
  get.dx_djs <- Deriv(f = get.dx.forgrad,
                      x = c(paste0("p", jsens)),
                      combine = "c",
                      cache.exp = FALSE)
  
  get.dx_djs_dxis_cbind <- Deriv(f = get.dx_djs,
                                 x = paste0("x", 1:n),
                                 combine = "cbind",
                                 cache.exp = FALSE)
  
  get.dx_djs_dxis_cbind_deparsed <- deparse(get.dx_djs_dxis_cbind)
  get.dx_djs_dxis_cbind_deparsed[1] <- paste0("get.dx_djs_dxis_cbind <- ", get.dx_djs_dxis_cbind_deparsed[1])
  eval(parse(text=paste0(get.dx_djs_dxis_cbind_deparsed, collapse = " \n")), envir = envir)
  
  
  eval(parse(text=paste0(paste0("get.dx_dj_di <- function(t,
                   y,
                   parms){\nreturn(get.dx_djs_dxis_cbind(",
                                paste(paste0("p",1:p, "=parms[", 1:p,"]"), collapse = ',\n'),
                                ",\n",
                                "t = t,\n",
                                paste(paste0("x",1:n, "=y[", 1:n,"]"), collapse = ',\n'),
                                "))}"), collapse = " \n")),
       envir = envir)
}


generate.get.dx_dj <- function(p,
                               n,
                               jsens,
                               ct.lst,
                               rct.lst,
                               envir = .GlobalEnv){
  generate.get.dx.forgrad(p = p,
                          n = n,
                          ct.lst = ct.lst,
                          rct.lst = rct.lst,
                          envir = environment())
  
  get.dx_djs_cbind <- Deriv(f = get.dx.forgrad,
                            x = c(paste0("p", jsens)),
                            combine = "cbind",
                            cache.exp = FALSE)
  
  get.dx_djs_cbind_deparsed <- deparse(get.dx_djs_cbind)
  get.dx_djs_cbind_deparsed[1] <- paste0("get.dx_djs_cbind <- ", get.dx_djs_cbind_deparsed[1])
  eval(parse(text=paste0(get.dx_djs_cbind_deparsed, collapse = " \n")), envir = envir)
  
  eval(parse(text=paste0(paste0("get.dx_dj <- function(t,
                   y,
                   parms){\nreturn(c(get.dx_djs_cbind(",
                                paste(paste0("p",1:p, "=parms[", 1:p,"]"), collapse = ',\n'),
                                ",\n",
                                "t = t,\n",
                                paste(paste0("x",1:n, "=y[", 1:n,"]"), collapse = ',\n'),
                                ")))}"), collapse = " \n")),
       envir = envir)
}

generate.get.dxdP_dj <- function(p,
                               n,
                               jsens,
                               envir = .GlobalEnv){
  generate.get.dxdP.forgrad(p = p,
                          n = n,
                          envir = environment())

  get.dxdP_djs_c <- Deriv(f = get.dxdP.forgrad,
                          x = c(paste0("p", jsens)),
                          combine = "c",
                          cache.exp = FALSE)

  get.dxdP_djs_c_deparsed <- deparse(get.dxdP_djs_c)
  get.dxdP_djs_c_deparsed[1] <- paste0("get.dxdP_djs_c <- ", get.dxdP_djs_c_deparsed[1])
  eval(parse(text=paste0(get.dxdP_djs_c_deparsed, collapse = " \n")), envir = envir)

  eval(parse(text=paste0(paste0("get.dxdP_dj <- function(t, y, parms){\nget.dxdP_djs_c(",
                                paste(paste0("p",1:p, "=parms[", 1:p,"]"), collapse = ',\n'),
                                ",\n",
                                "t = t,\n",
                                paste(paste0("x",1:n, "=y[", 1:n,"]"), collapse = ',\n'),
                                ")\n}"), collapse = " \n")),
       envir = envir)
}


generate.get.nl <- function(p,
                        jsens,
                        envir = .GlobalEnv){
  
  nl.string <- paste0("nl <- function(psi,\n",
                      ifelse(length(setdiff(1:p, jsens)) > 0, paste0(paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"), ",\n"), ""),
                      "y,
                   times,
                   V,
                   E,
                   hmax){
                      phi <- exp(psi)
                      theta <- head(phi, -2)
                      g <- tail(phi, 2)[1]
                      s2 <- tail(phi, 1)

                     x <- try(vode(y = get.x.steady(t = max(times),\n",
                      "parms = ",
                      ifelse(length(setdiff(1:p, jsens)) > 0, paste0("c(theta,\n",
                                                                     paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"),
                                                                     ")"), "theta"),
                      "),\ntimes = times,
                func = get.dx,
                jacfunc = get.Jf,
                mf = -21,
                hmax = hmax,
                parms = ",
                      ifelse(length(setdiff(1:p, jsens)) > 0, paste0("c(theta,\n",
                                                                     paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"),
                                                                     ")"), "theta"),
                      "),\n silent = TRUE)[,-1]
                mu <- x[,'O']*g*(V - E)
                n <- length(y)
                return(as.numeric(n*log(s2) +  t(y - mu) %*% ((y - mu)/s2)))\n}")
  
  eval(parse(text=nl.string), envir = envir)
}


generate.get.gnl <- function(p,
                         jsens,
                         envir = .GlobalEnv){
  
  gnl.string <- paste0("gnl <- function(psi,\n",
                       ifelse(length(setdiff(1:p, jsens)) > 0, paste0(paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"), ",\n"), ""),
                       "y,
                   times,
                   V,
                   E,
                   hmax){
                      phi <- exp(psi)
                      theta <- head(phi, -2)
                      g <- tail(phi, 2)[1]
                      s2 <- tail(phi, 1)

                     x <- try(vode(y = c(get.x.steady(t = max(times),\n",
                       "parms = ",
                       ifelse(length(setdiff(1:p, jsens)) > 0, paste0("c(theta,\n",
                                                                      paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"),
                                                                      ")"), "theta"),
                       "),
                      get.x.steady_djs(t = max(times),\n",
                       "parms = ",
                       ifelse(length(setdiff(1:p, jsens)) > 0, paste0("c(theta,\n",
                                                                      paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"),
                                                                      ")"), "theta"),
                       ")),
                times = times,
                func = get.dx.wg,
                jacfunc = get.Jf.wg,
                mf = -21,
                hmax = hmax,
                parms = ",
                       ifelse(length(setdiff(1:p, jsens)) > 0, paste0("c(theta,\n",
                                                                      paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"),
                                                                      ")"), "theta"),
                       "),\n silent = TRUE)[,-1]

                mu <- x[,'O']*g*(V - E)
                dmu_th <- x[,which(1:ncol(x) %% (ncol(x)/(length(theta) + 1)) == which(colnames(x) == 'O'))[-1]]*g*(V - E)
                dmu_g <- x[,'O']*(V - E)

                n <- length(y)
                nl_th <- (-2*t(dmu_th) %*% ((y - mu)/s2))*theta
                nl_g <- (-2*t(dmu_g) %*% ((y - mu)/s2))*g
                nl_s2 <- (n - t(y - mu) %*% ((y - mu)/s2))
                
                return(c(nl_th, nl_g, nl_s2))\n}")
  
  eval(parse(text=gnl.string), envir = envir)
  
}

generate.get.xP.steady_djs <- function(p,
                                      jsens,
                                      envir = .GlobalEnv){
  
  generate.get.xP.steady.forgrad(p = p,
                                envir = environment())
  
  get.xP.steady_djs_c <- Deriv(f = get.xP.steady.forgrad,
                              x = c(paste0("p", jsens)),
                              combine = "c",
                              cache.exp = FALSE)
  
  get.xP.steady_djs_c_deparsed <- deparse(get.xP.steady_djs_c)
  get.xP.steady_djs_c_deparsed[1] <- paste0("get.xP.steady_djs_c <- ", get.xP.steady_djs_c_deparsed[1])
  eval(parse(text=paste0(get.xP.steady_djs_c_deparsed, collapse = " \n")), envir = envir)
  
  eval(parse(text=paste0(paste0("get.xP.steady_djs <- function(t,
                   parms){\nget.xP.steady_djs_c(",
                                paste(paste0("p",1:p, "=parms[", 1:p,"]"), collapse = ',\n'),
                                ",\n",
                                "t = t)\n}"), collapse = " \n")),
       envir = envir)
}

#################################

# nl <- function(psi,
#                p17,
#                y,
#                times,
#                V,
#                E,
#                hmax){
#   phi <- exp(psi)
#   theta <- head(phi, -2)
#   g <- tail(phi, 2)[1]
#   s2 <- tail(phi, 1)
#   
#   x <- try(vode(y = get.x.steady(t = max(times),
#                                  parms = c(theta,
#                                            p17)),
#                 times = times,
#                 func = get.dx,
#                 jacfunc = get.Jf,
#                 mf = -21,
#                 hmax = hmax,
#                 parms = c(theta,
#                           p17)),
#            silent = TRUE)[,-1]
#   mu <- x[,'O']*g*(V - E)
#   n <- length(y)
#   return(as.numeric(n*log(s2) +  t(y - mu) %*% ((y - mu)/s2)))
# }
# 
# 
# gnl <- function(psi,
#                 p17,
#                 y,
#                 times,
#                 V,
#                 E,
#                 hmax){
#   phi <- exp(psi)
#   theta <- head(phi, -2)
#   g <- tail(phi, 2)[1]
#   s2 <- tail(phi, 1)
#   
#   x <- try(vode(y = c(get.x.steady(t = max(times),
#                                    parms = c(theta,
#                                              p17)),
#                       get.x.steady_djs(t = max(times),
#                                        parms = c(theta,
#                                                  p17))),
#                 times = times,
#                 func = get.dx.wg,
#                 jacfunc = get.Jf.wg,
#                 mf = -21,
#                 hmax = hmax,
#                 parms = c(theta,
#                           p17)),
#            silent = TRUE)[,-1]
#   
#   mu <- x[,'O']*g*(V - E)
#   dmu_th <- x[,which(1:ncol(x) %% (ncol(x)/(length(theta) + 1)) == which(colnames(x) == 'O'))[-1]]*g*(V - E)
#   dmu_g <- x[,'O']*(V - E)
#   
#   n <- length(y)
#   nl_th <- (-2*t(dmu_th) %*% ((y - mu)/s2))*theta
#   nl_g <- (-2*t(dmu_g) %*% ((y - mu)/s2))*g
#   nl_s2 <- (n - t(y - mu) %*% ((y - mu)/s2))
#   
#   return(c(nl_th, nl_g, nl_s2))
# }
# 
# 
# tril(diag(head(chs,-1) == "O"))
# 
# which(diag(head(chs,-1) == "O")[lower.tri(diag(head(chs,-1) == "O"), diag = TRUE)]) + length(head(chs,-1)) # 30

get.nl <- function(psi,
                   gamma,
                   times,
                   V,
                   E,
                   y,
                   hmax){
  
  phi <- exp(psi)
  
  theta <- head(phi, -3)
  gs <- phi[length(theta) + 1]
  s2 <- phi[length(theta) + 2]
  nu <- 1 + phi[length(theta) + 3]
  
  xP <- try(vode(y = get.xP.steady(t = max(times),
                                    parms = c(theta, gamma)), # nStates*(nStates - 1)/2 + nStates = 10
                 times = times,
                 func = get.dxdP,
                 jacfunc = get.Jf.xP,
                 mf = -21,
                 hmax = hmax, # max(tps),
                 parms = c(theta, gamma),
                 rtol = 1e-6),
            silent = TRUE)[,-1]
  
  mu <- xP[,"O"]*gs*nu*(V - E)
  s <- gs^2*nu*(V - E)^2*xP[,30]*diff(times)[1] + s2
  
  nl <- as.numeric(sum(log(s)) +  t(y - mu) %*% ((y - mu)/s))
  return(nl)
}


get.gnl <- function(psi,
                    gamma,
                    times,
                    V,
                    E,
                    y,
                    hmax){

  phi <- exp(psi)

  theta <- head(phi, -3)
  gs <- phi[length(theta) + 1]
  s2 <- phi[length(theta) + 2]
  nu <- 1 + phi[length(theta) + 3]

  xP.wg <- try(vode(y = c(get.xP.steady(t = max(times),
                                       parms = c(theta, gamma)),
                          get.xP.steady_djs(t = max(times),
                                           parms = c(theta, gamma))), #
                    times = times,
                    func = get.dxdP.wg,
                    jacfunc = get.Jf.xP.wg,
                    mf = -21,
                    hmax = hmax, # max(tps),
                    parms = c(theta, gamma),
                    rtol = c(rep(1e-6, 7),
                             rep(1e-6, 28),
                             rep(1e-6, 560))),
               silent = TRUE)[,-1]

  mu <- xP.wg[,"O"]*gs*nu*(V - E)
  s <- gs^2*nu*(V - E)^2*xP.wg[,30]*diff(times)[1] + s2

  dmu_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/17) == 5)[-1]]*gs*nu*(V - E) # which(chs == "O") = 5; length(theta) = 17
  dmu_gs <- xP.wg[,"O"]*nu*(V - E)
  dmu_nu <- xP.wg[,"O"]*gs*(V - E)

  ds_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/17) == 30)[-1]]*gs^2*nu*(V - E)^2*diff(times)[1]
  ds_gs <- 2*gs*nu*(V - E)^2*xP.wg[,30]*diff(times)[1]
  ds_s2 <- 1
  ds_dnu <- gs^2*(V - E)^2*xP.wg[,30]*diff(times)[1]

  gnl_theta <- as.numeric(colSums(ds_theta/s)
                          -as.numeric(2*t(dmu_theta) %*% ((y - mu)/s))
                          -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_theta*(y - mu))))*theta

  gnl_gs <- as.numeric(sum(ds_gs/s)
                       -as.numeric(2*t(dmu_gs) %*% ((y - mu)/s))
                       -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_gs*(y - mu))))*gs

  gnl_s2 <- as.numeric(sum(1/s) - t(y - mu) %*% ((y - mu)/s^2))*s2
  gnl_nu <- as.numeric(sum(ds_dnu/s)
                       -as.numeric(2*t(dmu_nu) %*% ((y - mu)/s))
                       -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_dnu*(y - mu))))*exp(tail(psi,1))

  return(c(gnl_theta,
           gnl_gs,
           gnl_s2,
           gnl_nu))

}


###############

get.nl.nus2 <- function(psi,
                        g,
                        gamma,
                        times,
                        V,
                        E,
                        y,
                        xP.wg){
  
  phi <- exp(psi)
  s2 <- phi[1]
  nu <- 1 + phi[2]
  
  mu <- xP.wg[,"O"]*g*(V - E)
  s <- g^2*(V - E)^2*xP.wg[,30]/nu*diff(times)[1] + s2
  
  nl <- as.numeric(sum(log(s)) +  t(y - mu) %*% ((y - mu)/s))
  return(nl)
}

get.gnl.nus2 <- function(psi,
                         g,
                         gamma,
                         times,
                         V,
                         E,
                         y,
                         xP.wg){
  
  phi <- exp(psi)
  s2 <- phi[1]
  nu <- 1 + phi[2]
  
  mu <- xP.wg[,"O"]*g*(V - E)
  s <- g^2*(V - E)^2*xP.wg[,30]/nu*diff(times)[1] + s2
  
  dmu_s2 <- 0
  ds_dnu <- -(g^2*(V - E)^2*xP.wg[,30]/nu^2)*diff(times)[1]
  
  gnl_s2 <- as.numeric(sum(1/s) - t(y - mu) %*% ((y - mu)/s^2))*s2
  gnl_nu <- as.numeric(sum(ds_dnu/s)
                       -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_dnu*(y - mu))))*exp(tail(psi,1))
  
  return(c(gnl_s2,
           gnl_nu))
}

##########################

get.x.em <- function(theta,
                     nu,
                     times,
                     chs.lst){
  
  Y <- Y.mean <- matrix(data = NA,
                        nrow = length(times),
                        ncol = length(chs.lst) - 1)
  colnames(Y) <- colnames(Y.mean) <- head(chs.lst,-1)
  Y.mean[1,] <- Y[1,] <- get.logit(get.x.steady(t = max(times),
                                                parms = theta))
  
  wPidx <- which(lower.tri(matrix(data = NA,
                                  nrow = length(chs.lst) - 1,
                                  ncol = length(chs.lst) - 1), diag = TRUE), arr.ind = TRUE)
  P <- matrix(NA, length(chs.lst) - 1, length(chs.lst) - 1)
  
  N <- length(chs.lst) - 1
  gi <- gj <- matrix(0,N,N)
  
  
  disc.times <- sapply(c(250,
                         300,
                         500,
                         1500,
                         2000,
                         3000,
                         6500,
                         7000),
                       function(t){which(times == t)})
  
  nStates <- length(chs.lst) - 1
  xP <- vode(y = c(get.xP.steady(t = max(times),
                                parms = theta)),
             times = times,
             func = get.dxdP,
             jacfunc = get.Jf.xP,
             mf = -21,
             hmax = .1, # max(tps),
             parms = theta,
             rtol = 1e-6)[,-1]
  
  # xP.wg <- vode(y = c(get.x.steady(t = max(times),
  #                                  parms = theta),
  #                     rep(0, 28), # nStates*(nStates - 1)/2 + nStates
  #                     get.x.steady_djs(t = max(times),
  #                                      parms = theta),
  #                     rep(0, 448)),
  #               times = times,
  #               func = get.dxdP.wg,
  #               jacfunc = get.Jf.xP.wg,
  #               mf = -21,
  #               hmax = max(times), # max(tps)
  #               parms = theta,
  #               rtol = 1e-6)[,-1]
  # xP <- xP.wg[,1:35]
  
  for (k in 2:length(times)) {
    cat(paste0("\t t = ", times[k], "\n"))
    Dt <- times[k] - times[k-1]
    
    Jg.x <- get.Jg.x(x = get.logistic(Y.mean[k-1,]))
    Nablag.x <- get.d2g.x(get.logistic(Y.mean[k-1,]))
    Y.mean[k,] <- get.logit(xP[k, 1:N])
    Sigma.t <- get.Sigma(x = xP[k, 1:N],
                         lambda = get.lambda(theta = theta,
                                             V = get.voltage(t = times[k-1])))/nu
    
    P[rbind(wPidx,wPidx[,c(2,1)])] <- xP[k,-c(1:N)]
    P.t <- P/nu
    Sigma.t <- P.t
    Sigma.tr.t <- Jg.x %*% (Sigma.t) %*% t(Jg.x)
    
    if((k - 1) %in% disc.times){
      Y[k-1,] <- Y.mean[k-1,]
    }
    
    tY <- vode(y = get.logistic(Y[k-1,]),
               times = times[c(k-1,k)],
               func = get.dx,
               jacfunc = get.Jf,
               mf = -21,
               hmax = .1, # max(tps),
               parms = theta,
               rtol = 1e-6)
    
    if(sum(tY[2,-1] < 0)>0){
      tY[2,-1] <- get.logistic(Y.mean[k,])
    }
    
    Y[k,] <- (get.logit(tY[2,-1]) + .5*sapply(1:N, function(j){
      gj[j,j] <-  Nablag.x[j]
      tr(Sigma.t %*% gj)
    })*Dt + rmvnorm(n = 1,
                    mean = rep(0, length(Y[k-1,])),
                    sigma = Dt*Sigma.tr.t)[1,])
    
    b <- 1
    maxit <- 1 # 10000
    nit <- 0
    q.tol <- 0.027
    
    while ((sum(Y[k,] > Y.mean[k,] + .5*sapply(1:N, function(j){
      gj[j,j] <-  Nablag.x[j]
      tr(P.t %*% gj)
    })*Dt + qnorm(1-q.tol)*sqrt(diag(Dt*(Jg.x %*% (P.t) %*% t(Jg.x))))) > 0 |
    sum(Y[k,] < Y.mean[k,] + .5*sapply(1:N, function(j){
      gj[j,j] <-  Nablag.x[j]
      tr(P.t %*% gj)
    })*Dt - qnorm(1-q.tol)*sqrt(diag(Dt*(Jg.x %*% (P.t) %*% t(Jg.x))))) > 0 |
    # sum(get.logistic(Y[k,])) > b)
    ifelse(is.na(sum(get.logistic(Y[k,])) > b), TRUE, sum(get.logistic(Y[k,])) > b))
    & nit < maxit) {
      Y[k,] <- (get.logit(tY[2,-1]) + .5*sapply(1:N, function(j){
        gj[j,j] <-  Nablag.x[j]
        tr(Sigma.t %*% gj)
      })*Dt + rmvnorm(n = 1,
                      mean = rep(0, length(Y[k-1,])),
                      sigma = Dt*Sigma.tr.t)[1,])
      nit <- nit + 1
    }
    
    if((sum(Y[k,] > Y.mean[k,] + .5*sapply(1:N, function(j){
      gj[j,j] <-  Nablag.x[j]
      tr(P.t %*% gj)
    })*Dt + qnorm(1-q.tol)*sqrt(diag(Dt*(Jg.x %*% (P.t) %*% t(Jg.x))))) > 0 |
    sum(Y[k,] < Y.mean[k,] + .5*sapply(1:N, function(j){
      gj[j,j] <-  Nablag.x[j]
      tr(P.t %*% gj)
    })*Dt - qnorm(1-q.tol)*sqrt(diag(Dt*(Jg.x %*% (P.t) %*% t(Jg.x))))) > 0 |
    # sum(get.logistic(Y[k,])) > b)){
    ifelse(is.na(sum(get.logistic(Y[k,])) > b), TRUE, sum(get.logistic(Y[k,])) > b))){
      Y[k,] <- Y.mean[k,]
    }
    
    if(tY[2,"time"] != times[k]){
      Y[k,] <- Y.mean[k,]
    }
    
    if(sum(is.nan(Y[k,])) > 0 |
       sum(is.infinite(Y[k,])) > 0 |
       sum(is.na(Y[k,])) > 0){
      Y[k,] <- Y.mean[k,]
    }
    
    if(sum(is.nan(get.logistic(Y[k,]))) > 0 |
       sum(is.infinite(get.logistic(Y[k,]))) > 0 |
       sum(is.na(get.logistic(Y[k,]))) > 0){
      Y[k,] <- Y.mean[k,]
    }
  }
  return(Y)
}

get.pred <- function(psi,
                     gamma,
                     times,
                     V,
                     E,
                     hmax){
  
  phi <- exp(psi)
  
  theta <- head(phi, -3)
  gs <- phi[length(theta) + 1]
  s2 <- phi[length(theta) + 2]
  nu <- 1 + phi[length(theta) + 3]
  
  xP <- try(vode(y = get.xP.steady(t = max(times),
                                   parms = c(theta, gamma)), # nStates*(nStates - 1)/2 + nStates = 10
                 times = times,
                 func = get.dxdP,
                 jacfunc = get.Jf.xP,
                 mf = -21,
                 hmax = hmax, # max(tps),
                 parms = c(theta, gamma),
                 rtol = 1e-6),
            silent = TRUE)[,-1]
  
  mu <- xP[,"O"]*gs*nu*(V - E)
  s <- gs^2*nu*(V - E)^2*xP[,30]*diff(times)[1] + s2
  
  return(cbind(mu, s))
}

CMAES.control = cmaes.control()
########
CMAES.control$options$StopOnWarnings = FALSE
CMAES.control$options$TolFun = 1e-8
CMAES.control$options$TolX = 1e-8
CMAES.control$options$DispModulo=10
CMAES.control$options$MaxIter = 1000
########