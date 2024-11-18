
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
                      cache.exp = TRUE)

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

get.dxdP <- function(t = 0,
                     y,
                     parms){

  lambda <- get.lambda(theta = parms,
                       V = get.voltage(t = t))

  return(list(dxdP = (get.Jfy.y(lambda) %*% y + get.bfy(lambda))[,1]))
}

get.Jf.xP <- function(t = 0,
                      y,
                      parms){

  return(get.Jfy.y(get.lambda(theta = parms,
                              V = get.voltage(t = t))))
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
                          cache.exp = TRUE)

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


#' @keywords internal
generate.get.dxdP.forgrad <- function(p,
                                      n,
                                      envir = .GlobalEnv){

  get.dxdP.forgrad.string <- paste0("get.dxdP.forgrad <- function(",
                                    paste0(paste(paste0("p", 1:p), collapse = ",\n"), ",\n"),
                                    "t,\n",
                                    paste(paste0("x", 1:n), collapse = ",\n"),
                                    "){\n",
                                    paste0("parms <- c(",
                                           paste(paste0("p", 1:p), collapse = ",\n"),
                                           ")\n"),
                                    paste0("y <- c(",
                                           paste(paste0("x", 1:n), collapse = ",\n"),
                                           ")\n"),
                                    "lambda <- get.lambda(theta = parms,\n V = get.voltage(t = t))\n",
                                    "c(get.Jfy.y(lambda) %*% y + get.bfy(lambda))\n}")

  eval(parse(text=get.dxdP.forgrad.string), envir = envir)
}


generate.get.dxdP.wg <- function(n,
                                 p,
                                 envir = .GlobalEnv){

  get.dxdP.wg.string <- paste0("get.dxdP.wg <- function(t = 0,
                        y,
                        parms){
  lambda <- get.lambda(theta = parms,
                       V = get.voltage(t = t))
  Jfy.y <- get.Jfy.y(lambda)
       return(list(dx = c(Jfy.y %*% y[1:", n,"] + get.bfy(lambda),
                     Jfy.y %*% matrix(y[", n+1, ":", n*(p+1),"], ncol = ", p,") + get.dxdP_dj(t = t,
                                                                         y = y[1:", n,"],
                                                                         parms = parms))))\n}")

  eval(parse(text=get.dxdP.wg.string), envir = envir)
}


get.Jf.xP.wg <- function(t = 0,
                         y,
                         parms){

  Jfull.xp[bdiag.idxs.xp] <- get.Jfy.y(get.lambda(theta = parms,
                                                  V = get.voltage(t = t)))

  Jfull.xp[Jh.idxs.xp] <- get.Jfy.y.dj(t = t,
                                       parms = parms)
  return(Jfull.xp)
}


generate.bdiagmat.xp <- function(n,
                                 p,
                                 envir = .GlobalEnv){

  bdiag.idxs.xp.string <- paste0("bdiag.idxs.xp <- cbind(c(matrix(1:",
                                 n*(p+1),", nrow = ", n, ")[,rep(1:", p+1, ", each = ", n, ")]),\n",
                                 "rep(1:", n*(p+1), ", each = ", n, "))")

  Jh.idxs.xp.string <- paste0("Jh.idxs.xp <- cbind(rep(",
                              n+1, ":", n*(p+1), ", times = ", n, "),\n",
                              "rep(", 1, ":", n, ", each = ", n*((p+1)-1), "))")

  Jfull.xp.string <- paste0("Jfull.xp <- matrix(", 0, ", ", n*(p+1), ", ", n*(p+1), ")")


  eval(parse(text=bdiag.idxs.xp.string), envir = envir)
  eval(parse(text=Jh.idxs.xp.string), envir = envir)
  eval(parse(text=Jfull.xp.string), envir = envir)

}

nl.try <- function(psi,
                   gamma,
                   times,
                   V,
                   E,
                   y,
                   hmax){
  nl.res <- try(quiet(get.nl(psi = psi,
                             gamma = gamma,
                             times = times,
                             V = V,
                             E = E,
                             y = y,
                             hmax = hmax)), silent = TRUE)
  if(inherits(nl.res,'try-error') |
     is.nan(nl.res) |
     is.infinite(nl.res)){
    return(1e8)
  }
  return(nl.res)
}


gnl.try <- function(psi,
                    gamma,
                    times,
                    V,
                    E,
                    y,
                    hmax){
  gnl.res <- try(quiet(get.gnl(psi = psi,
                               gamma = gamma,
                               times = times,
                               V = V,
                               E = E,
                               y = y,
                               hmax = hmax)), silent = TRUE)
  if(inherits(gnl.res,'try-error') |
     sum(is.nan(gnl.res)) > 0 |
     sum(is.infinite(gnl.res)) > 0){
    return(rep(1e8, length(psi)))
  }
  return(gnl.res)
}

generate.get.xP.steady <- function(ct.lst,
                                  envir = .GlobalEnv){
  
  get.xP.steady.string <- paste0("get.xP.steady <- function(t = 0,
                         parms){\n",
                                "lambda <- get.lambda(theta = parms,\n V = get.voltage(t = t))\n",
                                "xP.s <- (- solve(get.Jfy.y(lambda)) %*% get.bfy(lambda))[,1]\n",
                                "names(xP.s)[1:", length(ct.lst) - 1,"] <- c(",
                                paste(sapply(head(ct.lst,-1), function(s){paste0("'", s, "'")}), collapse = ",\n"), ")\n",
                                "return(xP.s)\n}")
  
  eval(parse(text=get.xP.steady.string), envir = envir)
  
}

generate.get.xP.steady.forgrad <- function(p,
                                          envir = .GlobalEnv){
  
  get.xP.steady.forgrad.string <- paste0("get.xP.steady.forgrad <- function(",
                                        paste(paste0("p", 1:p), collapse = ",\n"),
                                        ",\nt){\n",
                                        "parms <- c(",
                                        paste(paste0("p", 1:p), collapse = ",\n"),
                                        ")\n",
                                        "lambda <- get.lambda(theta = parms,\n V = get.voltage(t = t))\n",
                                        "Vl <- get.Vlambda(lambda)\n",
                                        "(- solve(get.Jfy.y(lambda)) %*% get.bfy(lambda))[,1]\n}")
  
  eval(parse(text=get.xP.steady.forgrad.string), envir = envir)
}

generate.get.xP.steady_djs <- function(p,
                                      jsens,
                                      envir = .GlobalEnv){
  
  generate.get.xP.steady.forgrad(p = p,
                                envir = environment())
  
  get.xP.steady_djs_c <- Deriv(f = get.xP.steady.forgrad,
                              x = c(paste0("p", jsens)),
                              combine = "c",
                              cache.exp = TRUE)
  
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
  s <- gs^2*nu*(V - E)^2*xP[,9]*diff(times)[1] + s2

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
                                           parms = c(theta, gamma))),
                    times = times,
                    func = get.dxdP.wg,
                    jacfunc = get.Jf.xP.wg,
                    mf = -21,
                    hmax = hmax, # max(tps),
                    parms = c(theta, gamma),
                    rtol = c(rep(1e-6, 4),
                             rep(1e-6, 10),
                             rep(1e-6, 112))),
               silent = TRUE)[,-1]

  mu <- xP.wg[,"O"]*gs*nu*(V - E)
  s <- gs^2*nu*(V - E)^2*xP.wg[,9]*diff(times)[1] + s2

  dmu_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 2)[-1]]*gs*nu*(V - E)
  dmu_gs <- xP.wg[,"O"]*nu*(V - E)
  dmu_nu <- xP.wg[,"O"]*gs*(V - E)

  ds_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 9)[-1]]*gs^2*nu*(V - E)^2*diff(times)[1]
  ds_gs <- 2*gs*nu*(V - E)^2*xP.wg[,9]*diff(times)[1]
  ds_s2 <- 1
  ds_dnu <- gs^2*(V - E)^2*xP.wg[,9]*diff(times)[1]

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
# get.nl <- function(psi,
#                    gamma,
#                    times,
#                    V,
#                    E,
#                    y,
#                    hmax){
# 
#   phi <- exp(psi)
# 
#   theta <- head(phi, -3)
#   gs <- phi[length(theta) + 1]
#   s2 <- phi[length(theta) + 2]
#   nu <- 1 + phi[length(theta) + 3]
# 
#   xP <- try(vode(y = c(get.x.steady(t = max(times),
#                                     parms = c(theta, gamma)),
#                        rep(0, 10)), # nStates*(nStates - 1)/2 + nStates = 10
#                  times = times,
#                  func = get.dxdP,
#                  jacfunc = get.Jf.xP,
#                  mf = -21,
#                  hmax = hmax, # max(tps),
#                  parms = c(theta, gamma),
#                  rtol = 1e-6),
#             silent = TRUE)[,-1]
# 
#   mu <- xP[,"O"]*gs*nu*(V - E)
#   s <- gs^2*nu*(V - E)^2*xP[,9]*diff(times)[1] + s2
# 
#   nl <- as.numeric(sum(log(s)) +  t(y - mu) %*% ((y - mu)/s))
#   return(nl)
# }
# 
# 
# get.gnl <- function(psi,
#                     gamma,
#                     times,
#                     V,
#                     E,
#                     y,
#                     hmax){
# 
#   phi <- exp(psi)
# 
#   theta <- head(phi, -3)
#   gs <- phi[length(theta) + 1]
#   s2 <- phi[length(theta) + 2]
#   nu <- 1 + phi[length(theta) + 3]
# 
#   xP.wg <- try(vode(y = c(get.x.steady(t = max(times),
#                                        parms = c(theta, gamma)),
#                           rep(0, 10), # nStates*(nStates - 1)/2 + nStates
#                           get.x.steady_djs(t = max(times),
#                                            parms = c(theta, gamma)),
#                           rep(0, 80)),
#                     times = times,
#                     func = get.dxdP.wg,
#                     jacfunc = get.Jf.xP.wg,
#                     mf = -21,
#                     hmax = hmax, # max(tps),
#                     parms = c(theta, gamma),
#                     rtol = c(rep(1e-6, 4),
#                              rep(1e-6, 10),
#                              rep(1e-6, 112))),
#                silent = TRUE)[,-1]
# 
#   mu <- xP.wg[,"O"]*gs*nu*(V - E)
#   s <- gs^2*nu*(V - E)^2*xP.wg[,9]*diff(times)[1] + s2
# 
#   dmu_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 2)[-1]]*gs*nu*(V - E)
#   dmu_gs <- xP.wg[,"O"]*nu*(V - E)
#   dmu_nu <- xP.wg[,"O"]*gs*(V - E)
# 
#   ds_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 9)[-1]]*gs^2*nu*(V - E)^2*diff(times)[1]
#   ds_gs <- 2*gs*nu*(V - E)^2*xP.wg[,9]*diff(times)[1]
#   ds_s2 <- 1
#   ds_dnu <- gs^2*(V - E)^2*xP.wg[,9]*diff(times)[1]
# 
#   gnl_theta <- as.numeric(colSums(ds_theta/s)
#                           -as.numeric(2*t(dmu_theta) %*% ((y - mu)/s))
#                           -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_theta*(y - mu))))*theta
# 
#   gnl_gs <- as.numeric(sum(ds_gs/s)
#                        -as.numeric(2*t(dmu_gs) %*% ((y - mu)/s))
#                        -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_gs*(y - mu))))*gs
# 
#   gnl_s2 <- as.numeric(sum(1/s) - t(y - mu) %*% ((y - mu)/s^2))*s2
#   gnl_nu <- as.numeric(sum(ds_dnu/s)
#                        -as.numeric(2*t(dmu_nu) %*% ((y - mu)/s))
#                        -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_dnu*(y - mu))))*exp(tail(psi,1))
# 
#   return(c(gnl_theta,
#            gnl_gs,
#            gnl_s2,
#            gnl_nu))
# 
# }
# get.gnl <- function(psi,
#                     gamma,
#                     times,
#                     V,
#                     E,
#                     y,
#                     hmax){
#   
#   phi <- exp(psi)
#   
#   theta <- head(phi, -3)
#   gs <- phi[length(theta) + 1]
#   s2 <- phi[length(theta) + 2]
#   nu <- 1 + phi[length(theta) + 3]
#   
#   dx_th.idx <- tail(which(1:126 %% (126/9) %in% 1:4), 32)
#   dP_th.idx <- tail(which(1:126 %% (126/9) %in% c(0,5:14)), 80)  
#   xP_dxP_0 <- rep(NA, 126)
#   xP_dxP_0[1:4] <- get.x.steady(t = max(times),
#                                 parms = c(theta, gamma))
#   xP_dxP_0[5:14] <- rep(0, 10)
#   xP_dxP_0[dx_th.idx] <- get.x.steady_djs(t = max(times),
#                                           parms = c(theta, gamma))
#   xP_dxP_0[dP_th.idx] <- rep(0, 80)
#   
#   xP.wg <- try(vode(y = xP_dxP_0,
#                     times = times,
#                     func = get.dxdP.wg,
#                     jacfunc = get.Jf.xP.wg,
#                     mf = -21,
#                     hmax = hmax, # max(tps),
#                     parms = c(theta, gamma),
#                     rtol = c(rep(1e-6, 4),
#                              rep(1e-6, 10),
#                              rep(1e-6, 112))),
#                silent = TRUE)[,-1]
#   
#   mu <- xP.wg[,"O"]*gs*nu*(V - E)
#   s <- gs^2*nu*(V - E)^2*xP.wg[,9]*diff(times)[1] + s2
#   
#   dmu_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 2)[-1]]*gs*nu*(V - E)
#   dmu_gs <- xP.wg[,"O"]*nu*(V - E)
#   dmu_nu <- xP.wg[,"O"]*gs*(V - E)
#   
#   ds_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 9)[-1]]*gs^2*nu*(V - E)^2*diff(times)[1]
#   ds_gs <- 2*gs*nu*(V - E)^2*xP.wg[,9]*diff(times)[1]
#   ds_s2 <- 1
#   ds_dnu <- gs^2*(V - E)^2*xP.wg[,9]*diff(times)[1]
#   
#   gnl_theta <- as.numeric(colSums(ds_theta/s)
#                           -as.numeric(2*t(dmu_theta) %*% ((y - mu)/s))
#                           -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_theta*(y - mu))))*theta
#   
#   gnl_gs <- as.numeric(sum(ds_gs/s)
#                        -as.numeric(2*t(dmu_gs) %*% ((y - mu)/s))
#                        -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_gs*(y - mu))))*gs
#   
#   gnl_s2 <- as.numeric(sum(1/s) - t(y - mu) %*% ((y - mu)/s^2))*s2
#   gnl_nu <- as.numeric(sum(ds_dnu/s)
#                        -as.numeric(2*t(dmu_nu) %*% ((y - mu)/s))
#                        -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_dnu*(y - mu))))*exp(tail(psi,1))
#   
#   return(c(gnl_theta,
#            gnl_gs,
#            gnl_s2,
#            gnl_nu))
#   
# }

get.fim <- function(psi,
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

  xP.wg <- try(vode(y = c(get.x.steady(t = max(times),
                                       parms = c(theta, gamma)),
                          rep(0, 10), # nStates*(nStates - 1)/2 + nStates
                          get.x.steady_djs(t = max(times),
                                           parms = c(theta, gamma)),
                          rep(0, 80)),
                    times = times,
                    func = get.dxdP.wg,
                    jacfunc = get.Jf.xP.wg,
                    mf = -21,
                    hmax = hmax, # max(tps),
                    parms = c(theta, gamma),
                    rtol = c(rep(1e-6, 4),
                             rep(1e-6, 10),
                             rep(1e-6, 112))),
               silent = TRUE)[,-1]

  mu <- xP.wg[,"O"]*gs*nu*(V - E)
  s <- gs^2*nu*(V - E)^2*xP.wg[,9]*diff(times)[1] + s2

  dmu_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 2)[-1]]*gs*nu*(V - E)
  dmu_gs <- xP.wg[,"O"]*nu*(V - E)
  dmu_s2 <- rep(0, length(dmu_gs))
  dmu_nu <- xP.wg[,"O"]*gs*(V - E)

  dmu <- cbind(dmu_theta, dmu_gs, dmu_s2, dmu_nu)
  dmu <- dmu*matrix(data = phi,
                    ncol = length(psi),
                    nrow = nrow(dmu),
                    byrow = TRUE)

  ds_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 9)[-1]]*gs^2*nu*(V - E)^2*diff(times)[1]
  ds_gs <- 2*gs*nu*(V - E)^2*xP.wg[,9]*diff(times)[1]
  ds_s2 <- 1
  ds_dnu <- gs^2*(V - E)^2*xP.wg[,9]*diff(times)[1]

  ds <- cbind(ds_theta, ds_gs, ds_s2, ds_dnu)
  ds <- ds*matrix(data = phi,
                  ncol = length(psi),
                  nrow = nrow(ds),
                  byrow = TRUE)

  fim <- t(dmu) %*% (dmu/s) + .5*t(ds) %*% (ds/s^2)


  return(fim)

}

# get.pred <- function(psi,
#                      gamma,
#                      times,
#                      V,
#                      E,
#                      hmax){
# 
#   phi <- exp(psi)
# 
#   theta <- head(phi, -3)
#   gs <- phi[length(theta) + 1]
#   s2 <- phi[length(theta) + 2]
#   nu <- 1 + phi[length(theta) + 3]
# 
#   xP <- try(vode(y = c(get.x.steady(t = max(times),
#                                     parms = c(theta, gamma)),
#                        rep(0, 10)), # nStates*(nStates - 1)/2 + nStates = 10
#                  times = times,
#                  func = get.dxdP,
#                  jacfunc = get.Jf.xP,
#                  mf = -21,
#                  hmax = hmax, # max(tps),
#                  parms = c(theta, gamma),
#                  rtol = 1e-6),
#             silent = TRUE)[,-1]
# 
#   mu <- xP[,"O"]*gs*nu*(V - E)
#   s <- gs^2*nu*(V - E)^2*xP[,9]*diff(times)[1] + s2
# 
#   return(cbind(mu, s))
# }

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
  s <- gs^2*nu*(V - E)^2*xP[,9]*diff(times)[1] + s2
  
  return(cbind(mu, s))
}

get.Jg.x <- function(x){
  diag(1/(x*(1-x)))
}

get.d2g.x <- function(x){
  (2*x - 1)/(x^2*(1-x)^2)
}

tr <- function(A){
  sum(diag(A))
}

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
  
  gi <- gj <- matrix(0,4,4)
  
  
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
             hmax = 1, # max(tps),
             parms = theta,
             rtol = 1e-6)[,-1]
  # xP <- vode(y = c(get.x.steady(t = max(times),
  #                               parms = theta),
  #                  rep(0, nStates*(nStates - 1)/2 + nStates)),
  #            times = times,
  #            func = get.dxdP,
  #            jacfunc = get.Jf.xP,
  #            mf = -21,
  #            hmax = 1, # max(tps),
  #            parms = theta,
  #            rtol = 1e-6)[,-1]
  
  for (k in 2:length(times)) {
    cat(paste0("\t t = ", times[k], "\n"))
    Dt <- times[k] - times[k-1]
    
    Jg.x <- get.Jg.x(x = get.logistic(Y.mean[k-1,]))
    Nablag.x <- get.d2g.x(get.logistic(Y.mean[k-1,]))
    Y.mean[k,] <- get.logit(xP[k, 1:4])
    Sigma.t <- get.Sigma(x = xP[k, 1:4],
                         lambda = get.lambda(theta = theta,
                                             V = get.voltage(t = times[k-1])))/nu
    
    P[rbind(wPidx,wPidx[,c(2,1)])] <- xP[k,-c(1:4)]
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
    
    Y[k,] <- (get.logit(tY[2,-1]) + .5*sapply(1:4, function(j){
      gj[j,j] <-  Nablag.x[j]
      tr(Sigma.t %*% gj)
    })*Dt + rmvnorm(n = 1,
                    mean = rep(0, length(Y[k-1,])),
                    sigma = Dt*Sigma.tr.t)[1,])
    
    b <- 1
    maxit <- 10000
    nit <- 0
    q.tol <- 0.027
    
    while ((sum(Y[k,] > Y.mean[k,] + .5*sapply(1:4, function(j){
      gj[j,j] <-  Nablag.x[j]
      tr(P.t %*% gj)
    })*Dt + qnorm(1-q.tol)*sqrt(diag(Dt*(Jg.x %*% (P.t) %*% t(Jg.x))))) > 0 |
    sum(Y[k,] < Y.mean[k,] + .5*sapply(1:4, function(j){
      gj[j,j] <-  Nablag.x[j]
      tr(P.t %*% gj)
    })*Dt - qnorm(1-q.tol)*sqrt(diag(Dt*(Jg.x %*% (P.t) %*% t(Jg.x))))) > 0 |
    sum(get.logistic(Y[k,])) > b)
    & nit < maxit) {
      Y[k,] <- (get.logit(tY[2,-1]) + .5*sapply(1:4, function(j){
        gj[j,j] <-  Nablag.x[j]
        tr(Sigma.t %*% gj)
      })*Dt + rmvnorm(n = 1,
                      mean = rep(0, length(Y[k-1,])),
                      sigma = Dt*Sigma.tr.t)[1,])
      nit <- nit + 1
    }
    
    if((sum(Y[k,] > Y.mean[k,] + .5*sapply(1:4, function(j){
      gj[j,j] <-  Nablag.x[j]
      tr(P.t %*% gj)
    })*Dt + qnorm(1-q.tol)*sqrt(diag(Dt*(Jg.x %*% (P.t) %*% t(Jg.x))))) > 0 |
    sum(Y[k,] < Y.mean[k,] + .5*sapply(1:4, function(j){
      gj[j,j] <-  Nablag.x[j]
      tr(P.t %*% gj)
    })*Dt - qnorm(1-q.tol)*sqrt(diag(Dt*(Jg.x %*% (P.t) %*% t(Jg.x))))) > 0 |
    sum(get.logistic(Y[k,])) > b)){
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

get.y.em <- function(times,
                     X,
                     gs,
                     nu,
                     E,
                     s2){
  V <- sapply(times, function(t){get.voltage(t)})
  y <- c(X[,"O"]*gs*nu*(V - E) + rnorm(length(V), 0, sqrt(s2)))
  return(y)
}


generate.get.Sigma <- function(ct.lst,
                               rct.lst,
                               envir = .GlobalEnv){
  
  Q.string <- matrix(data = "0",
                     nrow = length(ct.lst) - 1,
                     ncol = length(ct.lst) - 1)
  rownames(Q.string) <- colnames(Q.string) <- head(ct.lst,-1)
  
  for (i in head(ct.lst,-1)) {
    for (j in head(ct.lst,-1)) {
      
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
        
        rgts <- sapply(rgts, function(rgt){
          if(rgt == tail(ct.lst,1)){
            rgt <- "(1 - sum(x))"
          }else{
            rgt <- paste0("x['", rgt, "']")
          }
          return(rgt)
        })
        
        if(i == j){
          Q.string[i,j] <- paste(paste0(rgts, "*lambda['", rct.lst[wri & wrj], "']"), collapse = " + ")
        }else{
          Q.string[i,j] <- paste0("-", paste(paste0(rgts, "*lambda['", rct.lst[wri & wrj], "']"), collapse = " - "))
        }
      }
    }
  }
  
  Q.string <- paste0("get.Sigma <- function(x, lambda){\nmatrix(data = c(",
                     paste(Q.string, collapse = ",\n"),
                     "), nrow = ", length(ct.lst) - 1, ", ncol = ", length(ct.lst) - 1, ")\n}")
  
  eval(parse(text=Q.string), envir = envir)
  
}

generate.get.mu <- function(ct.lst,
                            rct.lst,
                            envir = .GlobalEnv,
                            flickering = FALSE){
  
  if(flickering){
    rct.lst <- c(rct.lst,
                 "F->O",
                 "O->F")
  }
  
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
  
  b <- Vth[-ncol(Vth),ncol(Vth)]
  
  get.dx.string <- paste0("get.mu <- function(t = 0,
                            x,
                            parms){\n
  lambda <- get.lambda(theta = parms,
                       V = get.voltage(t = t))\n
  Vl <- get.Vlambda(lambda)\n",
                          "mu <- Vl %*% x + c(",paste(b, collapse = ",\n"),")\n",
                          "return(mu[,1])\n}")
  
  eval(parse(text=get.dx.string), envir = envir)
}

nl.nus2.try <- function(psi,
                        g,
                        gamma,
                        times,
                        V,
                        E,
                        y,
                        xP.wg){
  nl.res <- try(get.nl.nus2(psi = psi,
                            g = g,
                            gamma = gamma,
                            times = times,
                            V = V,
                            E = E,
                            y = y,
                            xP.wg = xP.wg), silent = TRUE)
  if(inherits(nl.res,'try-error') |
     is.nan(nl.res) |
     is.infinite(nl.res)){
    return(1e8)
  }
  return(nl.res)
}

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
  s <- g^2*(V - E)^2*xP.wg[,9]/nu*diff(times)[1] + s2
  
  nl <- as.numeric(sum(log(s)) +  t(y - mu) %*% ((y - mu)/s))
  return(nl)
}


gnl.nus2.try <- function(psi,
                         g,
                         gamma,
                         times,
                         V,
                         E,
                         y,
                         xP.wg){
  gnl.res <- try(get.gnl.nus2(psi = psi,
                              g = g,
                              gamma = gamma,
                              times = times,
                              V = V,
                              E = E,
                              y = y,
                              xP.wg = xP.wg), silent = TRUE)
  if(inherits(gnl.res,'try-error') |
     sum(is.nan(gnl.res)) > 0 |
     sum(is.infinite(gnl.res)) > 0){
    return(rep(1e8, length(psi)))
  }
  return(gnl.res)
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
  s <- g^2*(V - E)^2*xP.wg[,9]/nu*diff(times)[1] + s2
  
  dmu_s2 <- 0
  ds_dnu <- -(g^2*(V - E)^2*xP.wg[,9]/nu^2)*diff(times)[1]
  
  gnl_s2 <- as.numeric(sum(1/s) - t(y - mu) %*% ((y - mu)/s^2))*s2
  gnl_nu <- as.numeric(sum(ds_dnu/s)
                       -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_dnu*(y - mu))))*exp(tail(psi,1))
  
  return(c(gnl_s2,
           gnl_nu))
}