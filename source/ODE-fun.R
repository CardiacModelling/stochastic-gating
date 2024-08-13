
#' @keywords internal
generate.lambda <- function(rct.lst,
                            constr.lst = NULL,
                            negSign = TRUE,
                            envir = .GlobalEnv,
                            flickering = FALSE){
  
  if(!is.null(constr.lst)){
    free_fixed <- t(sapply(constr.lst, function(constr){
      constr <- as.vector(unlist(strsplit(constr, split = "=", fixed = TRUE)))
      pattern <- gsub(pattern = "\\\\\\'\\\\\\]",
                      replacement = "",
                      x = gsub(pattern = "lambda\\\\\\[\\\\\\'",
                               replacement = "",
                               x = constr[1]))
      replacement <- gsub(pattern = "\\\\\\'\\\\\\]",
                          replacement = "",
                          x = gsub(pattern = "lambda\\\\\\[\\\\\\'",
                                   replacement = "",
                                   x = constr[2]))
      
      return(c(pattern, replacement))
    }, simplify = "array"))
    
    rct.lst <- rct.lst[c(which(rct.lst %in% free_fixed[,2]),
                         which(!(rct.lst %in% c(free_fixed[,2], free_fixed[,1]))),
                         which(rct.lst %in% free_fixed[,1]))]
  }
  
  nLambda <- length(rct.lst)
  js <- which(1:(2*nLambda) %% 4 == 1)
  
  l_string <- c(sapply(js, function(j){
    c(paste0("theta[", j, "]*exp(theta[", j+1, "]*V)"),
      paste0("theta[", j+2, "]*exp(",ifelse(negSign,"-",""),"theta[", j+3, "]*V)"))
  }))
  names(l_string) <- rct.lst
  
  
  for (constr in constr.lst) {
    constr <- as.vector(unlist(strsplit(constr, split = "=", fixed = TRUE)))
    pattern <- gsub(pattern = "\\\\\\'\\\\\\]",
                    replacement = "",
                    x = gsub(pattern = "lambda\\\\\\[\\\\\\'",
                             replacement = "",
                             x = constr[1]))
    replacement <- gsub(pattern = "\\\\\\'\\\\\\]",
                        replacement = "",
                        x = gsub(pattern = "lambda\\\\\\[\\\\\\'",
                                 replacement = "",
                                 x = constr[2]))
    
    l_string[pattern] <- l_string[replacement]
  }
  
  
  if(flickering){
    ntheta <- 2*length(unique(l_string))
    wToOpen <- which(sapply(names(l_string), function(name_i){unlist(strsplit(name_i, split = "->", fixed = TRUE))[2]}) == "O")
    wToFlickering <- which(sapply(names(l_string), function(name_i){unlist(strsplit(name_i, split = "->", fixed = TRUE))[2]}) == "F")
    
    l_string[wToOpen] <- paste(paste0("theta[",ntheta + 1,"]"), l_string[wToOpen], sep = "*")
    l_string[wToFlickering] <- paste(paste0("(1 - theta[",ntheta + 1,"])"), l_string[wToFlickering], sep = "*")
    l_string <- c(l_string, "F->O" = paste0("theta[",ntheta + 1,"]"), "O->F" = paste0("(1 - theta[",ntheta + 1,"])"))
    
  }
  
  l_string <- sapply(1:length(l_string), function(r){
    paste0(paste0("'", names(l_string)[r],"'"), " = ", l_string[r])
  })
  
  l_string <- paste(l_string, collapse = ",\n")
  
  get.lambda.string <- paste0("get.lambda <- function(theta, V){\nc(",
                              l_string,
                              ")\n}")
  
  eval(parse(text=get.lambda.string), envir = envir)
}


#' @keywords internal
generate.get.Vlambda <- function(ct.lst,
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
  
  Vth <- Vth[-ncol(Vth),]
  
  Vth <- t(sapply(1:nrow(Vth), function(i){
    sapply(1:(ncol(Vth)-1), function(j){
      paste(ifelse(Vth[i,j] == "0", ifelse(Vth[i,ncol(Vth)] == "0", "0", ""), Vth[i,j]),
            ifelse(Vth[i,ncol(Vth)] == "0", "", Vth[i,ncol(Vth)]),
            sep = ifelse(Vth[i,ncol(Vth)] == "0", "", " -"))
    })
  }))
  
  get.Vlambda.string <- paste("get.Vlambda <- function(lambda){
  matrix(data = c(",
                              paste(as.vector(unlist(lapply(split(as.vector(t(Vth)), ceiling(seq_along(1:(nrow(Vth)^2))/nrow(Vth))), function(lst){
                                paste(lst, collapse = ",")
                              }))), collapse = ",\n"),
                              "),\n byrow = TRUE, nrow = ", nrow(V) - 1, ", ncol = ", nrow(V) - 1,
                              ")\n}", sep = "")
  
  eval(parse(text=get.Vlambda.string), envir = envir)
}


generate.get.dx <- function(ct.lst,
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
  
  get.dx.string <- paste0("get.dx <- function(t = 0,
                            y,
                            parms){\n
  lambda <- get.lambda(theta = parms,
                       V = get.voltage(t = t))\n
  Vl <- get.Vlambda(lambda)\n",
                          "dx <- Vl %*% y + c(",paste(b, collapse = ",\n"),")\n",
                          "dx <- list(dx = dx[,1])\n
  return(dx)\n}")
  
  eval(parse(text=get.dx.string), envir = envir)
}

generate.get.dx.wg <- function(n,
                               p,
                               ct.lst,
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
  
  
  get.dx.wg.string <- paste0("get.dx.wg <- function(t = 0,
                      y,
                      parms){\n",
                             "x <- head(y, ", n, ")\n",
                             "s <- tail(y, ", n*p, ")\n",
                             "lambda <- get.lambda(theta = parms,\n V = get.voltage(t = t))\n",
                             "Vl <- get.Vlambda(lambda)\n",
                             "dx <- list(dx = c(Vl %*% x + c(",
                             paste(b, collapse = ",\n"), "),\n",
                             "c(Vl %*% matrix(s, ncol = ", p, ")) + get.dx_dj(t = t,
                                                              y = x,
                                                              parms = parms)))\n",
                             "return(dx)\n}")
  
  eval(parse(text=get.dx.wg.string), envir = envir)
}

generate.get.x.steady <- function(p,
                                  ct.lst,
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
  
  get.x.steady.string <- paste0("get.x.steady <- function(t = 0,
                         parms){\n",
                                "lambda <- get.lambda(theta = parms,\n V = get.voltage(t = t))\n",
                                "Vl <- get.Vlambda(lambda)\n",
                                "x.s <- (-solve(Vl)%*%c(",
                                paste(b, collapse = ",\n"),
                                "))[,1]\n",
                                "names(x.s) <- c(",
                                paste(sapply(head(ct.lst,-1), function(s){paste0("'", s, "'")}), collapse = ",\n"), ")\n",
                                "return(x.s)\n}")
  
  eval(parse(text=get.x.steady.string), envir = envir)
  
}

generate.get.x.steady.forgrad <- function(p,
                                          ct.lst,
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
  
  get.x.steady.forgrad.string <- paste0("get.x.steady.forgrad <- function(",
                                        paste(paste0("p", 1:p), collapse = ",\n"),
                                        ",\nt){\n",
                                        "parms <- c(",
                                        paste(paste0("p", 1:p), collapse = ",\n"),
                                        ")\n",
                                        "lambda <- get.lambda(theta = parms,\n V = get.voltage(t = t))\n",
                                        "Vl <- get.Vlambda(lambda)\n",
                                        "(-solve(Vl)%*%c(",
                                        paste(b, collapse = ",\n"),
                                        "))[,1]\n}")
  
  eval(parse(text=get.x.steady.forgrad.string), envir = envir)
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
                              cache.exp = TRUE)
  
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

generate.get.Jf.wg <- function(n,
                               envir = .GlobalEnv){
  
  
  get.Jf.wg.string <- paste0("get.Jf.wg <- function(t = 0,
                      y,
                      parms){\n",
                             "x <- head(y, ", n, ")\n",
                             "Vl <- get.Vlambda(get.lambda(theta = parms,\n V = get.voltage(t = t)))\n",
                             "Jfull[bdiag.idxs] <- Vl\n
  Jfull[Jh.idxs] <- get.dx_dj_di(t = t,
                                 y = x,
                                 parms = parms)\n",
                             "return(Jfull)\n}")
  
  eval(parse(text=get.Jf.wg.string), envir = envir)
  
}

generate.get.dx.forgrad <- function(p,
                                    n,
                                    ct.lst,
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
  
  fun.string <- paste0("get.dx.forgrad <- function(",
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
                       "lambda <- get.lambda(theta = parms,
                       V = get.voltage(t = t))
  Vl <- get.Vlambda(lambda)
  c(Vl %*% y + matrix(data = c(",paste(b, collapse = ",\n"),"), ncol = 1))\n}")
  
  eval(parse(text=fun.string), envir = envir)
}

get.V <- function(ct.lst = ct.lst, rct.lst = rct.lst){
  V <- matrix(data = 0, nrow = length(ct.lst), ncol = length(rct.lst))
  rownames(V) <- ct.lst
  colnames(V) <- rct.lst
  
  for (r in rct.lst) {
    rgts_prod <- as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))
    V[rgts_prod[1],r] <- -1
    V[rgts_prod[2],r] <- 1
  }
  return(V)
}

get.voltage <- function(t){
  if(t >= 0 & t <= 250){
    v <- -80
  }
  if(t > 250 & t <= 300){
    v <- -120
  }
  if(t > 300 & t <= 500){
    v <- -80
  }
  if(t > 500 & t <= 1500){
    v <- 40
  }
  if(t > 1500 & t <= 2000){
    v <- -120
  }
  if(t > 2000 & t <= 3000){
    v <- -80
  }
  if(t > 3000 & t <= 6500){
    v <- get.sinewave(t = t,
                      t0 = 2500)
  }
  if(t > 6500 & t <= 7000){
    v <- -120
  }
  
  if(t > 7000){
    v <- -80
  }
  
  v
}

get.sinewave <- function(t, t0){
  
  A1 <- 54
  A2 <- 26
  A3 <- 10
  
  w1 <- 0.007
  w2 <- 0.037
  w3 <- 0.19
  
  v <- (-30 + A1*sin(w1*(t-t0))
        +  A2*sin(w2*(t-t0))
        +  A3*sin(w3*(t-t0)))
  
  return(v)
}

get.Jf <- function(t = 0,
                   y,
                   parms){
  
  Vl <- get.Vlambda(get.lambda(theta = parms,
                               V = get.voltage(t = t)))
  
  return(Vl)
}

generate.bdiagmat <- function(n,
                              p,
                              envir = .GlobalEnv){
  
  bdiag.idxs.string <- paste0("bdiag.idxs <- cbind(c(matrix(1:",
                              n*(p+1),", nrow = ", n, ")[,rep(1:", p+1, ", each = ", n, ")]),\n",
                              "rep(1:", n*(p+1), ", each = ", n, "))")
  
  Jh.idxs.string <- paste0("Jh.idxs <- cbind(rep(",
                           n+1, ":", n*(p+1), ", times = ", n, "),\n",
                           "rep(", 1, ":", n, ", each = ", n*((p+1)-1), "))")
  
  Jfull.string <- paste0("Jfull <- matrix(", 0, ", ", n*(p+1), ", ", n*(p+1), ")")
  
  
  eval(parse(text=bdiag.idxs.string), envir = envir)
  eval(parse(text=Jh.idxs.string), envir = envir)
  eval(parse(text=Jfull.string), envir = envir)
  
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
                      cache.exp = TRUE)
  
  get.dx_djs_dxis_cbind <- Deriv(f = get.dx_djs,
                                 x = paste0("x", 1:n),
                                 combine = "cbind",
                                 cache.exp = TRUE)
  
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
                            cache.exp = TRUE)
  
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

generate.nl_try <- function(p,
                            jsens,
                            envir = .GlobalEnv){
  
  
  nl_try.string <- paste0("nl_try <- function(psi,\n",
                          ifelse(length(setdiff(1:p, jsens)) > 0, paste0(paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"), ",\n"), ""),
                          "y,
                   times,
                   V,
                   E,
                   hmax){
  nl.res <- try(nl(psi = psi,\n",
                          ifelse(length(setdiff(1:p, jsens)) > 0, paste0(paste(paste0("p",setdiff(1:p, jsens), " = p", setdiff(1:p, jsens)), collapse = ",\n"), ",\n"), ""),
                          "y = y,
                   times = times,
                   V = V,
                   E = E,
                   hmax = hmax), silent = TRUE)
  if(inherits(nl.res,'try-error') |
     is.nan(nl.res) |
     is.infinite(nl.res)){
    return(1e8)
  }
  return(nl.res)
}")
  
  eval(parse(text=nl_try.string), envir = envir)
  
}

generate.gnl_try <- function(p,
                             jsens,
                             envir = .GlobalEnv){
  
  
  gnl_try.string <- paste0("gnl_try <- function(psi,\n",
                           ifelse(length(setdiff(1:p, jsens)) > 0, paste0(paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"), ",\n"), ""),
                           "y,
                   times,
                   V,
                   E,
                   hmax){
  gnl.res <- try(gnl(psi = psi,\n",
                           ifelse(length(setdiff(1:p, jsens)) > 0, paste0(paste(paste0("p",setdiff(1:p, jsens), " = p", setdiff(1:p, jsens)), collapse = ",\n"), ",\n"), ""),
                           "y = y,
                   times = times,
                   V = V,
                   E = E,
                   hmax = hmax), silent = TRUE)
  if(inherits(gnl.res,'try-error') |
     sum(is.nan(gnl.res)) > 0 |
     sum(is.infinite(gnl.res)) > 0){
    return(rep(1e8, length(psi)))
  }
  return(gnl.res)
}")
  
  eval(parse(text=gnl_try.string), envir = envir)
  
}


get.logistic <- function(x){
  return(exp(x)/(exp(x) + 1))
}

get.dx.logistic <- function(x){
  return(exp(x)/(exp(x) + 1)^2)
}

get.logit <- function(p){
  return(log(p/(1-p)))
}

generate.nl <- function(p,
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


generate.gnl <- function(p,
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

generate.fim <- function(p,
                         jsens,
                         envir = .GlobalEnv){
  
  fim.string <- paste0("fim <- function(psi,\n",
                       # paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"),
                       ifelse(length(setdiff(1:p, jsens)) > 0, paste0(paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"), ",\n"), ""),
                       #",\ny,
                       "y,
                   times,
                   V,
                   E,
                   hmax){
                   phi <- exp(psi)
                   theta <- head(phi, -1)
                   g <- tail(phi, 1)

                     x <- try(vode(y = c(get.x.steady(t = max(times),\n",
                       # "parms = c(theta,",
                       "parms = ",
                       # paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"),
                       ifelse(length(setdiff(1:p, jsens)) > 0, paste0("c(theta,\n",
                                                                      paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"),
                                                                      ")"), "theta"),
                       "),
                      get.x.steady_djs(t = max(times),\n",
                       # "parms = c(theta,",
                       "parms = ",
                       # paste(paste0("p",setdiff(1:p, jsens)), collapse = ",\n"),
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
                dmu <- cbind(dmu_th, dmu_g)

                dmu <- dmu * matrix(data = exp(psi),
                      ncol = length(psi),
                      nrow = nrow(dmu),
                      byrow = TRUE)

                return(t(dmu) %*% dmu)\n}")
  
  eval(parse(text=fim.string), envir = envir)
  
}

