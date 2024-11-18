n2c.nl.comp <- function(psi, mu, s2, V, E){
  i <- psi[1]
  N <- 1 + exp(psi[2])

  s2.fit <- i*mu*(V-E) - mu^2/N
  nl <- as.numeric(-2*t(s2.fit)%*%s2 + t(s2.fit)%*%s2.fit)
  return(nl)
}

n2c.gnl.comp <- function(psi, mu, s2, V, E){
  i <- psi[1]
  N <- 1 + exp(psi[2])

  s2.fit <- i*mu*(V-E) - mu^2/N

  gnl <- rep(NA,2)
  gnl[1] <- as.numeric(-2*t(mu*(V-E))%*%s2 + 2*t(s2.fit)%*%(mu*(V-E)))
  gnl[2] <- as.numeric(-2*t(mu^2/N^2)%*%s2 + 2*t(s2.fit)%*%mu^2/N^2)*exp(psi[2])

  return(gnl)
}

