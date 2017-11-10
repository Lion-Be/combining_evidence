# ----------------------------------------------------------------------------------
# Combining Evidence for Inequality Constrained Hypotheses Using Multiple Studies
# Simulation Study - Defining Parameters
# Lion Behrens
# ----------------------------------------------------------------------------------

# Loading relevant packages (my code)
library(MASS)   # to call mvrnorm
library(Matrix) # to transform vcov matrix to nearest positive definite matrix
library(peperr) # to bootstrap
library(Bain)   # to construct BFs
library(dplyr)  # to shift y to first column in std function

# Loading relevant packages (Gorica code)
library(easypackages)
list.of.packages <- c("base", "blavaan", "boot", "coda", "foreign", "FRACTION", "lavaan", "lme4",
                      "MASS", "matrixcalc", "mvtnorm", "nlme", "quadprog", "R2OpenBUGS")
libraries(list.of.packages)

# Defining relevant functions for Gorica weights analysis
gorica <-
  function(object, ..., iter=100000){
    UseMethod("gorica")
  }


gorica_penalty <-
  function(object, iter=100000, mc.cores=1){
    
    if (!(inherits(object, "ormle") )) stop("object needs to be of class ormle")
    if (all(object$constr == 0) & object$nec == 0){
      penalty <- length(object$est)
    } else {
      if (iter < 1) stop("No of iterations < 1")
      
      est<-object$est 
      K<-length(est)
      covmtrx <- object$covmtrx
      constr<-object$constr
      rhs=object$rhs
      nec=object$nec
      
      Z <- rmvnorm(n=iter, mean=rep(0, K), sigma=covmtrx)
      Dmat2=2*ginv(covmtrx)
      
      nact <- apply(Z, 1, function(z){ 
        
        dvec2=2*(z%*%ginv(covmtrx)) 
        solveQP2= solve.QP(Dmat2,dvec2,t(constr),rhs,meq =nec,factorized = FALSE)
        if (solveQP2$iact[1] == 0) return(0) else return(length(solveQP2$iact))
      })
      
      dimsol <- K - nact
      LP <- sapply(1:K, function(x) sum(x == (dimsol)))/iter
      penalty <- sum((1:K)*LP[])
      
    }
    
    return(penalty)
    
  }



gorica.ormle <-
  function(object, ..., iter=100000){
    if (!inherits(object, "ormle") & !inherits(object, "list")) stop("object needs to be of class ormle or a list of ormle objects")
    if (iter < 1) stop("No of iterations < 1")
    if (inherits(object, "ormle")) objlist <- list(object, ...) else objlist <- object
    isorlm <- sapply(objlist, function(x) inherits(x, "ormle"))
    orlmlist <- objlist[isorlm]  
    Call <- match.call()
    Call$iter <- NULL
    if (inherits(object, "ormle")) names(orlmlist) <- as.character(Call[-1L])[isorlm]
    loglik <- -2*sapply(orlmlist, function(x) x$logLik)
    penalty <- 2*sapply(orlmlist, function(x) gorica_penalty(x, iter=iter))
    gorica <- loglik + penalty
    delta <- gorica - min(gorica)
    gorica_weights <- exp(-delta/2) / sum(exp(-delta/2))
    data.frame(misfit=loglik,complexity=penalty,gorica=gorica,gorica_weights=round(gorica_weights,4))
  }


gorica.list <- function(object, ..., iter=100000){
  if (all(sapply(object, class) == "ormle")) out <- gorica.ormle(object, iter=iter)
  return(out)
}


ormle <-
  function(est,covmtrx,constr,rhs,nec){
    K=length(est)
    covmtrx = as.matrix(covmtrx)
    Dmat = 2*ginv(covmtrx)
    dvec = 2*(est%*% ginv(covmtrx))
    solveQP = solve.QP(Dmat, dvec = dvec, t(constr), rhs, meq = nec, factorized = FALSE)
    tildeQ = solveQP$solution
    restrictedest=solveQP$solution
    names(restrictedest)=names(est)
    loglik =as.numeric( ( -K/2*log(2*pi) )-( 0.5*log(det(covmtrx) ) )-( 0.5* t(est- tildeQ)%*%ginv(covmtrx)%*% (est-tildeQ)) )
    
    out <- list(est=est, covmtrx=covmtrx, constr=constr, rhs=rhs, nec=nec, logLik=loglik,restrictedest=restrictedest)
    class(out) <- "ormle"
    return(out)
    
  }

print.ormle <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\n$est\n")
  print(x$est) 
  cat("\n")
  cat("\n$restrictedest\n")
  print(x$restrictedest) 
  cat("\n")
  invisible(x)
}

# Define standardization function
std <- function(data) {
  
  # Shift y to first column
  data <- data %>%
    select(y, everything())
  
  # Standardize all other except the first/dependent variable
  for (i in 2:length(data)) {
    data[,i] <- scale(data[,i])
  }
  as.data.frame(data)
}


# ---------------------------------- #
# Specifying simulation parameters   #
# ---------------------------------- #

# -------------
# Hypotheses
# -------------

n.hypos <- 2 # simulation code generalized for up to 4 hypotheses


  # Bain: Specifying hypotheses
  
    # H1: coef1>coef2>coef3>coef4 (true)
    ERr1<-matrix(0,0,0)  
    IRr1<-matrix(c(0,1,-1,0,0,0,
                   0,0,1,-1,0,0,
                   0,0,0,1,-1,0),nrow=3,ncol=6, byrow = TRUE)
    
    # H2: coef1>coef2>coef3<coef4 (slightly off)
    ERr2<-matrix(0,0,0)  
    IRr2<-matrix(c(0,1,-1,0,0,0,
                   0,0,1,-1,0,0,
                   0,0,0,-1,1,0),nrow=3,ncol=6, byrow = TRUE)
    
    # H3: coef1<coef2>coef3>coef4 (slightly off)
    ERr3<-matrix(0,0,0)  
    IRr3<-matrix(c(0,-1,1,0,0,0,
                   0,0,1,-1,0,0,
                   0,0,0,1,-1,0),nrow=3,ncol=6, byrow = TRUE)


  # Gorica: Specifying hypotheses
  
    # H1: coef1>coef2>coef3>coef4 (true)
    constr1<-matrix(c(0,1,-1,0,0,
                      0,0,1,-1,0,
                      0,0,0,1,-1),nrow=3,ncol=5, byrow = TRUE)
    rhs1 <- rep(0, 3) # right hand side values
    nec1 <- 0 # number of equality constrains
    
    # H2: coef1>coef2>coef3<coef4 (slightly off)
    constr2<-matrix(c(0,1,-1,0,0,
                      0,0,1,-1,0,
                      0,0,0,-1,1),nrow=3,ncol=5, byrow = TRUE)
    rhs2 <- rep(0, 3) # right hand side values
    nec2 <- 0 # number of equality constrains
    
    # H3: coef1<coef2>coef3>coef4 (slightly off)
    constr3<-matrix(c(0,-1,1,0,0,
                      0,0,1,-1,0,
                      0,0,0,1,-1),nrow=3,ncol=5, byrow = TRUE)
    rhs3 <- rep(0, 3) # right hand side values
    nec3 <- 0 # number of equality constrains
    
    # Hu: coef1, coef2, coef3, coef4, coef5 (unconstrained)
    constr <- matrix(c(rep(0, 5)), nrow = 1, ncol = 5, byrow = TRUE)
    rhs <- rep(0, 1)
    nec <- 0



# ---------------------------------------------
# Variance/Covariance Matrices in Populations
# ---------------------------------------------

# Population 1 (columns are X1, X2, X3, X4 without intercept)  
sigma.1 <- matrix(c(1,    0.25,  0.25,  0.25,   
                    0.25,   1,     0.25,  0.25,    
                    0.25,   0.25,  1,     0.25,     
                    0.25,   0.25,  0.25,   1),nrow=4,ncol=4)
mu.1 <- rep(0,4)

# Population 2  
sigma.2 <- matrix(c(1,    0.25,  0.25,  0.25,   
                    0.25,   1,     0.25,  0.25,    
                    0.25,   0.25,  1,     0.25,     
                    0.25,   0.25,  0.25,   1),nrow=4,ncol=4)
mu.2 <- rep(0,4)

# Population 3 
sigma.3 <- matrix(c(1,    0.25,  0.25,  0.25,   
                    0.25,   1,     0.25,  0.25,    
                    0.25,   0.25,  1,     0.25,     
                    0.25,   0.25,  0.25,   1),nrow=4,ncol=4)
mu.3 <- rep(0,4)

# Population 4 
sigma.4 <- matrix(c(1,    0.25,  0.25,  0.25,   
                    0.25,   1,     0.25,  0.25,    
                    0.25,   0.25,  1,     0.25,     
                    0.25,   0.25,  0.25,   1),nrow=4,ncol=4)
mu.4 <- rep(0,4)



# -------------------
# Regression models
# -------------------

model.1 <- "lm(y~V1+V2+V3+V4, data=boot.sample1)"
model.2 <- "lm(y~V1+V2+V3+V4, data=boot.sample2)"
model.3 <- "lm(y~V1+V2+V3+V4, data=boot.sample3)"
model.4 <- "lm(y~V1+V2+V3+V4, data=boot.sample4)"



# -----------------------
# Regression parameters
# -----------------------

# Population 1
n.coef1 <- 4 # without intercept
intercept.1 <- 0
ratios.1 <- c(4, 2, 1, 0.5) # ratios between true beta coefficients
f2.1 <- 0.35 # r2 = f2/(1+f2) 
ss.1 <- 500 # sample size

# Population 2
n.coef2 <- 4
intercept.2 <- 0
ratios.2 <- c(4, 2, 1, 0.5) # ratios between true beta coefficients
f2.2 <- 0.35 # r2 = f2/(1+f2) 
ss.2 <- 500 # sample size

# Population 3
n.coef3 <- 4
intercept.3 <- 0
ratios.3 <- c(4, 2, 1, 0.5) # ratios between true beta coefficients
f2.3 <- 0.35 # r2 = f2/(1+f2) 
ss.3 <- 500 # sample size

# Population 4
n.coef4 <- 4
intercept.4 <- 0
ratios.4 <- c(4, 2, 1, 0.5) # ratios between true beta coefficients
f2.4 <- 0.35 # r2 = f2/(1+f2) 
ss.4 <- 500 # sample size




# ------------------------------------
# Simulation and Bootstap Iterations
# ------------------------------------

gorica.iterations <- 10000

# Population 1
# Simulation Iterations
n.iter.1 <- 5
# Bootstrap Iterations
n.iter1 <- 100

# Population 2
# Simulation Iterations
n.iter.2 <- 5
# Bootstrap Iterations
n.iter2 <- 100

# Population 3
# Simulation Iterations
n.iter.3 <- 5
# Bootstrap Iterations
n.iter3 <- 100

# Population 4
# Simulation Iterations
n.iter.4 <- 5
# Bootstrap Iterations
n.iter4 <- 100


