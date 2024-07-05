#-----------------------------------------------------------------------------------
# On the relationship between life expectancy, modal age at death, and
# the threshold age of the life table entropy
#
#                                                                          
# Chiara MICHELETTI - Francisco VILLAVICENCIO 
#
# July 2024
#-----------------------------------------------------------------------------------

## Functions needed for the paper

## Function to compute the Exponential Integral
CalcE1 <- function(z, n = 150) {
  # Summation
  S <- (-1)^(1:n)*z^(1:n) / ((1:n)*factorial(1:n))
  S <- sum(S)
  # E1
  E1 <- digamma(1) - log(z) - S
  return(E1)
}

## Function phi(z)
CalcPhi <- function(z) {
  phi <- exp(z)*CalcE1(z)
  return(phi)
}

## Derivative phi'(z)
CalcPhiPrime <- function(z) {
  phi <- CalcPhi(z) - 1/z
  return(phi)
}

## Gompertz (modal) remaining life expectancy
CalcEx <- function(b, M, x) {
  # z = m(x)/beta
  z <- exp(b*(x - M))
  # Remaining life expectancy
  ex <- CalcPhi(z) / b
  # Output
  return(ex)
}

## Newton's method
NewtonInv <- function(y, tol = 1e-6, max_iter = 100) {
  
  # Initial guess
  if (y < 0.05) {
    z <- 40
  } else if (y < 0.1) {
    z <- 10
  } else if (y < 0.3) {
    z <- 5
  } else if (y < 0.5) {
    z <- 1.5
  } else if (y <= 1) {
    z <- .9
  } else if (y < 2) {
    z <- .2
  } else z <- 0.01
  # Iterations
  iter <- 0
  # Newton method
  while (abs(CalcPhi(z) - y) > tol & iter < max_iter) {
    zNew <- z - (CalcPhi(z) - y) / CalcPhiPrime(z)
    if (zNew > 0) {
      z <- zNew
      iter <- iter + 1
    } else break()
  }
  if (iter == max_iter) {
    warning("Maximum number of iterations reached")
  }
  # Output
  return(z)
  
}

## Function to deal with mortality rates at age 0
## From the paper "Dynamics of life expectancy and life span equality" (Aburto et al. 2020)
## Code available at https://github.com/jmaburto/Dynamics_Code
AKm02a0        <- function(m0, sex = "m"){
  sex <- rep(sex, length(m0))
  ifelse(sex == "m", 
         ifelse(m0 < .0230, {0.14929 - 1.99545 * m0},
                ifelse(m0 < 0.08307, {0.02832 + 3.26201 * m0},.29915)),
         # f
         ifelse(m0 < 0.01724, {0.14903 - 2.05527 * m0},
                ifelse(m0 < 0.06891, {0.04667 + 3.88089 * m0}, 0.31411))
  )
}


# Function to compute h := -log(\Bar{H}), the transformation of the life table entropy 
# (to get equality instead of inequality), from the the life table death rates mx
## From the paper "Dynamics of life expectancy and life span equality" (Aburto et al. 2020)
## Code available at https://github.com/jmaburto/Dynamics_Code
h.frommx <- function(mx, sex){
  i.openage <- length(mx)
  OPENAGE   <- i.openage - 1
  RADIX     <- 1
  ax        <- mx * 0 + .5
  ax[1]     <- AKm02a0(m0 = mx[1], sex = sex)
  qx        <- mx / (1 + (1 - ax) * mx)
  qx[i.openage]       <- ifelse(is.na(qx[i.openage]), NA, 1)
  ax[i.openage]       <- 1 / mx[i.openage]                   
  px 				    <- 1 - qx
  px[is.nan(px)]      <- 0
  lx 			        <- c(RADIX, RADIX * cumprod(px[1:OPENAGE]))
  dx 				    <- lx * qx
  Lx 				    <- lx - (1 - ax) * dx
  Lx[i.openage ]	    <- lx[i.openage ] * ax[i.openage ]
  Tx 				    <- c(rev(cumsum(rev(Lx[1:OPENAGE]))),0) + Lx[i.openage]
  ex 				    <- Tx / lx
  l <- length(ex)
  v <- (sum(dx[-l]* (ex[-l] + ax[-l]*(ex[-1]-ex[-l]) )) + ex[l])
  k <- v/ex[1]
  eq <- -log(k)
  return(eq)
}

## Function to compute the threshold age of the life table entropy for every country-year
## From the paper "Dynamics of life expectancy and life span equality" (Aburto et al. 2020)
## Code available at https://github.com/jmaburto/Dynamics_Code
threshold.decomp.function <- function(DT = .SD, func2 = h.frommx, Age = 20:100){
  mat        <- acast(DT, Age~Year, value.var="mx")  # matrice age x years
  count      <- DT$PopName2[1]
  Sex        <- DT$Sex1[1]
  
  registerDoParallel(4)
  threshold.ages     <- unlist(foreach(i=1:dim(mat)[2],.packages = 'DemoDecomp',.export = 
                                         c('AKm02a0')) %dopar% {
                                           if (as.integer(colnames(mat)[i]) >= 1960){Age <- 45:100}
                                           if (count == 'ISL' & Sex == 'Female' & as.integer(colnames(mat)[i]) %in% c(1981,1982,1991,2016)){Age <- 70:100}
                                           perturbation <- horiuchi(func=func2, pars1=mat[,i] , pars2=mat[,i]*.99, N=35,sex=DT$Sex[1])
                                           f <- approxfun(Age,perturbation[(Age+1)],method = "linear",rule=2 )
                                           tryCatch({
                                             a <- uniroot(function(x) f(x)-0,interval = c(Age[1],Age[length(Age)]))$root
                                             a},
                                             error=function(x){
                                               a <- -1
                                               a})})
  stopImplicitCluster()
  y          <- do.call(rbind, lapply(threshold.ages, data.frame, stringsAsFactors=FALSE))
  y
}

## Function to compute the Gompertz parameters alpha and beta
## Age range default to 30:90
param.gomp <- function(DT = .SD, Age = 30:90){
  
  mat.dx <- acast(DT, Age~Year, value.var="Deaths")[Age+1,] # +1 because we start from age 0
  mat.ex <- acast(DT, Age~Year, value.var="Exposures")[Age+1,]
  
  registerDoParallel(4)
  
  par.gomp     <- matrix(unlist(foreach(i=1:dim(mat.ex)[2]) %dopar% {
    #if ( as.integer(colnames(mat.ex)[i]) >= 1960){Age <- 45:100}
    
    X <- cbind(1, Age)
    Off <- log(mat.ex[,i])
    
    a <- coefficients(glm(mat.dx[,i] ~  -1 + X + offset(Off),
                        family = poisson(link = "log")))[1]
    
    b <- coefficients(glm(mat.dx[,i] ~  -1 + X + offset(Off),
                        family = poisson(link = "log")))[2]
    
    pars <- c(exp(a),b)
    
    pars}), byrow = T, ncol=2)
  
  stopImplicitCluster()
  y              <- data.frame(par.gomp, stringsAsFactors = FALSE)
  names(y)       <- c('Alpha', 'Beta')
  row.names(y)   <- colnames(mat.dx)
  y
}


## Function to compute the Gompertz mode from alpha and beta
gompertz.modal.age <- function(alpha, beta){
  M <- (log(beta/alpha))/beta
  return(M)
}

## Function to compute the approximation of the threshold age of life table entropy in terms of the mode
approx.ah <- function(M, beta){
  gamma<- -digamma(1) # Euler-Mascheroni constant
  M - gamma / beta
}

## Function to compute the mode of a density
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

#-----------------------------------------------------------------------------------

## GGplot aesthetics
custom_theme <- function() {
  theme_bw() +
    theme(
      axis.text = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      plot.title = element_text(size = 18, hjust = .1),
      axis.title.y = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.text = element_text(color = 'black'),
      axis.line.x = element_line(color = "black", size = .5),
      axis.line.y = element_line(color = "black", size = .5),
      axis.ticks = element_line(colour = "black", size = .8),
      axis.ticks.length=unit(0.25, "cm"),
      strip.placement = 'outside',
      strip.background = element_rect(colour = "white", fill = "white"),
      legend.position = 'none'
    )
}
