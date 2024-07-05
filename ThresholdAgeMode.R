#------------------------------------------------------------------------------#
# On the relationship between life expectancy, modal age at death, and
# the threshold age of the life table entropy
#
#                                                                          
# Chiara MICHELETTI - Francisco VILLAVICENCIO 
#
# July 2024
#------------------------------------------------------------------------------#

#----------#
# PACKAGES #                      
#----------#

# Clear working directory
rm(list = ls())

# Libraries
suppressPackageStartupMessages({
  # install.packages('viridisLite')
  library(viridisLite)
  # install.packages('data.table')
  library(data.table)
  # install.packages('foreach')
  library(foreach)
  # install.packages('doParallel')
  library(doParallel)
  # install.packages('DemoDecomp')
  library(DemoDecomp)
  # install.packages('tidyverse')
  library(tidyverse)
  # install.packages('reshape2')
  library(reshape2)
  # install.packages('patchwork')
  library(patchwork)
})

# Avoid scientific notation
options(scipen = 999)

# Necessary functions
source('Functions.R')


#------------------------------------------------------------------------------#
#                               PROOF OF (1) 
#------------------------------------------------------------------------------#

# Modal age at death
M <- 66.06 

# Gompertz rate of ageing
b <- 0.061 

# Multiplier
round(exp(exp(-b*M)), 4)

# Summation
round(-(CalcE1(z = exp(-b*M)) - digamma(1) - b*M), 4)

# Approximation
M + digamma(1)/b
CalcEx(b = b, M = M, x = 0)


#------------------------------------------------------------------------------#
#                   FIGURE 1: LIFE EXPECTANCY AND MODE                      
#------------------------------------------------------------------------------#

#------------#
# PARAMETERS #
#------------#

# Modal age at death
M <- seq(30, 100, .05)

# Gompertz rate of ageing
b <- seq(.05, .16, by = 0.001)


#------#
# PLOT #
#------#

# Categories
cuts <- c(0, 0.05, 0.1, .5, 1, 10)

# Color palette
colPalette <- rev(viridis(n = length(cuts) - 1, option = 'viridis'))

# Output file
PDF <- F
if (PDF) pdf('Figure1.pdf', width = 7, height = 6.5, family = 'Helvetica')

# LAYOUT
layout(cbind(c(1, 0, 4), c(3, 2, 4)),
       widths = c(0.11, 1), heights = c(1, 0.08, 0.07))

# Y-axis
par(mar = rep(0, 4))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = '', ylab = '')
text(0.3, 0.5, expression(paste('Gompertz rate of ageing ', beta)), cex = 1.8, srt = 90)

# X-aixs
par(mar = rep(0, 4))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = '', ylab = '')
text(0.5, 0.4, 'Gompertz modal age at death M', cex = 1.8)

param <- cbind(rep(M, each = length(b)), rep(b, length(M)))
dat <- apply(param, 1, function(y) {
  e0 <- CalcEx(b = y[2], M = y[1], x = 0)
  er <- abs(e0 - y[1] - digamma(1)/y[2])
  return(er)
})

# Store results
results <- as.data.frame(cbind(param, dat))
names(results) <- c('M', 'beta', 'E')

# Re-shape data for plotting
dat <- matrix(dat, nrow = length(M), byrow = T)

# Plot
par(las = 1, mar = c(1.5, 1.5, 1.5, 1.5))
image(dat, breaks = cuts, col = colPalette, axes = F)

# X-axis labels
ticks <- seq(min(M), max(M), 10)
axis(1, at = seq(1/length(M), 1 - 1/length(M), length.out = length(ticks)),
     labels = seq(min(M), max(M), 10), cex.axis = 1.2)

# Y-axis labels
ticks <- seq(0, 1, length.out = 1 + round(100*(max(b) - min(b))))
axis(2, at = ticks,
     labels = format(seq(min(b), max(b), length.out = length(ticks)), nsmall = 2),  
     cex.axis = 1.2)

# Legend
par(mar = c(0, 1, 0, .5))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = '', ylab = '')
textwidths <- c(0, .15, .35, .53, .7) / 0:(length(cuts) - 2)
textwidths[1] <- 0

legend('left',
       legend = sapply(c(bquote(epsilon < .(cuts[2])),
                         bquote(.(cuts[2]) < {epsilon < .(cuts[3])}),
                         bquote(.(cuts[3]) < {epsilon < .(cuts[4])}),
                         bquote(.(cuts[4]) < {epsilon < .(cuts[5])}),
                         bquote(epsilon > .(cuts[5]))),
                       as.expression),
       border = NA, fill = colPalette, bty = 'n', cex = 1.5,
       ncol = length(cuts) - 1, text.width = textwidths)

# Close device
if (PDF) dev.off()


#-------------------#
# IN-TEXT ESTIMATES #                      
#-------------------#

round(100*length(which(results$E < 0.5)) / nrow(results))
round(100*length(which(results$E > 1)) / nrow(results))
p <- which(results$M >= 60 & results$beta >= 0.08)
round(100*length(which(results$E[p] < 0.1)) / length(p))
round(max(results$E[p]), 2)
round(CalcEx(b = 0.08, M = 60, x = 0), 1)

# Remove everything but functions
rm(list = setdiff(ls(), lsf.str()))


#------------------------------------------------------------------------------#
#                              PROOF OF (2a) 
#------------------------------------------------------------------------------#

# Modal age
M <- 30:100

# Rate of ageing
b <- 5:16 / 100


#--------------------#
# DOMAIN OF PHI^(-1) #
#--------------------#

# Data frame with Gompertz parameters
dat <- data.frame(M = rep(M, each = length(b)), beta = rep(b, length(M)))

# Life expectancy at birth
dat$e0 <- apply(dat[, c(1, 2)], 1, 
                function(x) {
                  ex <- CalcEx(b = x[2], M = x[1], x = 0)
                  return(ex)
                })

# Quotient w = beta*e_o / (beta*e_o + 1) 
dat$w <- dat$beta*dat$e0 / (dat$beta*dat$e0 + 1)
floor(1000*range(dat$w)[1])/1000                  
ceiling(1000*range(dat$w)[2])/1000                  

# Restrict to modern mortality settings
dat2 <- dat[dat$beta >= 0.08 & dat$M >= 60, ]
floor(1000*range(dat2$w)[1])/1000                  
ceiling(1000*range(dat2$w)[2])/1000
rm(dat2)


#-----------------#
# NEWTON'S METHOD #
#-----------------#

# Generate input values for y
input_values <- dat$w[which(dat$w >= 0.81 & dat$w <= 0.94)]

# Calculate inverses using Newton's method
inverses <- sapply(input_values, function(y) NewtonInv(y))

# Range of -log (Equation 7)
floor(1000*range(-log(inverses))[1])/1000                  
ceiling(1000*range(-log(inverses))[2])/1000  

# Maximum error
ceiling(1000*max(abs(-digamma(1)+log(inverses))))/1000


#-----------------------------------#
# TABLE A-1 (APPENDIX): Test EQ.(7) #
#-----------------------------------#

# Quotient w = e0 / (beta*e0 + 1)
dat$w <- dat$beta*dat$e0 / (dat$beta*dat$e0 + 1)

# Keep only those corresponding to beta>= 0.08 and M>=60
dat <- dat[which(dat$w >= .81 & dat$w <= .94), ]

# Calculate inverse of \phi in w
dat$inverse <- sapply(dat$w, function(y) NewtonInv(y))

# -log of the inverse
dat$log <- -log(dat$inverse)

# Tidy up
dat$e0       <- as.numeric(format(round(dat$e0, 2), nsmall = 2))
dat$w        <- as.numeric(format(round(dat$w, 3), nsmall = 3))
dat$inverse  <- as.numeric(format(round(dat$inverse, 3), nsmall = 3))
dat$log      <- as.numeric(format(round(dat$log, 3), nsmall = 3))
dat$err      <- as.numeric(format(round(abs(digamma(1) + dat$log), 3), nsmall = 3))

# Min and maximum absolute error
range(dat$err)

# Table A-1
dat <- dat[which(dat$M %in% seq(40,95,5) & abs(dat$log + digamma(1)) < 0.03), ]
dat

# Remove everything but functions
rm(list = setdiff(ls(), lsf.str()))


#------------------------------------------------------------------------------#
#          FIGURE 2: THRESHOLD AGE, MODE AND LIFE EXPECTANCY                      
#------------------------------------------------------------------------------#

#------------#
# PARAMETERS #
#------------#

# Ages
x <- 0:110

# Modal ages
M <- 30:100

# Rate of ageing beta
b <- seq(.05, .16, by = 0.01)

# Data frame to store results
results <- data.frame()


#------#
# PLOT #
#------#

# Colours
col.eo  <- "#ffbb44" 
col.m   <- "#859b6c" 
col.a.h <- "#004f63" 

# Output file
PDF <- F
if (PDF) pdf('Figure2.pdf', width = 10, height = 12, family = 'Helvetica')

# LAYOUT
layout(cbind(c(rep(1, 4), 0, 15),
             c(3, 6, 9, 12, 2, 15),
             c(4, 7, 10, 13, 2, 15),
             c(5, 8, 11, 14, 2, 15)),
       widths = c(0.2, rep(1, 3)), heights = c(rep(1, 4), 0.25, 0.15))

# Y-axis
par(mar = rep(0, 4))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = '', ylab = '')
text(0.4, 0.5, 'Years (age)', cex = 2.2, srt = 90)

# X-aixs
par(mar = rep(0, 4))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = '', ylab = '')
text(0.5, 0.35, 'Gompertz modal age at death M', cex = 2)

# Plot
for (i in 1:length(b)) {
  
  # THRESHOLD AGE
  param <- cbind(rep(M, each = length(x)), rep(x, length(M)))
  dat <- apply(param, 1, function(y) {
    if (abs(y[2] - y[1]) <= 20) {
      ex <- CalcEx(b = b[i], M = y[1], x = y[2])
      e0 <- CalcEx(b = b[i], M = y[1], x = 0)
      out <- ex - (e0 / (b[i]*e0 + 1))
    } else out <- NA
    return(out)
  })
  
  # Approximate with SPLINES
  dat <- cbind(param, dat)
  thresh <- tapply(dat[, 3], INDEX = dat[, 1], function(y) {
    age <- which(!is.na(y))
    fit <- spline(x = age - 1, y = y[age], n = 10000)
    xi <- fit$x 
    yi <- fit$y 
    xn <- xi[which.min(abs(yi))] 
    return(xn)
  })
  
  # LIFE EXPECTANCY
  e0 <- sapply(M, function(y) {
    e0 <- CalcEx(b = b[i], M = y, x = 0)
    return(e0)
  })
  
  # Store results
  results <- rbind(results, 
                   data.frame(beta = rep(b[i], length(M)), M = M, 
                              e0 = e0, aH = thresh))
  
  # Y-axis
  yAxis <- F
  
  # X-axis
  if (i %in% 10:12) {
    xAxis <- T
  } else xAxis <- F
  
  # Margins
  if (i %in% c(1, 4, 7, 10)) {
    par(las = 1, mar = c(1.5, 2.2, 1.8, 0))  
    yAxis <- T
  } else if (i %in% c(2, 5, 8, 11)) {
    par(las = 1, mar = c(1.5, 1.6, 1.8, 0.6))  
  } else {
    par(las = 1, mar = c(1.5, 1, 1.8, 1.2))
  }
  
  # Plot
  plot(NA, xlim = range(M), ylim = range(M), type = 'l', lty = 2, 
       xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  
  # Modal age at death
  abline(a = 0, b = 1, col = col.m, lwd = 3, lty = 3)
  
  # Threshold age
  lines(M, thresh, col = col.a.h, lwd = 3, lty = 1)
  
  # Life expectancy at birth
  lines(M, e0, col = col.eo, lwd = 3, lty = 4)
  
  # Title
  title(bquote(beta~'='~.(format(b[i], nsmall = 2))), adj = 0.02, line = .8, cex.main = 1.7)
  
  # Y-axis labels
  if (yAxis) {
    axis(2, at = seq(min(M), max(M), 10), cex.axis = 1.2)
  } else {
    axis(2, at = seq(min(M), max(M), 10), labels = NA)
  }
  
  # X-axis labels
  if (xAxis) {
    axis(1, at = seq(min(M), max(M), 10), cex.axis = 1.2)
    axis(1, at = max(M), cex.axis = 1.2, tick = F)
  } else axis(1, at = seq(min(M), max(M), 10), labels = NA)
  
}

# Legend
par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = '', ylab = '')
legend('center', 
       legend = c('Modal age at death', 'Threshold age', 'Life expectancy at birth'),
       col = c(col.m, col.a.h, col.eo), lwd = c(2.5, 2.5, 2.5), lty = c(3, 1, 4),
       seg.len = 3, ncol = 3, bty = 'n', cex = 2, text.width = c(.2,.15, .2))

# Close device
if (PDF) dev.off()

# Remove everything but functions
rm(list = setdiff(ls(), lsf.str()))


#------------------------------------------------------------------------------#
#  FIGURE 3: THRESHOLD AGE, MODE AND LIFE EXPECTANCY SWEDISH FEMALES 1900-2019                     
#------------------------------------------------------------------------------#

# Swedish female life tables from 1900 to 2021 - Human Mortality Database
load("Data/Sweden_female_life_tables.RData")

# NOTE: The modal age at death (column 'mode') was calculated using the
#  'MortalitySmooth' R package. This package is no longer available in CRAN,
#  but can be downloaded from Tim Riffe's GitHub using the 'devtools' R package.
#  https://github.com/timriffe/MortalitySmooth

# Females
df.LT <- LT.SWE
df.LT$Sex1 <- df.LT$Sex
df.LT$PopName2 <- df.LT$PopName

# Compute the threshold ages for Sweden
Threshold.ages <- df.LT[, c(a.h = threshold.decomp.function(.SD, func2 = h.frommx)), by = list(Sex, PopName)]

# Modify structure of the dataset to add the year
names(Threshold.ages) <- c('Sex', 'PopName','a.h')
year.label <- df.LT[, list(Year = unique(Year)), by = list(Sex,PopName)]
df.a.h <- cbind(Threshold.ages, year.label[,3])
rm(year.label, df.LT, Threshold.ages)

# Extract life expectancy at birth and modal age at death (max of the age-at-death distribution)
df.eo <- LT.SWE[Age == 0, ]
rm(LT.SWE)


#------#
# PLOT #
#------#

# Graphic parameters
# Colours
col.eo  <- "#ffbb44" 
col.m   <- "#859b6c" 
col.a.h <- "#004f63" 

# Labels
colors <- c('Modal age at death' = col.m, 'Threshold age' = col.a.h, 'Life expectancy at birth' = col.eo)    # Colours
lntp   <- c('Modal age at death' = 3,     'Threshold age' = 1,       'Life expectancy at birth' = 4)         # Linetypes
sz     <- c('Modal age at death' = 1.2,   'Threshold age' = 1.2,     'Life expectancy at birth' = 1.2)       # Linewidths

ggplot() +
  geom_line(data = df.eo, aes(x = Year, y = mode, colour = 'Modal age at death', lty = 'Modal age at death', linewidth = 'Modal age at death'))+
  geom_line(data = df.a.h, aes(x = Year, y = a.h, colour = 'Threshold age', lty = 'Threshold age', linewidth = 'Threshold age'))+
  geom_line(data = df.eo, aes(x = Year, y = ex, colour = 'Life expectancy at birth', lty = 'Life expectancy at birth', linewidth = 'Life expectancy at birth'))+
  ylab('Years (age)')+
  xlab('Year (time)')+
  custom_theme()+
  theme(legend.position = c(.75,.14),
        legend.text = element_text(size = 18),
        legend.key.width = unit(2.2, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.spacing = unit(-0.6,"cm")) +
  scale_linetype_manual(values = lntp, breaks = c('Modal age at death', 'Threshold age', 'Approximation', 'Life expectancy at birth'))+
  scale_colour_manual(values = colors)+
  scale_linewidth_manual(values = sz)+
  scale_x_continuous(limits = c(1900, 2020),
                     breaks = seq(1900, 2020, 20))+
  labs(colour = '',
       linetype = '',
       linewidth = '') +
  guides(colour = 'none',
         linewidth = 'none',
         linetype = guide_legend(override.aes = list(colour = colors, linewidth=sz),
                                 label.position = 'left',
                                 ncol = 1))

# Remove everything but functions
rm(list = setdiff(ls(), lsf.str()))


#------------------------------------------------------------------------------#
#          FIGURE 4: THRESHOLD AGE, MODE AND LIFE EXPECTANCY      
#     Females of England and Wales, Italy, Japan and USA, 1900-2021 
#------------------------------------------------------------------------------#

#-------------#
# LIFE TABLES # 
#-------------#

# HUMAN MORTALITY DATABASE https://www.mortality.org

# English and Walish Females from 1900
load("Data/EnglandWales_female_life_tables.RData")  
# Italian Females from 1900
load("Data/Italy_female_life_tables.RData")
# Japanese Females from 1900
load("Data/Japan_female_life_tables.RData")   
# US Females from 1900
load("Data/USA_female_life_tables.RData")             

# Combine the datasets
df.LT <- LT.ENW |> 
  bind_rows(LT.ITA,
            LT.JPN,
            LT.USA)

# Modify the dataset in order to be able to apply the function 'threshold.decomp.function'
df.LT$Sex1 <- df.LT$Sex
df.LT$PopName2 <- df.LT$PopName
rm(LT.ENW, LT.ITA, LT.JPN, LT.USA)


#----------------#
# THRESHOLD AGES #
#----------------#

# Compute the threshold age
Threshold.ages <- df.LT[, c(a.h = threshold.decomp.function(.SD, func2 = h.frommx)), by = list(Sex, PopName)]

# Modify structure of the dataset to add the year
names(Threshold.ages) <- c('Sex', 'PopName','a.h')
year.label <- df.LT[, list(Year = unique(Year)), by = list(Sex, PopName)]

# Dataset with the threshold ages
df.a.h <- cbind(Threshold.ages, year.label[,3])
rm(Threshold.ages)


#---------------------#
# GOMPERTZ PARAMETERS #
#---------------------#

# Compute the Gompertz parameters alpha and beta with 
#  three different age ranges: 30:90, 40:90 and 50:90

# Age range 30 to 90
alpha.beta.30 <- df.LT[, param.gomp(.SD, Age = 30:90), by = PopName]
alpha.beta.30 <- cbind(alpha.beta.30, year.label[,3])

# Age range 40 to 90
alpha.beta.40 <- df.LT[, param.gomp(.SD, Age = 40:90), by = PopName]
alpha.beta.40 <- cbind(alpha.beta.40, year.label[,3])

# Age range 50 to 90
alpha.beta.50 <- df.LT[, param.gomp(.SD, Age = 50:90), by = PopName]
alpha.beta.50 <- cbind(alpha.beta.50, year.label[,3])

# Construct a unique dataset with alpha, beta 
# and the Gompertz modal age at death as M = [log(beta) - log(alpha)] / beta
# obtained with  the three different age ranges
df.m.g <- alpha.beta.30 |>
  filter(Year <= max(df.a.h$Year) & Year >= min(df.a.h$Year)) |> 
  left_join(alpha.beta.40, by = c('Year', 'PopName'), suffix = c('.30', '.40')) |>
  left_join(alpha.beta.50, by = c('Year', 'PopName')) |>
  rename(Alpha.50 = Alpha,
         Beta.50  = Beta) |> 
  relocate(Year, .after=PopName) |> 
  group_by(PopName, Year) |> 
  mutate(mode.30 = gompertz.modal.age(alpha = Alpha.30, beta = Beta.30),      # Gompertz mode obtained with parameters estimated with age range 30:90
         mode.40 = gompertz.modal.age(alpha = Alpha.40, beta = Beta.40),      # Gompertz mode obtained with parameters estimated with age range 40:90
         mode.50 = gompertz.modal.age(alpha = Alpha.50, beta = Beta.50)) |>   # Gompertz mode obtained with parameters estimated with age range 50:90
  ungroup() |> 
  data.table()

rm(alpha.beta.30, alpha.beta.40, alpha.beta.50, year.label)


#-------#
# PLOTS #
#-------#

# Countries
countries <- unique(df.LT$PopName)
countries.names <- c("England and Wales", "Italy", "Japan", "United States")

# Empty lists to fill up with the plots of the countries with the three different Gompertz starting ages
all.plots  <- list()
plot.30    <- list()
plot.40    <- list()
plot.50    <- list()

for(i in 1:length(countries)){
  
  df.a     <- df.a.h[PopName == countries[i], ]                                                           # Threshold Age
  df.eo    <- df.LT[PopName == countries[i] & Age == 0 & Year<=max(df.a$Year) & Year>=min(df.a$Year)]     # Life expectancy at birth
  df.m     <- df.m.g[PopName == countries[i], ]                                                           # Modal age at death
  
  gamma <- - digamma(1) # Euler-Mascheroni constant
  
  # Approximation of threshold age with Eq. (2a)
  form.30 <- data.frame(Approx = approx.ah(df.m$mode.30, df.m$Beta.30), Year = unique(df.a$Year))
  form.40 <- data.frame(Approx = approx.ah(df.m$mode.40, df.m$Beta.40), Year = unique(df.a$Year))
  form.50 <- data.frame(Approx = approx.ah(df.m$mode.50, df.m$Beta.50), Year = unique(df.a$Year))
  
  # Graphic parameters
  # Colours
  col.eo  <- "#ffbb44"
  col.m   <- "#859b6c"
  col.a.h <- "#004f63"
  col.app <- "#ce4441" 
  
  # Labels
  colors <- c('Threshold age' = col.a.h, 'Approximation' = col.app, 'Life expectancy at birth' = col.eo)     # Colours
  lntp   <- c('Threshold age' = 1,       'Approximation' = 2,       'Life expectancy at birth' = 4)          # Linetypes
  sz     <- c('Threshold age' = 1.2,     'Approximation' = 1.2,     'Life expectancy at birth' = 1.2)        # Linewidths
  
  # Gompertz age range 30:90
  df30 <- ggplot() +
    ggtitle('Gompertz starting age: 30')+
    geom_line(data = df.a, aes(x = Year, y = a.h, colour ='Threshold age', lty ='Threshold age', linewidth ='Threshold age'))+
    geom_line(data = form.30, aes(x = Year, y = Approx, colour = 'Approximation', lty ='Approximation', linewidth ='Approximation'))+
    geom_line(data = df.eo, aes(x = Year, y = ex, colour = 'Life expectancy at birth', lty = 'Life expectancy at birth', linewidth = 'Life expectancy at birth'))+
    custom_theme() +
    scale_y_continuous(expression("Years"))+
    scale_x_continuous(expression(" "),
                       limits = c(NA, 2020),
                       breaks = seq(1900, 2020, 20))+
    scale_linetype_manual(breaks = c('Threshold age', 'Approximation', 'Life expectancy at birth'), values = lntp)+
    scale_colour_manual(values = colors)+
    scale_linewidth_manual(values = sz)+
    labs(colour = '',
         linetype = '')+
    guides(colour  =  'none',
           linewidth = 'none',
           linetype = guide_legend(override.aes = list(colour = colors,
                                                       linewidth = c(1,1,1))))
  # Gompertz age range 40:90
  df40 <- ggplot() +
    ggtitle('Gompertz starting age: 40')+
    geom_line(data = df.a, aes(x = Year, y = a.h, colour = 'Threshold age', lty = 'Threshold age', linewidth = 'Threshold age'))+
    geom_line(data = form.40, aes(x = Year, y = Approx, colour = 'Approximation', lty = 'Approximation', linewidth = 'Approximation'))+
    geom_line(data = df.eo, aes(x = Year, y = ex, colour = 'Life expectancy at birth', lty = 'Life expectancy at birth', linewidth = 'Life expectancy at birth'))+
    custom_theme()+
    scale_y_continuous(expression("Years"))+
    scale_x_continuous(expression(" "),
                       limits = c(NA, 2020),
                       breaks = seq(1900, 2020, 20))+
    scale_linetype_manual(breaks = c('Threshold age', 'Approximation', 'Life expectancy at birth'), values = lntp)+
    scale_colour_manual(values = colors)+
    scale_linewidth_manual(values = sz)+
    labs(colour = '',
         linetype = '')+
    guides(colour = 'none',
           linewidth = 'none',
           linetype = guide_legend(override.aes = list(colour = colors,
                                                       linewidth = c(1,1,1))))
  # Gompertz age range 50:90
  df50 <- ggplot() +
    ggtitle('Gompertz starting age: 50')+
    geom_line(data = df.a, aes(x = Year, y = a.h, colour = 'Threshold age', lty = 'Threshold age', linewidth = 'Threshold age'))+
    geom_line(data = form.50, aes(x = Year, y = Approx, colour = 'Approximation', lty ='Approximation', linewidth = 'Approximation'))+
    geom_line(data = df.eo, aes(x = Year, y = ex, colour = 'Life expectancy at birth', lty = 'Life expectancy at birth', linewidth = 'Life expectancy at birth'))+
    custom_theme() +
    scale_y_continuous(expression("Years"))+
    scale_x_continuous(expression(" "),
                       limits = c(NA, 2020),
                       breaks = seq(1900, 2020, 20))+
    scale_linetype_manual(breaks = c('Threshold age', 'Approximation', 'Life expectancy at birth'), values = lntp)+
    scale_colour_manual(values = colors)+
    scale_linewidth_manual(values = sz)+
    labs(colour = '',
         linetype = '')+
    guides(colour = 'none',
           linewidth = 'none',
           linetype = guide_legend(override.aes = list(colour = colors,
                                                       linewidth = c(1,1,1))))
  plot.30[[i]] <- df30
  plot.40[[i]] <- df40
  plot.50[[i]] <- df50
  
  all.plots[[i]] <- (df30 + df40 + df50) +
    plot_annotation(title = countries.names[i],
                    theme = theme(plot.title = element_text(size = 25))) +
    plot_layout(guides = 'collect') &
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 22),
          legend.position = 'bottom',
          legend.key.width = unit(2, 'cm'),
          legend.key.height = unit(1, 'cm'))
  
}

# Final Plots
all.plots[[1]] # England and Wales
all.plots[[2]] # Italy
all.plots[[3]] # Japan
all.plots[[4]] # USA

# Remove everything but functions
rm(list = setdiff(ls(), lsf.str()))


#------------------------------------------------------------------------------#
#         FIGURE 5: CENTRAL DISTRIBUTION OF \hat{a}^H - a^H     
#------------------------------------------------------------------------------#

# Load the data
load("Data/Density_plots.RData")
str(df.density)
# y.approx   =   \hat{a}^H - a^H                         # Figure 5
# y.obs      =   M_obs - (gamma/beta) - a^H              # Figure A-1

#---------#
# FEMALES #
#---------#

df.dens <- df.density |> filter(Sex == 'Female')

# Extract the extremes range of the 95% central distribution of a^H - \hat{a}^H
# 1900-1950
xlims.a <- quantile(df.dens[Year<1951, y.approx], probs = c(.025,.975)) 
# 1951-2000
xlims.b <- quantile(df.dens[Year>1950 & Year<2001, y.approx], probs = c(.025,.975)) # we plot 95% of the data 1951-2000

# Graphic parameters
coll <- c("#ffbb44", "#eb7926", "#62929a", "#859b6c")
brks.y <- seq(0, 1.5, .5)

# DENSITIES 1900-1950
density.plot.40a <- df.dens  |>  
  filter(Year<1951, between(y.approx, xlims.a[1], xlims.a[2])) |> 
  ggplot(aes(y.approx)) +
  geom_density(alpha = 0.6, fill = coll[1], linewidth = .8, colour = coll[1]) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, linetype='dashed', linewidth = .5) +
  custom_theme() +
  scale_x_continuous(limits = c(-2,4), breaks = seq(-2,4,2), labels = rep('', 4)) +
  scale_y_continuous(limits = c(0,1.8)) +
  labs(x='',
       y='')

# DENSITIES 1951-2000
density.plot.40b <- df.dens |> 
  filter(Year>1950 & Year<2001, between(y.approx, xlims.b[1], xlims.b[2])) |> 
  ggplot(aes(y.approx)) +
  geom_density(alpha = 0.6, fill = coll[2], linewidth = .8, colour = coll[2]) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 'dashed', linewidth = .5) +
  custom_theme() +
  scale_x_continuous(limits = c(-2,4), breaks = seq(-2,4,2), labels = rep('', 4)) +
  scale_y_continuous(limits = c(0,1.8), breaks = brks.y, labels = rep('', length(brks.y))) +
  labs(x='',
       y='')


#-------#
# MALES #
#-------#

df.dens.m <- df.density |> filter(Sex == 'Male')

# Extract the extremes range of the 95% central distribution of a^H - \hat{a}^H
# 1900-1950
xlims.m.a <- quantile(df.dens.m[Year<1951, y.approx], probs = c(.025,.975))
# 1951-2000
xlims.m.b <- quantile(df.dens.m[Year>1950 & Year<2001, y.approx], probs = c(.025,.975))

# DENSITIES 1900-1950
density.plot.50a.m <- df.dens.m |> 
  filter(Year<1951, between(y.approx, xlims.m.a[1], xlims.m.a[2])) |> 
  ggplot(aes(y.approx)) +
  geom_density(alpha = 0.6, fill = coll[3], linewidth = .8, colour = coll[3]) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, linetype ='dashed', linewidth =.5) +
  custom_theme()+
  scale_x_continuous(limits = c(-2,4), breaks = seq(-2,4,2)) +
  scale_y_continuous(limits = c(0,1.8)) +
  labs(x='',
       y='') 

# DENSITIES 1951-2000
density.plot.50b.m <- df.dens.m |> 
  filter(Year>1950 & Year<2001, between(y.approx, xlims.m.b[1], xlims.m.b[2])) |> 
  ggplot(aes(y.approx)) +
  geom_density(alpha = 0.6, fill = coll[4], linewidth = .8, colour = coll[4]) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 'dashed', linewidth = .5) +
  custom_theme()+
  scale_x_continuous(limits = c(-2,4), breaks = seq(-2,4,2)) +
  scale_y_continuous(limits = c(0,1.8), breaks = brks.y, labels = rep('', length(brks.y))) +
  labs(x='',
       y='') 


# Final Plot
(density.plot.40a + density.plot.40b) /        # Females            
  (density.plot.50a.m + density.plot.50b.m)    # Males


#-------------------#
# IN-TEXT ESTIMATES #                      
#-------------------#

## Means and modes of the densities

# Females mean
mean(df.dens[Year<1951, y.approx])                            # -0.2744423
mean(df.dens[Year>1950 & Year<2001, y.approx])                #  0.2470369

# Males mean
mean(df.dens.m[Year<1951, y.approx])                          #  0.5056341
mean(df.dens.m[Year>1950 & Year<2001, y.approx])              # -0.2244362

# Females mode
estimate_mode(df.dens[Year<1951, y.approx])                   # -0.1087201
estimate_mode(df.dens[Year>1950 & Year<2001, y.approx])       #  0.3383297

# Males mode
estimate_mode(df.dens.m[Year<1951, y.approx])                 #  0.2605978
estimate_mode(df.dens.m[Year>1950 & Year<2001, y.approx])     #  0.04200034

# Remove everything but functions
rm(list = setdiff(ls(), lsf.str()))


#------------------------------------------------------------------------------#
#   FIGURE 6: Age-at-death distribution of Swedish females in 1900 and 2000     
#------------------------------------------------------------------------------#

# Swedish female life tables from 1900 to 2021 - Human Mortality Database
load("Data/Sweden_female_life_tables.RData") 

# NOTE: The modal age at death (column 'mode') was calculated using the
#  'MortalitySmooth' R package. This package is no longer available in CRAN,
#  but can be downloaded from Tim Riffe's GitHub using the 'devtools' R package.
#  https://github.com/timriffe/MortalitySmooth

# Select years of interest and modify the dataset to compute threshold ages
df.LT <- LT.SWE |> 
  filter(Year %in% c(1900, 2000)) |> 
  mutate(Sex1 = Sex,
         PopName2 = PopName)
rm(LT.SWE)

# Compute the threshold age
a.h.SWE <- df.LT[, c(a.h = threshold.decomp.function(.SD, func2 = h.frommx)), by = list(Sex, PopName)]
names(a.h.SWE) <- c('Sex', 'PopName', 'a.h')
a.h.SWE <- a.h.SWE |> mutate(Year = c(1900, 2000))


#-------#
# PLOTS #
#-------#

# Graphic parameters
# Colours
col.eo  <- "#ffbb44"
col.m   <- "#859b6c"
col.a.h <- "#004f63"

# Labels
colors  <- c('Threshold age' = col.a.h, 'Modal age at death' = col.m, 'Life expectancy at birth' = col.eo)    # Colours
lntp    <- c('Threshold age' = 1,       'Modal age at death' = 3,     'Life expectancy at birth' = 4)         # Linetypes

# Plot
df.LT |> 
  ggplot(aes(x = Age, y = dx/10^5)) + 
  geom_line(linewidth = 1.2, colour = 'darkgrey') +
  facet_wrap(.~Year) +
  geom_vline(data = filter(a.h.SWE, Year == 1900), aes(xintercept = a.h, colour = 'Threshold age', lty = 'Threshold age'), linewidth = 1.2) +
  geom_vline(data = filter(a.h.SWE, Year == 2000), aes(xintercept = a.h, colour = 'Threshold age', lty = 'Threshold age'), linewidth = 1.2) +
  geom_vline(data = filter(df.LT, Year == 1900), aes(xintercept = ex[Age == 0], colour = 'Life expectancy at birth', lty = 'Life expectancy at birth'), linewidth = 1.2) +
  geom_vline(data = filter(df.LT, Year == 2000), aes(xintercept = ex[Age == 0], colour = 'Life expectancy at birth', lty = 'Life expectancy at birth'), linewidth = 1.2) +
  geom_vline(data = filter(df.LT, Year == 1900), aes(xintercept = mode[Age == 0], colour = 'Modal age at death', lty = 'Modal age at death'), linewidth = 1.2) +
  geom_vline(data = filter(df.LT, Year == 2000), aes(xintercept = mode[Age == 0], colour = 'Modal age at death', lty = 'Modal age at death'), linewidth = 1.2) +
  custom_theme() +
  scale_y_continuous(breaks = seq(0, 0.075, 0.025),
                     labels = c(0, 0.025, 0.050, 0.075))+
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 20),
        legend.key.height = unit(1.2, 'cm'),
        strip.text.x = element_text(hjust = 0, size = 22)) +
  labs(x = 'Years (age)',
       colour = '',
       linetype = '') +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = lntp, breaks = c('Threshold age', 'Modal age at death', 'Life expectancy at birth')) +
  guides(colour = 'none',
         linewidth = 'none',
         linetype = guide_legend(override.aes = list(colour = colors, linewidth = 1.2)))

# Remove everything but functions
rm(list = setdiff(ls(), lsf.str()))


#------------------------------------------------------------------------------#
#      FIGURE A-1: CENTRAL DISTRIBUTION OF M_obs - gamma/beta - a^H      
#------------------------------------------------------------------------------#

# Load the data
load("Data/Density_plots.RData")
str(df.density)
# y.approx   =   \hat{a}^H - a^H                         # Figure 5
# y.obs      =   M_obs - gamma/beta - a^H                # Figure A-1

#---------#
# FEMALES #
#---------#

df.dens <- df.density |> filter(Sex == 'Female')

# Extract the range of the 95% central distribution of M_obs - gamma/beta - a^H 
# 1900-1950
xlims.a <- quantile(df.dens[Year<1951, y.obs], probs = c(.025,.975))
# 1951-2000
xlims.b <- quantile(df.dens[Year>1950 & Year<2001, y.obs], probs = c(.025,.975))

# Graphic parameters
coll <- c("#ffbb44", "#eb7926", "#62929a", "#859b6c")
brks.y <- seq(0, 1.5, .5)

# DENSITIES 1900-1950
density.plot.40a <- df.dens  |>  
  filter(Year<1951, between(y.obs, xlims.a[1], xlims.a[2])) |> 
  ggplot(aes(y.obs)) +
  geom_density(alpha = 0.6, fill = coll[1], linewidth = .8, colour = coll[1]) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 'dashed', linewidth = .5) +
  custom_theme() +
  scale_x_continuous(limits = c(-2,4), breaks = seq(-2,4,2), labels = rep('', 4)) +
  scale_y_continuous(limits = c(0,1.8)) +
  labs(x='',
       y='')

# DENSITIES 1951-2000
density.plot.40b <- df.dens |> 
  filter(Year>1950 & Year<2001, between(y.obs, xlims.b[1], xlims.b[2])) |> 
  ggplot(aes(y.obs)) +
  geom_density(alpha = 0.6, fill = coll[2], linewidth = .8, colour = coll[2]) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 'dashed', linewidth = .5) +
  custom_theme() +
  scale_x_continuous(limits = c(-2,4), breaks = seq(-2,4,2), labels = rep('', 4)) +
  scale_y_continuous(limits = c(0,1.8), breaks = brks.y, labels = rep('', length(brks.y))) +
  labs(x='',
       y='')


#-------#
# MALES #
#-------#

df.dens.m <- df.density |> filter(Sex == 'Male')

# Extract the range of the 95% central distribution of a^H - \hat{a}^H
# 1900-1950
xlims.m.a <- quantile(df.dens.m[Year<1951, y.obs], probs = c(.025,.975))
# 1951-2000
xlims.m.b <- quantile(df.dens.m[Year>1950 & Year<2001, y.obs], probs = c(.025,.975))

# DENSITIES 1900-1950
density.plot.50a.m <- df.dens.m |> 
  filter(Year<1951, between(y.obs, xlims.m.a[1], xlims.m.a[2])) |> 
  ggplot(aes(y.obs)) +
  geom_density(alpha = 0.6, fill = coll[3], linewidth = .8, colour = coll[3]) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 'dashed', linewidth = .5) +
  custom_theme()+
  scale_x_continuous(limits = c(-2,4), breaks = seq(-2,4,2)) +
  scale_y_continuous(limits = c(0,1.8)) +
  labs(x='',
       y='') 

# DENSITIES 1951-2000
density.plot.50b.m <- df.dens.m |> 
  filter(Year>1950 & Year<2001, between(y.obs, xlims.m.b[1], xlims.m.b[2])) |> 
  ggplot(aes(y.obs)) +
  geom_density(alpha = 0.6, fill = coll[4], linewidth = .8, colour = coll[4]) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 'dashed', linewidth = .5) +
  custom_theme()+
  scale_x_continuous(limits = c(-2,4), breaks = seq(-2,4,2)) +
  scale_y_continuous(limits = c(0,1.8), breaks = brks.y, labels = rep('', length(brks.y))) +
  labs(x='',
       y='') 

# Final Plot
(density.plot.40a + density.plot.40b) /        # Females 
  (density.plot.50a.m + density.plot.50b.m)    # Males

# Remove everything but functions
rm(list = setdiff(ls(), lsf.str()))


#-------------------------------------------------------------------------------#
# FIGURE A-2: THRESHOLD AGE, MODE AND LIFE EXPECTANCY SPANISH FEMALES 1908-2021                     
#-------------------------------------------------------------------------------#

# Spanish female life tables from 1908 to 2020 - Human Mortality Database
load("Data/Spain_female_life_tables.RData")

# NOTE: The modal age at death (column 'mode') was calculated using the
#  'MortalitySmooth' R package. This package is no longer available in CRAN,
#  but can be downloaded from Tim Riffe's GitHub using the 'devtools' R package.
#  https://github.com/timriffe/MortalitySmooth

# Females
df.LT <- LT.ESP
df.LT$Sex1 <- df.LT$Sex
df.LT$PopName2 <- df.LT$PopName

# Compute the threshold age
Threshold.ages <- df.LT[, c(a.h = threshold.decomp.function(.SD, func2 = h.frommx)), by = list(Sex, PopName)]

# Modify structure of the dataset to add the year
names(Threshold.ages) <- c('Sex', 'PopName','a.h')
year.label <- df.LT[, list(Year = unique(Year)), by = list(Sex,PopName)]
df.a.h <- cbind(Threshold.ages, year.label[,3])

# Compute the Gompertz parameters to use in the approximation
alpha.beta.40 <- df.LT[, param.gomp(.SD, Age = 40:90), by = PopName]
alpha.beta.40 <- cbind(alpha.beta.40, year.label[,3])

# Compute the approximation with Eq.(2a)
form.40 <- alpha.beta.40 |>
  filter(Year <= max(df.a.h$Year) & Year >= min(df.a.h$Year)) |> 
  mutate(mode = gompertz.modal.age(alpha = Alpha, beta = Beta),
         Approx = mode + digamma(1)/Beta) |> 
  data.table()
rm(year.label, df.LT, Threshold.ages, alpha.beta.40)

# Extract life expectancy at birth and modal age at death
df.eo <- LT.ESP[Age == 0, ]
rm(LT.ESP)


#------#
# PLOT #
#------#

# Graphic parameters
# Colours
col.eo  <- "#ffbb44"
col.m   <- "#859b6c"
col.a.h <- "#004f63"
col.app <- "#ce4441"

# Labels
colors  <- c('Modal age at death' = col.m, 'Threshold age' = col.a.h, 'Approximation' = col.app,  'Life expectancy at birth' = col.eo)    # Colours
lntp    <- c('Modal age at death' = 3,     'Threshold age' = 1,       'Approximation' = 2,        'Life expectancy at birth' = 4)         # Linetypes
sz      <- c('Modal age at death' = 1.2,   'Threshold age' = 1.2,     'Approximation' = 1.2,      'Life expectancy at birth' = 1.2)       # Linewidths

ggplot() +
  geom_line(data = df.eo, aes(x = Year, y = mode, colour = 'Modal age at death', lty = 'Modal age at death', linewidth = 'Modal age at death'))+
  geom_line(data = df.a.h, aes(x = Year, y = a.h, colour = 'Threshold age', lty = 'Threshold age', linewidth='Threshold age'))+
  geom_line(data = df.eo, aes(x = Year, y = ex, colour = 'Life expectancy at birth', lty = 'Life expectancy at birth', linewidth = 'Life expectancy at birth'))+
  geom_line(data = form.40, aes(x = Year, y = Approx, colour = 'Approximation', lty = 'Approximation', linewidth = 'Approximation'))+
  ylab('Years (age)')+
  xlab('Year (time)')+
  custom_theme()+
  theme(legend.position = c(.75,.14),
        legend.text = element_text(size = 18),
        legend.key.width = unit(2.2, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.margin = unit(-0.6,"cm")) +
  scale_linetype_manual(values = lntp, breaks = c('Modal age at death', 'Threshold age', 'Approximation', 'Life expectancy at birth'))+
  scale_colour_manual(values = colors)+
  scale_linewidth_manual(values = sz)+
  scale_x_continuous(limits = c(1900, 2020),
                     breaks = seq(1900, 2020, 20))+
  labs(colour = '',
       linetype = '',
       linewidth = '') +
  guides(colour = 'none',
         linewidth = 'none',
         linetype = guide_legend(override.aes = list(colour = colors, linewidth = sz),
                                 label.position = 'left',
                                 ncol = 1))
