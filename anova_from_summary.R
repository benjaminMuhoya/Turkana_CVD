
setwd("/Users/bm0211/RegEx/Turkana_CVD")
## Function to Run ANOVA on several groups ##################
FNONCT <- function(x,df1,df2,prob, interval=c(0,10000), my.tol=0.000001){
		temp <- function(ncp) pf(x,df1,df2,ncp) - prob
		return(uniroot(temp, interval, tol = my.tol)$root)
}
######################################
power.f <- function(u, n, delta, sig.level=.05){
  fc <- qf(p=sig.level,df1=u,df2=(u+1)*(n-1),lower.tail=FALSE)
  lamda <- (delta^2)*(n*(u+1))
  v <- (u+1)*(n-1)

  z1b <- (sqrt(2*(u+lamda)-((u+2*lamda)/(u+lamda)))-
  sqrt((2*v-1)*((u*fc)/v)))/
  sqrt(((u*fc)/v)+((u+2*lamda)/(u+lamda)))
  output <- pnorm(z1b)
  return(output)
}
library(gtools)
###################################################################
ind.oneway.second <- function(m, sd, n, unbiased=TRUE, contr=NULL, sig.level=.05, digits=3){
##if orthogonal
if(length(n)==1){
     n <- rep(n, length(m))
}
#if biased standard deviation
if(unbiased==FALSE){
         sd <- ssd2sd(n,sd)
}
##(a) anova table
k  <- length(m)           #number of groups
Xg  <- sum(n*m)/sum(n)

dfb <- k - 1              #degree of freedom
dfw   <- sum(n) - k       #degree of freedom

MSb <- sum(n * (m - Xg)^2)/(k-1)  #MS between
MSw <-  sum((n-1)*sd^2)/dfw       #MS within
SSb <- dfb * MSb
SSw <- dfw * MSw
SSt <- SSb + SSw

f.value <- MSb/MSw                #f value

anova.table  <- data.frame(matrix(NA,ncol=4, nrow=3))
rownames(anova.table) <- c("Between (A)", "Within", "Total")
colnames(anova.table) <- c("SS", "df", "MS", "F")
anova.table$SS <- c(SSb, SSw,SSt)
anova.table$df <- c(dfb, dfw, dfb+dfw)
anova.table$MS <- c(MSb,MSw,NA)
anova.table$F  <- c(f.value, NA,NA)
class(anova.table) <- c("anova", "data.frame")
anova.table <- round(anova.table, digits)
##(b) omnibus effect size eta and omega squared
#eta square
etasq <- SSb / SSt
delta.lower <- delta.upper <- numeric(length(etasq))
delta.lower <- try(FNONCT(f.value, dfb, dfw, prob=1-sig.level/2), silent=TRUE)
delta.upper <- try(FNONCT(f.value, dfb, dfw, prob=sig.level/2), silent=TRUE)
if(is.character(delta.lower)){
  delta.lower <- 0
}

etasq.lower <- delta.lower / (delta.lower + dfb + dfw + 1)
etasq.upper <- delta.upper / (delta.upper + dfb + dfw + 1)

#omega square
omegasq <- (SSb - dfb * MSw)/(SSt + MSw)
sosb_L  <- SSt * etasq.lower
msw_L   <- (SSt - sosb_L)/dfw
omegasq.lower <- (sosb_L - (dfb*msw_L))/(SSt+msw_L)

sosb_U  <- SSt * etasq.upper
msw_U   <- (SSt - sosb_U)/dfw
omegasq.upper <- (sosb_U - (dfb*msw_U))/(SSt+msw_U)

omnibus.es <- round(c(etasq=etasq, etasq.lower=etasq.lower, etasq.upper=etasq.upper),
                      digits)
##(c) raw contrasts
temp  <- combinations(k,2)
cont1 <- matrix(0, nrow=nrow(temp),ncol=k)
cont1.lab <- rep(0,nrow(temp))
#in case did not specify contrasts
for(i in 1:nrow(temp)){
     cont1[i, temp[i,1]] <- 1
     cont1[i, temp[i,2]] <- -1
     cont1.lab[i] <- paste(temp[i,1],"-",temp[i,2], sep="")
     rownames(cont1) <- cont1.lab
}
#in case specify contrasts
if(!is.null(contr)){
      if(is.vector(contr)){
                cont1 <- t(as.matrix(contr))
      }else{
                cont1 <- contr
      }
}
#F test for contrasts
psi <-  colSums(t(cont1)  * as.vector(m))                #raw contrasts
SSpsi <- (psi^2)/colSums(t(cont1^2) / as.vector(n))

nmat <- matrix(n, nrow=nrow(cont1), ncol=length(n), byrow=TRUE)
psi.std <- sqrt(MSw * rowSums(cont1 ^ 2/nmat))

psi.lower <- psi + psi.std * qt(sig.level/2, dfw)
psi.upper <- psi + psi.std * qt(sig.level/2, dfw, lower.tail=FALSE)

raw.contrasts <- round(data.frame(mean.diff=psi, lower=psi.lower, upper=psi.upper, std=psi.std), digits)
rownames(raw.contrasts) <- rownames(cont1)
##(d) standardized contrasts
gpsi <- psi/sqrt(MSw)       #effect size
gpsi.std <- sqrt(rowSums(cont1 ^ 2/nmat))
gpsi.lower <- gpsi + gpsi.std * qt(sig.level/2, dfw)
gpsi.upper <- gpsi + gpsi.std * qt(sig.level/2, dfw, lower.tail=FALSE)
standardized.contrasts <- round(data.frame(es=gpsi, lower=gpsi.lower, upper=gpsi.upper, std=gpsi.std), digits)
rownames(standardized.contrasts) <- rownames(cont1)
##(e) statistical power
c.delta <- c(.10, .25, .4)
criterion.power <- round(power.f(sig.level=sig.level, u=dfb, n=sum(n)/k,delta=c.delta), digits)
names(criterion.power) <- c("small", "medium", "large")
##(e) output
output <- list(anova.table=anova.table, omnibus.es=omnibus.es, raw.contrasts=raw.contrasts, standardized.contrasts = standardized.contrasts, power=criterion.power)
return(output)
}
###############################################################################################
##Fake Data to test the Function
mean <- c(90,85,92,100,102,106)
sd <- c(9.035613,11.479667,9.760268,7.662572,9.830258,9.111457)
no_of_participants <- c(9,9,9,9,9,9)
ind.oneway.second (mean, sd, no_of_participants) ##This this the function and compares Max of 7 groups at once

## HDL SECTION #############################################################
## My Data seven regions chosen, listed in order as the occur in this list 
## Pastoralist, Ghana, Inuits, India, Japan, Spain, Hispanic in America
m <- c(60.97, 46.3906, 60.5041, 39, 52.8344, 52.7, 48.73)
sd <- c(15.1,14.6886, 17.9685, 12, 14.345, 13.94, 13.17)
n <- c(800, 3317, 5405, 2042, 18077, 51462, 12730)
ind.oneway.second (m, sd, n)

## Second : 
## Pastoralist, Mexico 2016, Aboriginals, TSI, Samoa, Tibetian, White
m <- c(60.97, 43.9, 32.0183, 38.67, 52.5457, 56.4582, 54.138)
sd <- c(15.1,12.2, 2.1965, 4.1251, 11.7545, 11.601, 15.468)
n <- c(800,586, 621, 304, 3451, 587, 46788)
ind.oneway.second (m, sd, n)

##Third: 
## Pas, UK Biobank-Healthy, Australia ANHF study, NHANES-III, Tsimane, Singapore, Thailand, Nairobi, Turkana_Overall, Shuar 
m <- c(61.46, 57, 50.271, 50.9, 36.8, 50.271, 46.404, 48.7242, 58.33, 51.2616)
sd <- c(15.7,14.3, 15.468, 15.1, 8.9, 11.601, 11.601, 18.1749, 15.5, 15.2328)
n <- c(800, 340470, 8924, 10843, 356, 3223, 4876, 2003, 2300, 315)
ind.oneway.second (m, sd, n)

## LDL Section ###############################################################
##My Data seven regions chosen, listed in order as the occur in this list 
## Turkana Pastoralist, Kenya Nairobi, Ghana, Inuits, Australia ANHF study, NHANES III, White
m <- c(67.36, 117.1701, 101.9239, 129.7211, 135.345, 127, 146.946)
sd <- c(30.5, 44.0838, 36.9395, 39.5453, 38.67, 38.5, 42.537, 34.803, 42.537)
n <- c(800, 2003, 3317, 5405, 8924, 10843, 46788)
ind.oneway.second (m, sd, n)

## Second : 
## Pastoralist, Singapore, Thailand, Japan, TSI, Aboriginals, Tsimane
m <- c(67.37, 135.345, 135.345, 119.3713, 141.5526, 152.1973, 70.6)
sd <- c(30.5, 34.803, 42.537, 31.2704, 15.2073, 12.8488, 21.9)
n <- c(800, 3223, 4876, 18077, 304, 621, 356)
ind.oneway.second (m, sd, n)

##Third: 
## Pas,Turkana Overall, Mexico 2016, Hispanic in America, Tibetian, India
m <- c(67.37, 66.2, 91.9, 128.87, 124.1307, 91.1557)
sd <- c(30.5, 31.2, 25.4, 38.13, 35.5764, 32.2409)
n <- c(800,2300, 586, 12730, 587, 4084)
ind.oneway.second (m, sd, n)

##Fourth: Pas, Spain, Samoa, Shuar
m <- c(67.37,126.2, 149.2207, 96.054)
sd <- c(30.5, 34.33, 40.3461, 32.241)
n <- c(800,51462, 3451, 4084)
ind.oneway.second (m, sd, n)

## Chol Section
##My Data seven regions chosen, listed in order as the occur in this list Turkana Pastoralist,
## Ghana, Inuits, India, Japan, Spain, Hispanic in America
m <- c(149, 167.7339, 220.0125, 158, 197.8987, 211.01, 205.78)
sd <- c(33.19, 43.1175, 46.2037, 40, 34.7322, 40.97, 45.66)
n <- c(800, 3317, 1420, 2042, 18077, 51462, 12730)
ind.oneway.second (m, sd, n)

## Second : 
## Pastoralist, Mexico 2016, Aboriginals, TSI, Samoa, Tibetian, White
m <- c(149, 159.8, 186.5027, 183.4993, 201.084, 197.6037, 224.286)
sd <- c(33.19, 31.9, 24.4277, 14.7794, 36.4395, 39.8301, 46.404)
n <- c(800, 586, 621, 304, 3451, 587, 46788)
ind.oneway.second (m, sd, n)

##Third: 
## Pas,UK Biobank-Healthy, Australia ANHF study, NHANES-III, Tsimane, Singapore, Thailand, Turkana Overall, Shuar
m <- c(149, 212.685, 203.3, 138, 204.951, 208.818, 148.15, 171.8648)
sd <- c(33.19, 42.537, 44, 29.2, 38.67, 46.404, 33.3, 42.0112)
n <- c(800, 8924, 10843, 356, 3223, 4876, 2300, 315)
ind.oneway.second (m, sd, n)
############
## Write the results in the Meta_analysis_data_file) then Graph
##########

library(ggplot2)
my_colors_P <- c("#FAD510", "#008080", "lightpink", "green","#FF2400",
  "cyan", "#90A959", "#9D858D", "#A4243B", "#6495ED", "#5B1A18", "#1B2021"
)
## Graphs ########################################################
result_anova = import("/Users/bm0211/RegEx/Turkana_CVD/Meta_analysis_Data_and_Results.xlsx", sheet="ResultsLDL") ## For LDL
head(result_anova)
names(result_anova)[1] = "Pairwise"
##
result_anova$Ancestry <- factor(result_anova$Ancestry)
BoxPlot_ANOVA_meta <-ggplot(result_anova, aes(x= reorder(Pairwise, es), y = es))+ geom_point(aes(color = factor(Ancestry)), size=8) +scale_color_manual(values = my_colors_P)+ geom_errorbar(aes(ymin=lower, ymax=upper, width=.2)) +
xlab("Comparison") + ylab("Standardized difference with lower and upper CI") + ggtheme + coord_flip() + geom_hline(yintercept=-1) + geom_hline(yintercept=0.05) + labs(title =  "LDL Distribution Differences", color = "Ancestry")
##Save file as PNG
ggsave("ANOVA_LDL.png", BoxPlot_ANOVA_meta, width = 30, height = 50, units = "cm", dpi = 300)

## Graphs ####################################################################################
result_anova = import("/Users/bm0211/RegEx/Turkana_CVD/Meta_analysis_Data_and_Results.xlsx", sheet="ANOVA ResultsHDL") ## For HDL
names(result_anova)[1] = "Pairwise"
##
result_anova$Ancestry <- factor(result_anova$Ancestry)
BoxPlot_ANOVA_meta <-ggplot(result_anova, aes(x= reorder(Pairwise, es), y = es))+ geom_point(aes(color = factor(Ancestry)), size=8) +scale_color_manual(values = my_colors_P)+ geom_errorbar(aes(ymin=lower, ymax=upper, width=.2)) +
  xlab("Comparison") + ylab("Standardized difference with lower and upper CI") + ggtheme + coord_flip() + geom_hline(yintercept=1) + geom_hline(yintercept=0.23) + labs(title =  "HDL Distribution Differences", color = "Ancestry")
##Save file as PNG
ggsave("ANOVA_HDL.png", BoxPlot_ANOVA_meta, width = 30, height = 50, units = "cm", dpi = 300)

## Graphs ###################################################################################
result_anova = import("/Users/bm0211/RegEx/Turkana_CVD/Meta_analysis_Data_and_Results.xlsx", sheet="CholResults") ## For Chol
names(result_anova)[1] = "Pairwise"
##
result_anova$Ancestry <- factor(result_anova$Ancestry)
BoxPlot_ANOVA_meta <-ggplot(result_anova, aes(x= reorder(Pairwise, es), y = es))+ geom_point(aes(color = factor(Ancestry)), size=8) +scale_color_manual(values = my_colors_P)+ geom_errorbar(aes(ymin=lower, ymax=upper, width=.2)) +
  xlab("Comparison") + ylab("Standardized difference with lower and upper CI") + ggtheme + coord_flip() + geom_hline(yintercept=-1) + geom_hline(yintercept=0.03) + labs(title =  "Cholesterol Distribution Differences", color = "Ancestry")
##Save file as PNG
ggsave("ANOVA_chol.png", BoxPlot_ANOVA_meta, width = 30, height = 50, units = "cm", dpi = 300)

