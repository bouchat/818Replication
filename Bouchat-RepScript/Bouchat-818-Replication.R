#Sarah Bouchat
#PS818-Maximum Likelihood Estimation
#Fall 2013
#Replication Project: Malesky & Schuler (2010)


setwd("/Users/sbouchat/Dropbox/WISC/Classes/2013-Fall/MLE/Replication")
library(foreign) 
library(xtable) 
library(MASS)
library(apsrtable)
library(sandwich)
library(lmtest)
library(sm)
library(pscl)
library(ggplot2)
library(boot)
library(RColorBrewer)
library(munsell)
library(stargazer)
library(texreg)
library(car)
library(coefplot)
library(rms)
library(arm)
library(ROCR)
library(separationplot)
library(bootstrap)
library(coda)
library(reshape)
# install.packages("coefplot2",
 # repos="http://www.math.mcmaster.ca/bolker/R",
 # type="source")
library(coefplot2)
library(glmmADMB)
library(lme4)
library(aod)
library(ZIM)

#########################################################################
# REPLICATION 
#########################################################################

###### LOAD DATA ###### 
apsr<-read.dta("apsr.dta")
summary(apsr)

# Exclude missing observation in row 185, where "party" is NA
apsr.n<-apsr[-185,]
attach(apsr.n)


# Function for robust clustered standard errors
clx.HC1<-function(fm,dfcw,cluster){
M <- length(unique(cluster))
N <- length(cluster)
dfc <- (M/(M-1))*((N-1)/(N-fm$rank))
u <- apply(estfun(fm),2,
function(x) tapply(x, cluster, sum))
vcovCL <- dfc*sandwich(fm, meatHC1=crossprod(u)/N)*dfcw
coeftest(fm, vcovCL)
}

###############################
# REPLICATION OF MODELS BY DV
###############################


## DV1: SPEAKNUM_COUNT ##

# H_1: Delegates with safe seats speak less
nb.dv1.1<-glm.nb(speaknum_count ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime, data=apsr.n, control=glm.control(maxit=1000))

nb.dv1.1.se<-clx.HC1(nb.dv1.1, 1, apsr.n$pci_id)


nb.dv1.2<-glm.nb(speaknum_count ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage, data=apsr.n, control=glm.control(maxit=1000))

nb.dv1.2.se<-clx.HC1(nb.dv1.2, 1, apsr.n$pci_id)


# H_4: Transfer recipients are less likely to complain
nb.dv1.3<-glm.nb(speaknum_count ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary, data=apsr.n, control=glm.control(maxit=1000))

nb.dv1.3.se<-clx.HC1(nb.dv1.3, 1, apsr.n$pci_id)


# H_5: Career FE
nb.dv1.4<-glm.nb(speaknum_count ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + careerfx, data=apsr.n, control=glm.control(maxit=1000))

nb.dv1.4.se<-clx.HC1(nb.dv1.4, 1, apsr.n$pci_id)


# Interaction of Fulltime*Centralnominated
nb.dv1.5<-glm.nb(speaknum_count ~ age+ male+ degree+ south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(career2), data=apsr.n, control=glm.control(maxit=1000))

nb.dv1.5.se<-clx.HC1(nb.dv1.5, 1, apsr.n$pci_id)

# Malesky and Schuler estimate the final model for each DV with simulations using "estsimp" from the clarify package in STATA, setting variables to their means

###############################################################################

## DV2: QUESTION COUNT ##

# H_1: Delegates with safe seats speak less
nb.dv2.1<-glm.nb(question_count ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime, data=apsr.n, control=glm.control(maxit=1000))

nb.dv2.1.se<-clx.HC1(nb.dv2.1, 1, apsr.n$pci_id)


nb.dv2.2<-glm.nb(question_count ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage, data=apsr.n, control=glm.control(maxit=1000))

nb.dv2.2.se<-clx.HC1(nb.dv2.2, 1, apsr.n$pci_id)


# H_4: Transfer recipients are less likely to complain
nb.dv2.3<-glm.nb(question_count ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary, data=apsr.n, control=glm.control(maxit=1000))

nb.dv2.3.se<-clx.HC1(nb.dv2.3, 1, apsr.n$pci_id)


# H_5: Career FE
nb.dv2.4<-glm.nb(question_count ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(career2), data=apsr.n)

nb.dv2.4.se<-clx.HC1(nb.dv2.4, 1, apsr.n$pci_id)


# Interaction of Fulltime*Centralnominated
nb.dv2.5<-glm.nb(question_count ~ age+ male+ degree+ south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(career2), data=apsr.n)

# The previous 2 models omit the glm.control in order to produce initial output. The model's variance covariance matrix is not positive definite. In order to achieve the same point estimates and standard errors, the glmmADMB package is also used, but the SEs should still not be relied upon.

nb.dv2.5.se<-clx.HC1(nb.dv2.5, 1, apsr.n$pci_id)

# glmmADMB defaults to a type-2 negative binomial, which counts failures til success. This data is better modeled as a type 1, which counts successes. 

glmm.nb.dv2.4.1<-glmmadmb(question_count ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(career2), data=apsr.n, family="nbinom1")

glmm.nb.dv2.5.1<-glmmadmb(question_count ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated*fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(career2), data=apsr.n, family="nbinom1")

###############################################################################

## DV3: CRITICIZE (crit_per)##

# H_1: Delegates with safe seats speak less
apsr.n.2<-na.omit(data.frame(ln_crit,age, male, south, party, degree, incumbent_dummy, centralnominated, fulltime, percentage, pop_stacked, transfer, percen_asphalted_roads, secondary, as.factor(career2), pci_id, vcp))

lm.dv3.1<-lm(ln_crit ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime, data=apsr.n.2)

lm.dv3.2<-lm(ln_crit ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage, data=apsr.n.2)

# H_4: Transfer recipients are less likely to complain
lm.dv3.3<-lm(ln_crit ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary, data=apsr.n.2)

# H_5: Career FE
lm.dv3.4<-lm(ln_crit ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(vcp), data=apsr.n.2)

# Interaction of Fulltime*Centralnominated
lm.dv3.5<-lm(ln_crit ~ age+ male+ degree+ south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(vcp), data=apsr.n.2)

#Robust SEs
lm.dv3.5$newse<-vcovHC(lm.dv3.5, type="HC1")
DV3<-coeftest(lm.dv3.5,lm.dv3.5$newse)

# Dropping 0 Questions 
ln_crit.n0<-subset(apsr.n, question_count!=0)

lm.dv3.n0<-lm(ln_crit.n0$ln_crit ~ ln_crit.n0$age+ ln_crit.n0$male+ ln_crit.n0$degree+ ln_crit.n0$south + ln_crit.n0$party + ln_crit.n0$incumbent_dummy + (ln_crit.n0$centralnominated*ln_crit.n0$fulltime) + ln_crit.n0$percentage + ln_crit.n0$pop_stacked + ln_crit.n0$transfer + ln_crit.n0$percen_asphalted_roads + ln_crit.n0$secondary, data=ln_crit.n0)

#Robust SEs
lm.dv3.n0$newse<-vcovHC(lm.dv3.n0, type="HC1")
DV3.n0<-coeftest(lm.dv3.n0,lm.dv3.n0$newse)

###############################################################################

## DV4: LOCAL ISSUES ##

# H_1: Delegates with safe seats speak less
apsr.n.3<-na.omit(data.frame(ln_local, age, male, south, party, degree, incumbent_dummy, centralnominated, fulltime, percentage, pop_stacked, transfer, percen_asphalted_roads, secondary, as.factor(career2), pci_id, vcp))

lm.dv4.1<-lm(ln_local ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary, data=apsr.n.3)

# H_5: Career FE
lm.dv4.2<-lm(ln_local ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(vcp), data=apsr.n.3)

# Interaction of Fulltime*Centralnominated
lm.dv4.3<-lm(ln_local ~ age+ male+ degree+ south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + vcp, data=apsr.n.3)

#Robust SEs
lm.dv4.3$newse<-vcovHC(lm.dv4.3, type="HC1")
DV4<-coeftest(lm.dv4.3,lm.dv4.3$newse)

# Dropping 0 Questions 
ln_local.n0<-subset(apsr.n, question_count!=0)

lm.dv4.n0<-lm(ln_local.n0$ln_local ~ ln_local.n0$age+ ln_local.n0$male+ ln_local.n0$degree+ ln_local.n0$south + ln_local.n0$party + ln_local.n0$incumbent_dummy + (ln_local.n0$centralnominated*ln_local.n0$fulltime) + ln_local.n0$percentage + ln_local.n0$pop_stacked + ln_local.n0$transfer + ln_local.n0$percen_asphalted_roads + ln_local.n0$secondary, data=ln_local.n0)

#Robust SEs
lm.dv4.n0$newse<-vcovHC(lm.dv4.n0, type="HC1")
DV4.n0<-coeftest(lm.dv4.n0,lm.dv4.n0$newse)

###############################################################################

## DV5: CONSTITUENCY ##

# H_1: Delegates with safe seats speak less
apsr.n.4<-na.omit(data.frame(ln_con, age, male, south, party, degree, incumbent_dummy, centralnominated, fulltime, percentage, pop_stacked, transfer, percen_asphalted_roads, secondary, as.factor(career2), pci_id, vcp))

lm.dv5.1<-lm(ln_con ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary, data=apsr.n.4)

# H_5: Career FE
lm.dv5.2<-lm(ln_con ~ age+ male+ degree+ south + party + incumbent_dummy + centralnominated + fulltime + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(vcp), data=apsr.n.4)

# Interaction of Fulltime*Centralnominated
lm.dv5.3<-lm(ln_con ~ age+ male+ degree+ south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(vcp), data=apsr.n.4)

#Robust SEs
lm.dv5.3$newse<-vcovHC(lm.dv5.3, type="HC1")
DV5<-coeftest(lm.dv5.3,lm.dv5.3$newse)

# Dropping 0 Questions 
ln_con.n0<-subset(apsr.n, question_count!=0)

lm.dv5.n0<-lm(ln_con.n0$ln_con ~ ln_con.n0$age+ ln_con.n0$male+ ln_con.n0$degree+ ln_con.n0$south + ln_con.n0$party + ln_con.n0$incumbent_dummy + (ln_con.n0$centralnominated*ln_con.n0$fulltime) + ln_con.n0$percentage + ln_con.n0$pop_stacked + ln_con.n0$transfer + ln_con.n0$percen_asphalted_roads + ln_con.n0$secondary, data=ln_con.n0)


#Robust SEs
lm.dv5.n0$newse<-vcovHC(lm.dv5.n0, type="HC1")
DV5.n0<-coeftest(lm.dv5.n0,lm.dv5.n0$newse)


#########################################################################
# EXTENSION
#########################################################################

# Demonstrating/visualizing justification for using a hurdle or zero-inflated model with the first dependent variable

ln.speaknum_count<-log(speaknum_count+.001)
ln.question_count<-log(question_count+.001)

p1<-ggplot(apsr.n, aes(x=ln.speaknum_count))+geom_density(colour="#5BC6E8")+ggtitle("Log DV1")+xlim(-10,10)+labs(x="Log Speech Count", y="Density")+theme_bw()
p2<-ggplot(apsr.n, aes(x=ln.question_count))+geom_density(colour="#5BC6E8")+ggtitle("Log DV2")+xlim(-10,10)+labs(x="Log Num. Questions", y="Density")+theme_bw()

multiplot(p1,p2, cols=2)

#### Estimating Hurdle and Zero-Inflated Models ####

hurdle.dv1.4star<-hurdle(speaknum_count~age + male + degree + south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(vcp), data=apsr.n, dist="negbin")
summary(hurdle.dv1.4star)

zero.dv1.4star<-zeroinfl(speaknum_count~age + male + degree + south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(vcp), data=apsr.n, dist="negbin")
summary(zero.dv1.4star)

aic.hurd<- -2*(-541.3) + 2*(33)
aic.zero<--2*(-540.9) + 2*(31)

## Zero-inflated is a better fit than the hurdle ##

#### For visualizing the hurdle/zero-inflated comparison ####
plotreg(hurdle.dv1.4star, custom.model.names="Hurdle", omit.coef="(intercept)|(theta)", mfrow=TRUE,cex=2)
plotreg(zero.dv1.4star, custom.model.names="Zero Inflated", omit.coef="(intercept)|(theta)", mfrow=TRUE,cex=2)

# Comparing fit to the data between the negative binomial with full career fixed effects and the zero-inflated with VCP: Models are indistinguishable

vuong(nb.dv1.5, zero.dv1.4star)


# Estimating a modified version of nb.dv1.5 with VCP dummy instead of career fixed effects (as.factor(career2)) for comparison with zero inflated model

nb.dv1.4star<-glm.nb(speaknum_count~age + male + degree + south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(vcp), data=apsr.n, control=glm.control(maxit=1000))

# New Vuong test suggests that the zero-inflated model is a better fit
vuong(nb.dv1.4star, zero.dv1.4star)



#### Re-binning career fixed effects to estimate fuller zero-inflated model ####

apsr.n$careerfx<-NA
apsr.n$careerfx[apsr.n$career==1 | apsr.n$career==2 | apsr.n$career==5 | apsr.n$career==12 | apsr.n$career==21] <-1
apsr.n$careerfx[apsr.n$career==3 | apsr.n$career==4 | apsr.n$career==6] <-2
apsr.n$careerfx[apsr.n$career==7 | apsr.n$career==9 | apsr.n$career==13 | apsr.n$career==15] <-3
apsr.n$careerfx[apsr.n$career==11 | apsr.n$career==14 | apsr.n$career==16 | apsr.n$career==18 | apsr.n$career==20] <-4
apsr.n$careerfx[apsr.n$career==8 | apsr.n$career==17 | apsr.n$career==19] <-5


apsr.n$careerfx2<-NA
apsr.n$careerfx2[apsr.n$career==1 | apsr.n$career==2 | apsr.n$career==21] <-1
apsr.n$careerfx2[apsr.n$career==3] <-2
apsr.n$careerfx2[apsr.n$career==4] <-3
apsr.n$careerfx2[apsr.n$career==5 | apsr.n$career==6 | apsr.n$career==7 | apsr.n$career==9 | apsr.n$career==13 | apsr.n$career==15] <-4
apsr.n$careerfx2[apsr.n$career==8 | apsr.n$career==17 | apsr.n$career==19] <-5
apsr.n$careerfx2[apsr.n$career==11 | apsr.n$career==14 | apsr.n$career==16 | apsr.n$career==18] <-6
apsr.n$careerfx2[apsr.n$career==12] <-7


apsr.n$careerfx3<-NA
apsr.n$careerfx3[apsr.n$career==1 | apsr.n$career==2 |apsr.n$career==3 | apsr.n$career==4 | apsr.n$career==5 |apsr.n$career==6 | apsr.n$career==12 | apsr.n$career==21] <-0
apsr.n$careerfx3[apsr.n$career==7 |apsr.n$career==8 |apsr.n$career==9 |apsr.n$career==11 |apsr.n$career==13 |apsr.n$career==14 |apsr.n$career==15 |apsr.n$career==16 | apsr.n$career==17 | apsr.n$career==18 | apsr.n$career==19 | apsr.n$career==20] <-1

nb.dv1.6<-glm.nb(speaknum_count ~ age+ male+ degree+ south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(careerfx), data=apsr.n, control=glm.control(maxit=1000))

zero.dv1.6<-zeroinfl(speaknum_count~age + male + degree + south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(careerfx), data=apsr.n, dist="negbin", model=TRUE)

aic.zero6<--2*(-535.1) + 2*(39)
aic.zero6

vuong(nb.dv1.6, zero.dv1.6)

nb.dv1.6.2<-glm.nb(speaknum_count ~ age+ male+ degree+ south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(careerfx2), data=apsr.n, control=glm.control(maxit=1000))

zero.dv1.6.2<-zeroinfl(speaknum_count~age + male + degree + south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + as.factor(careerfx2), data=apsr.n, dist="negbin", model=TRUE)

aic.zero6.2<--2*(-536.8) + 2*(43)
aic.zero6.2

nb.dv1.6.3<-glm.nb(speaknum_count ~ age+ male+ degree+ south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + careerfx3, data=apsr.n, control=glm.control(maxit=1000))

zero.dv1.6.3<-zeroinfl(speaknum_count~age + male + degree + south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + careerfx3, data=apsr.n, dist="negbin")
#### Within-sample leave-one-out cross-validation for negative binomial vs. zero-inflated models ####

## Models with no career fixed effects (no career 2, no VCP) ##
zero.dv1.star<-zeroinfl(speaknum_count~age + male + degree + south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary, data=apsr.n, dist="negbin")

cv.err1<-cv.glm(apsr.n, zero.dv1.star, K=492)

nb.dv1.star<-glm.nb(speaknum_count ~ age+ male+ degree+ south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary, data=apsr.n, control=glm.control(maxit=1000))

cv.err2<-cv.glm(apsr.n, nb.dv1.star, K=492)


## Models with re-binned career fixed effects ##

cv.err.nb6<-cv.glm(apsr.n, nb.dv1.6, K=492)
cv.err.z6<-cv.glm(apsr.n, zero.dv1.6, K=492)

cv.err.z6.3<-cv.glm(apsr.n, zero.dv1.6.3, K=492)


# Cross-validation by hand
model.data<-zero.dv1.6$model
colnames(model.data)<-c("speaknum_count", "age", "male", "degree", "south", "party", "incumbent_dummy", "centralnominated", "fulltime", "percentage", "pop_stacked", "transfer", "percen_asphalted_roads", "secondary", "careerfx")
seq<-c(1:203,205:249, 251:492)
for(i in seq){
	f<-model.data[-i,]
	zero<-zeroinfl(speaknum_count~age + male + degree + south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + careerfx, data=f, dist="negbin")
	apsr.n$crossval[i]<-predict(zero, newdata=model.data[i,], type="response")
}

apsr.n$crossval[204]<-speaknum_count[204]
apsr.n$crossval[250]<-speaknum_count[250]

cv.delta<-(1/492)*sum((speaknum_count - apsr.n$crossval)^2)
cv.delta.2<-(1/492)*sum(abs(speaknum_count -apsr.n$crossval))


model.data.nb<-nb.dv1.6$model
colnames(model.data.nb)<-c("speaknum_count", "age", "male", "degree", "south", "party", "incumbent_dummy", "centralnominated", "fulltime", "percentage", "pop_stacked", "transfer", "percen_asphalted_roads", "secondary", "careerfx")
seq<-c(1:203,205:249, 251:492)
for(i in seq){
	g<-model.data.nb[-i,]
	nb<-glm.nb(speaknum_count ~ age+ male+ degree+ south + party + incumbent_dummy + (centralnominated*fulltime) + percentage + pop_stacked + transfer + percen_asphalted_roads + secondary + careerfx, data=g, control=glm.control(maxit=1000))
	apsr.n$crossval.nb[i]<-predict(nb, newdata=model.data.nb[i,], type="response")
}

cv.delta.nb<-(1/492)*sum((speaknum_count - apsr.n$crossval.nb)^2)
cv.delta.nb.2<-(1/492)*sum(abs(speaknum_count -apsr.n$crossval.nb))

#########################################################################
# FIGURES & GRAPHICS 
#########################################################################
redblue.palette<-c("#FF0A28","#5BC6E8","#464747")
redblue.2<-c("#FF0A28","#5BC6E8")
blues.palette<-c("#5BC6E8", "#0FB2E4", "#464747")
bluegray.palette<-c("#0C7A39", "#4D9440", "#464747")


stargazer(nb.dv1.1.se, nb.dv1.2.se, nb.dv1.3.se, nb.dv1.4.se, nb.dv1.5.se, align=TRUE)

stargazer(nb.dv2.1.se, nb.dv2.2.se, nb.dv2.3.se, nb.dv2.4.se, nb.dv2.5.se, align=TRUE)


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

glmm.coefs<-glmm.nb.dv2.5.1$b
glmm.ses<-glmm.nb.dv2.5.1$stdbeta
glmm.pval<-c(.8125, .5626, .0158, .1531, .4814, .5511, .4609, .5771, .000000067, .9420, .9717, .0229, .4304, .4397, .3656, .0409, .0335, .5020, .0083, .1210, .7197, .0460, .0849, .0436, .5149, .0071, .9994, .9989, .0019, .0047, .0011)

dv1<-plotreg(nb.dv1.5.se, omit.coef="(Intercept)|(careerfx)")
dv2<-plotreg(nb.dv2.5, omit.coef="(Intercept)|(careerfx)", override.coef=glmm.coefs, override.se=glmm.ses, override.pval=glmm.pval)
dv3<-plotreg(DV3)
dv4<-plotreg(DV4)
dv5<-plotreg(DV5)

dvnames<-c("Num Speeches", "Num Questions", "Critical Questions", "Local", "Constituency") 

dvs<-c("Critical Questions", "Local", "Constituency")

varnames<-c("Intercept", "Age", "Male", "Highest Degree", "South VN", "Party", "Incumbent", "Central Nominated", "Full-Time", "% Votes 2007", "2007 Population (mil)", "Transfers", "% Asphalted Roads", "% Secondary Educated", "VCP", "Central Nominated*Full-Time")

varnames2<-c("Intercept", "Age", "Male", "Highest Degree", "South", "Party", "Incumbent", "Central Nom.", "Full-Time", "% Votes 2007", "2007 Pop. (mil)", "Transfers", "% Asphalt Roads", "% Secondary Edu", "VCP", "Centralnom*FT")

varnames3<-c("Intercept", "Age", "Male", "Highest Degree", "South", "Party", "Incumbent", "Central Nom.", "Full-Time", "% Votes 2007", "2007 Pop. (mil)", "Transfers", "% Asphalt Roads", "% Secondary Edu", "careerfx2", "careerfx3", "careerfx4", "careerfx5", "careerfx7", "careerfx8", "careerfx9", "careerfx11", "careerfx12", "careerfx13", "careerfx14", "careerfx15", "careerfx17", "careerfx18", "careerfx19", "careerfx21", "Centralnom*FT")

varnames.zero<-c("Count: Intercept", "Count: Age", "Count: Male", "Count: Highest Degree", "Count: South", "Count: Party", "Count: Incumbent", "Count: Central Nom.", "Count: Full-Time", "Count: % Votes 2007", "Count: 2007 Pop. (mil)", "Count: Transfers", "Count: % Asphalt Roads", "Count: % Secondary Edu", "Count: VCP", "Count: Centralnom*FT", "Count: Log(theta)", "Zero: Intercept", "Zero: Age", "Zero: Male", "Zero: Highest Degree", "Zero: South", "Zero: Party", "Zero: Incumbent", "Zero: Central Nom.", "Zero: Full-Time", "Zero: % Votes 2007", "Zero: 2007 Pop. (mil)", "Zero: Transfers", "Zero: % Asphalt Roads", "Zero: % Secondary Edu", "Zero: VCP", "Zero: Centralnom*FT")

varnames.zero.2<-c("Count: Intercept", "Count: Age", "Count: Male", "Count: Highest Degree", "Count: South", "Count: Party", "Count: Incumbent", "Count: Central Nom.", "Count: Full-Time", "Count: % Votes 2007", "Count: 2007 Pop. (mil)", "Count: Transfers", "Count: % Asphalt Roads", "Count: % Secondary Edu", "Careerfx2", "Careerfx3", "Careerfx4", "Careerfx5", "Count: Centralnom*FT", "Count: Log(theta)", "Zero: Intercept", "Zero: Age", "Zero: Male", "Zero: Highest Degree", "Zero: South", "Zero: Party", "Zero: Incumbent", "Zero: Central Nom.", "Zero: Full-Time", "Zero: % Votes 2007", "Zero: 2007 Pop. (mil)", "Zero: Transfers", "Zero: % Asphalt Roads", "Zero: % Secondary Edu", "Careerfx2", "Careerfx3", "Careerfx4", "Careerfx5", "Zero: Centralnom*FT")


plotreg(nb.dv1.5.se, omit.coef="(Intercept)|(career2)", custom.note="", cex=1.8, custom.model.names="Num. Speak", custom.coef.names=varnames3, lwd.inner=5, lwd.outer=3, signif.light = "#FF787E", signif.medium = "#FA3C4F", signif.dark = "#FF0A28")

plotreg(nb.dv2.5, omit.coef="(Intercept)|(career2)", override.coef=glmm.coefs, override.se=glmm.ses, override.pval=glmm.pval, custom.note="*Variance-Covariance Matrix not Positive Definite", cex=1.8, custom.model.names="Num. Questions*", custom.coef.names=varnames3, lwd.inner=5, lwd.outer=3, signif.light = "#FF787E", signif.medium = "#FA3C4F", signif.dark = "#FF0A28")

plotreg(DV3, omit.coef="(Intercept)|(career2)", custom.note="", cex=1.8, custom.model.names="Critical Questions", custom.coef.names=varnames2, lwd.inner=5, lwd.outer=3, signif.light = "#FF787E", signif.medium = "#FA3C4F", signif.dark = "#FF0A28")

plotreg(DV4, omit.coef="(Intercept)|(career2)", custom.note="", cex=1.8, custom.model.names="Local", custom.coef.names=varnames2, lwd.inner=5, lwd.outer=3, signif.light = "#FF787E", signif.medium = "#FA3C4F", signif.dark = "#FF0A28")

plotreg(DV5, omit.coef="(Intercept)|(career2)", custom.note="", cex=1.8, custom.model.names="Constituency", custom.coef.names=varnames2, lwd.inner=5, lwd.outer=3, vertical.lines=TRUE, signif.light = "#FF787E", signif.medium = "#FA3C4F", signif.dark = "#FF0A28")

plotreg(l=list(DV3, DV4, DV5), file="DVs3-5.pdf", omit.coef="(Intercept)|(career2)", custom.note="", cex=1.8, custom.model.names=dvs, custom.coef.names=varnames2, lwd.inner=5, lwd.outer=3)



###### ZERO INFLATED MODELS ######

plotreg(zero.dv1.4star, omit.coef="(Intercept)|(Log)", custom.note="", cex=1.8, custom.model.names="Zero-Inflated: Num. Speach", custom.coef.names=varnames.zero, lwd.inner=5, lwd.outer=3)

plotreg(zero.dv1.4star, omit.coef="(Intercept)|(Log)", custom.note="", cex=1.8, custom.model.names="Zero-Inflated: Num. Speach", custom.coef.names=varnames.zero, lwd.inner=5, lwd.outer=3, signif.light = "#FF787E", signif.medium = "#FA3C4F", signif.dark = "#FF0A28")

plotreg(zero.dv1.6, omit.coef="(Intercept)|(Log)|(Careerfx)", custom.note="", cex=1.8, custom.model.names="Zero-Inflated: Num. Speach, New Fixed Effects", custom.coef.names=varnames.zero.2, lwd.inner=5, lwd.outer=3, signif.light = "#FF787E", signif.medium = "#FA3C4F", signif.dark = "#FF0A28", vertical.lines=FALSE)

###### VISUALIZING THE DATA ######

ggplot(apsr.n, aes(speaknum_count, fill = fulltime)) + geom_histogram() + scale_x_log10() + facet_grid(fulltime ~ ., margins=TRUE, scales="free_y")+scale_fill_manual(values=redblue.palette, name="Full Time")+theme_bw()+labs(x="Num. Speak", y="Count")

ggplot(apsr.n, aes(speaknum_count, fill = centralnominated)) + geom_histogram() + scale_x_log10() + facet_grid(centralnominated ~ ., margins=TRUE, scales="free_y")+scale_fill_manual(values=redblue.palette, name="Central Nominated")+theme_bw()+labs(x="Num. Speak", y="Count")

##### PREDICTED OUTCOMES #####

X.nb1<-expand.grid(mean(age),mean(male), mean(degree), min(south), max(party), max(incumbent_dummy), min(centralnominated):max(centralnominated), min(fulltime):max(fulltime), mean(percentage), mean(pop_stacked), mean(transfer), mean(percen_asphalted_roads), mean(secondary),min(vcp):max(vcp))

X.nb1<-cbind(X.nb1,(X.nb1[,7]*X.dv1[,8]))

colnames(X.nb1)<-c("age", "male", "degree", "south", "party", "incumbent_dummy", "centralnominated", "fulltime", "percentage", "pop_stacked", "transfer", "percen_asphalted_roads", "secondary", "vcp", "centralnominated*fulltime")

X.nb1$phat<-predict(nb.dv1.4star, newdata=X.nb1, type="response")

ggplot(X.nb1, aes(x = factor(centralnominated), y = phat, fill = factor(fulltime))) + geom_bar(stat="identity", position=position_dodge())+ facet_wrap(~vcp) + labs(x = "Central Nominated", y = "Predicted Speaking")+ theme_bw()+scale_fill_manual(values=c("#FF0A28","#5BC6E8"), name="Full Time")+ggtitle("VCP")+ylim(c(0,2))



X.dv1<-expand.grid(mean(age),mean(male), mean(degree), min(south), max(party), max(incumbent_dummy), min(centralnominated):max(centralnominated), min(fulltime):max(fulltime), mean(percentage), mean(pop_stacked), mean(transfer), mean(percen_asphalted_roads), mean(secondary),min(vcp):max(vcp))

X.dv1<-cbind(X.dv1,(X.dv1[,7]*X.dv1[,8]))

colnames(X.dv1)<-c("age", "male", "degree", "south", "party", "incumbent_dummy", "centralnominated", "fulltime", "percentage", "pop_stacked", "transfer", "percen_asphalted_roads", "secondary", "vcp", "centralnominated*fulltime")

X.dv1$phat<-predict(zero.dv1.4star, newdata=X.dv1, type="response")

	  
ggplot(X.dv1, aes(x = factor(centralnominated), y = phat, fill = factor(fulltime))) + geom_bar(stat="identity", position=position_dodge())+ facet_wrap(~vcp) + labs(x = "Central Nominated", y = "Predicted Speaking")+ theme_bw()+scale_fill_manual(values=c("#FF0A28","#5BC6E8"), name="Full Time")+ggtitle("VCP")+ylim(c(0,2))











