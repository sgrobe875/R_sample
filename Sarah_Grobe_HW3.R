# Sarah Grobe
# Stat 229
# Logistic Regression (Individual)


# Import relevant packages
library(ResourceSelection)
library(dplyr)





#################################
##        Data cleaning        ##
#################################


# Import data
d<-read.csv("/Users/sarahgrobe/Desktop/STAT229/Data/Hypoxia.csv")

# create variables for numeric neurological outcomes
d$neur_num<-0
d$neur_num[d$neurological.outcome=="moderate"]<-1
d$neur_num[d$neurological.outcome=="severe"]<-2
d$neur_num[d$neurological.outcome=="death"]<-3

# bin neurological outcome into healthy = 0, not healthy = 1
d$neur_bin<-1
d$neur_bin[d$neurological.outcome=="healthy"]<-0


### Imputation ##########################################################################


# The following vars didn't require imputation:

# cooling

# neurological.outcome

# gender

# gestational.age

# delivery.mode




# copep vars ###############################

# vector of the time periods
time<-c(6,12,24,48,72,168)
time2<-time[1:5]

# vector of mean copep for each time period
mean_cop<-apply(d[,3:8], 2, mean, na.rm=TRUE)

# vector of means for each outcome group
control<-rep(-1,6)
cases<-rep(-1,6)
mean_cop_bygr<-cbind(control,cases)
for (i in (1:6)){
  mean_cop_bygr[i,1]<-tapply( d[,i+2],d$neur_bin, mean, na.rm=TRUE)[1]
  mean_cop_bygr[i,2]<-tapply( d[,i+2],d$neur_bin, mean, na.rm=TRUE)[2]
}

# 1 / mean of each group gives most linear trend:
plot(time,1/(mean_cop), pch=16, main="mean 1/copeptin over time")


#create separate linear models for different outcomes
cop_control_mod1<-lm(1/mean_cop_bygr[,1]~time)
cop_control_mod2<-lm(1/mean_cop_bygr[1:5,1]~time2)
cop_cases_mod1<-lm(1/mean_cop_bygr[,2]~time)
cop_cases_mod2<-lm(1/mean_cop_bygr[1:5,2]~time2)
par(mfrow=c(1,1))
plot(time, 1/(mean_cop_bygr[,1]),col='blue', pch=16,cex=1.2, main='mean 1/copeptin by outcome')
points(time, 1/(mean_cop_bygr[,2]), pch=18, cex=1.2, col='red')
abline(cop_control_mod1,col='blue', lwd=2)
abline(cop_control_mod2,col='blue',lty=2, lwd=2)
abline(cop_cases_mod1,col='red', lwd=2)
abline(cop_cases_mod2,col='red',lty=2, lwd=2)
legend("topleft", legend=c("Controls", "Cases"),
       col=c("blue", "red"), pch=c(16,18), cex=1.2)

# impute values from linear model 1/copeptin~time for each observation, 
# ignore 168h time point, since there are too many missing values
cop_vals<-d[1:85,3:8]
inv_cop<-t(1/(cop_vals))
inv_cop<-data.frame(cbind(time, inv_cop))

new.time <- data.frame(time=c(6,12,24,48,72,168))

for (k in (2:86)) {
  m<-lm(inv_cop[,k]~inv_cop$time)
  cop_vals$pred6[k-1]<-1/(predict(m, newdata=new.time)[1])
  cop_vals$pred12[k-1]<-1/(predict(m, newdata=new.time)[2])
  cop_vals$pred24[k-1]<-1/(predict(m, newdata=new.time)[3])
  cop_vals$pred48[k-1]<-1/(predict(m, newdata=new.time)[4])
  cop_vals$pred72[k-1]<-1/(predict(m, newdata=new.time)[5])
  cop_vals$slope[k-1]<-m$coefficients[2]
}

# impute as necessary
d$copep6.imp1 <- cop_vals$pred6
d$copep6.imp <- d$copeptin.6h..pmol.l.
d$copep6.imp[is.na(d$copep6.imp)] <- d$copep6.imp1[is.na(d$copep6.imp)]
d$copep6.imp1<-NULL

d$copep12.imp1 <- cop_vals$pred12
d$copep12.imp <- d$copeptin.12h..pmol.l.
d$copep12.imp[is.na(d$copep12.imp)] <- d$copep12.imp1[is.na(d$copep12.imp)]
d$copep12.imp1<-NULL

d$copep24.imp1 <- cop_vals$pred24
d$copep24.imp <- d$copeptin.24h..pmol.l.
d$copep24.imp[is.na(d$copep24.imp)] <- d$copep24.imp1[is.na(d$copep24.imp)]
d$copep24.imp1<-NULL

d$copep48.imp1 <- cop_vals$pred48
d$copep48.imp <- d$copeptin.48h..pmol.l.
d$copep48.imp[is.na(d$copep48.imp)] <- d$copep48.imp1[is.na(d$copep48.imp)]
d$copep48.imp1<-NULL

d$copep72.imp1 <- cop_vals$pred72
d$copep72.imp <- d$copeptin.72h..pmol.l.
d$copep72.imp[is.na(d$copep72.imp)] <- d$copep72.imp1[is.na(d$copep72.imp)]
d$copep72.imp1<-NULL


# NSE vars ############################

# vector of the time periods
time<-c(6,12,24,48,72,168)

# vector of means for each outcome group
control<-rep(-1,6)
cases<-rep(-1,6)
mean_nse_bygr<-cbind(control,cases)
for (i in (1:6)){
  mean_nse_bygr[i,1]<-tapply( d[,i+8],d$neur_bin, mean, na.rm=TRUE)[1]
  mean_nse_bygr[i,2]<-tapply( d[,i+8],d$neur_bin, mean, na.rm=TRUE)[2]
}

# create separate linear models for different outcomes
nse_control_mod1<-lm(mean_nse_bygr[,1]~time)
nse_control_mod2<-lm(mean_nse_bygr[1:5,1]~time2)
nse_cases_mod1<-lm(mean_nse_bygr[,2]~time)
nse_cases_mod2<-lm(mean_nse_bygr[1:5,2]~time2)
par(mfrow=c(1,1))
plot(time, mean_nse_bygr[,1],col='blue', pch=16,cex=1.5, main='mean NSE by outcome')
points(time, mean_nse_bygr[,2], pch=18, cex=1.5, col='red')
abline(nse_control_mod1,col='blue', lwd=2)
abline(nse_control_mod2,col='blue',lty=2, lwd=2)
abline(nse_cases_mod1,col='red', lwd=2)
abline(nse_cases_mod2,col='red',lty=2, lwd=2)
legend("topright", legend=c("Controls", "Cases"),
       col=c("blue", "red"), pch=c(16,18), cex=1.2)



# as with copep, use nse.mod2 for time periods 1-5, not enough data for time period 6

# impute as necessary
d$nse6.imp <- d$NSE.6h..ng.ml.
d$nse6.imp[is.na(d$NSE.6h..ng.ml.) & d$neur_bin==0] <- mean_nse_bygr[1,1]
d$nse6.imp[is.na(d$NSE.6h..ng.ml.) & d$neur_bin==1] <- mean_nse_bygr[1,2]

d$nse12.imp <- d$NSE.12h..ng.ml.
d$nse12.imp[is.na(d$NSE.12h..ng.ml.) & d$neur_bin==0] <-mean_nse_bygr[2,1]
d$nse12.imp[is.na(d$NSE.12h..ng.ml.) & d$neur_bin==1] <-mean_nse_bygr[2,2]

d$nse24.imp <- d$NSE.24h..ng.ml.
d$nse24.imp[is.na(d$NSE.24h..ng.ml.) & d$neur_bin==0] <- mean_nse_bygr[3,1]
d$nse24.imp[is.na(d$NSE.24h..ng.ml.) & d$neur_bin==1] <- mean_nse_bygr[3,2]

d$nse48.imp <- d$NSE.48h..ng.ml.
d$nse48.imp[is.na(d$NSE.48h..ng.ml.) & d$neur_bin==0] <- mean_nse_bygr[4,1]
d$nse48.imp[is.na(d$NSE.48h..ng.ml.) & d$neur_bin==1] <- mean_nse_bygr[4,2]

d$nse72.imp <- d$NSE.72h..ng.ml.
d$nse72.imp[is.na(d$NSE.72h..ng.ml.) & d$neur_bin==0] <- mean_nse_bygr[5,1]
d$nse72.imp[is.na(d$NSE.72h..ng.ml.) & d$neur_bin==1] <- mean_nse_bygr[5,2]







# pH 6 and pH 12 ##############################

# both vars are pretty tightly bounded; pH 6 from 6.8 to 7.5 and pH 12 from 6.91 to 7.52
# pH 6 has a slight skew left, but pH 12 is fairly normal

hist(d$pH.6h)
hist(d$pH.12h)

# because of all of this, I think it is reasonable to simply impute with means (use group means based on outcome)
d$pH6.imp <- d$pH.6h
d$pH6.imp[is.na(d$pH.6h) & d$Sample.. < 76 & d$neur_bin==0] <- tapply( d$pH.6h,d$neur_bin, mean, na.rm=TRUE)[1]
d$pH6.imp[is.na(d$pH.6h) & d$Sample.. < 76 & d$neur_bin==1] <- tapply( d$pH.6h,d$neur_bin, mean, na.rm=TRUE)[2]

d$pH12.imp <- d$pH.12h
d$pH12.imp[is.na(d$pH.12h) & d$Sample.. < 76 & d$neur_bin==0] <- tapply( d$pH.12h,d$neur_bin, mean, na.rm=TRUE)[1]
d$pH12.imp[is.na(d$pH.12h) & d$Sample.. < 76 & d$neur_bin==1] <- tapply( d$pH.12h,d$neur_bin, mean, na.rm=TRUE)[2]


# base excess at 6h and 12h ##############################

# distributions are asymmetrical
# impute group medians

hist(d$base.excess.6h)
hist(d$base.excess.12h)

d$base6.imp <- d$base.excess.6h
d$base6.imp[is.na(d$base.excess.6h) & d$Sample.. < 76 & d$neur_bin==0] <- tapply( d$base.excess.6h,d$neur_bin, median, na.rm=TRUE)[1]
d$base6.imp[is.na(d$base.excess.6h) & d$Sample.. < 76 & d$neur_bin==1] <- tapply( d$base.excess.6h,d$neur_bin, median, na.rm=TRUE)[2]

d$base12.imp <- d$base.excess.6h
d$base12.imp[is.na(d$base.excess.12h) & d$Sample.. < 76 & d$neur_bin==0] <- tapply( d$base.excess.12h,d$neur_bin, median, na.rm=TRUE)[1]
d$base12.imp[is.na(d$base.excess.12h) & d$Sample.. < 76 & d$neur_bin==1] <- tapply( d$base.excess.12h,d$neur_bin, median, na.rm=TRUE)[2]



### lactate 6 and 12 hrs ###################

# distributions are asymmetrical
# impute group medians

hist(d$lactate.6h)
hist(d$lactate.12h)

d$lac6.imp <- d$lactate.6h
d$lac6.imp[is.na(d$lactate.6h) & d$Sample.. < 76 & d$neur_bin==0] <- tapply( d$lactate.6h,d$neur_bin, median, na.rm=TRUE)[1]
d$lac6.imp[is.na(d$lactate.6h) & d$Sample.. < 76 & d$neur_bin==1] <- tapply( d$lactate.6h,d$neur_bin, median, na.rm=TRUE)[2]

d$lac12.imp <- d$lactate.12h
d$lac12.imp[is.na(d$lactate.12h) & d$Sample.. < 76 & d$neur_bin==0] <- tapply( d$lactate.12h,d$neur_bin, median, na.rm=TRUE)[1]
d$lac12.imp[is.na(d$lactate.12h) & d$Sample.. < 76 & d$neur_bin==1] <- tapply( d$lactate.12h,d$neur_bin, median, na.rm=TRUE)[2]





### Apgar 5 min ##########################

# create bins for special attention/not
# 0 = apgar >=7 (no special attention); 1 = apgar < 7 (special attention needed)
d$apgar5.bin<-1
d$apgar5.bin[d$Apgar.5min=="7"]<-0
d$apgar5.bin[d$Apgar.5min=="8"]<-0
d$apgar5.bin[d$Apgar.5min=="9"]<-0
d$apgar5.bin[d$Apgar.5min=="10"]<-0
d$apgar5.bin[d$Apgar.5min=="nd"]<-NA

# check against neurological outcome
table(d$neurological.outcome, d$apgar5.bin)
table(d$neur_bin, d$apgar5.bin)
# all categories more likely to be 1 (although likely due to sample size)

# impute accordingly
d$apgar5.bin[is.na(d$apgar5.bin)] <- 1
d$apgar5.bin[76:85]<-NA




### Apgar 10 min ##########################

# same procedure as apgar 5 min above

d$apgar10.bin<-1
d$apgar10.bin[d$Apgar.10min=="7"]<-0
d$apgar10.bin[d$Apgar.10min=="8"]<-0
d$apgar10.bin[d$Apgar.10min=="9"]<-0
d$apgar10.bin[d$Apgar.10min=="10"]<-0
d$apgar10.bin[d$Apgar.10min == ''] <- NA
d$apgar10.bin[d$Apgar.10min=="nd"]<-NA

table(d$neurological.outcome, d$apgar10.bin)
table(d$neur_bin, d$apgar10.bin)

d$apgar10.bin[is.na(d$apgar10.bin)] <- 1
d$apgar10.bin[76:85]<-NA





### target T ################################

# does not seem to be a clear relationship between target temp and other covariates
# as a result, it will be sufficient to replace missing values with group median, since 
# the distribution is skewed:
par(mfrow=c(1,2))
hist(d$target.T.reached..h.[d$neur_bin==0])
hist(d$target.T.reached..h.[d$neur_bin==1])

d$target.T.imp <- d$target.T.reached..h.
d$target.T.imp[is.na(d$target.T.reached..h.) & d$Sample.. < 76 & d$neur_bin==0]<- tapply(d$target.T.reached..h.,d$neur_bin,median, na.rm = TRUE)[1]
d$target.T.imp[is.na(d$target.T.reached..h.) & d$Sample.. < 76 & d$neur_bin==1]<- tapply(d$target.T.reached..h.,d$neur_bin,median, na.rm = TRUE)[2]






####### Loop to check number of NA values in first 75 rows #####################
j <- 1 
counter <- 0
while (j < 76) {
  if (is.na(d$apgar5.bin[j]) | d$apgar5.bin[j] == '') {
    counter <- counter + 1
  }
  j <- j + 1
}
print(counter)
####################################################################################




# delete observations with only 1 copeptin time point

d$cop_slope<-cop_vals$slope
d<-d[!is.na(d$cop_slope),]

# designate categorical variables as factors

d$apgar10.bin<-as.factor(d$apgar10.bin)
d$apgar5.bin<-as.factor(d$apgar5.bin)
d$gender..1.male.<-as.factor(d$gender..1.male.)
d$neur_bin<-as.factor(d$neur_bin)
d$neur_num<-as.factor(d$neur_num)

# Create data subset of only variables with data in all 85 rows ##############################

d1<-d[,29:40]

# Create data subset of only first 75 rows, all subjects had cooling ##############################

d2<-d[1:75,-(2:21)]
d2<-d2[1:75,-(5:7)]





##############################################
##      Building the final group model      ##
##############################################




# Test each covariate for significance in a model with only one explanatory variable
# generate list of p-values
pvals <- list()
for(i in names(d2)[-1]){
  m<- glm(neur_bin ~ get(i),family=binomial,data= d2)
  pvals[[i]] <-summary(m)$coefficients[2,4]
}
pvals

# Use only variables with p-value<0.2 and variables initially hypothesized to predict the outcome (copeptin and nse)
# designate this set of explanatories as the maxmod (model with maximum number of parameters)
# perform step wise selection to determine best model to minimizine AIC and BIC

minmod<-glm(neur_bin~1, family= binomial, data=d2)
maxmod<-glm(neur_bin~copep6.imp+copep12.imp+copep24.imp+copep48.imp+copep72.imp+nse6.imp+nse12.imp+
              nse24.imp+nse48.imp+nse72.imp+apgar5.bin+gender..1.male., data=d2, family="binomial")

(d2_m_AIC<-stepAIC(maxmod, direction= "both")) 

# BIC is to restrictive --> no explanatories at all in the model

bic.pen <- log(nrow(d2))
(d2_m_BIC<-stepAIC(maxmod, direction= "both", k=bic.pen))

# call model obtained by minimizing AIC --> call this the base model

summary(d2_m_AIC)

basem<- glm(neur_bin ~copep24.imp + copep48.imp + copep72.imp + 
              nse6.imp + nse24.imp + apgar5.bin + gender..1.male. ,family=binomial,data= d2)
summary(basem)$coefficients




# repeat check of all other variables for significance when added to this base model 

pvals2 <- list()
for(i in names(d2)[-1]){
  m<- glm(neur_bin ~copep24.imp + copep48.imp + copep72.imp + 
            nse6.imp + nse24.imp + apgar5.bin + gender..1.male.+ get(i),family=binomial,data= d2)
  pvals2[[i]] <-summary(m)$coefficients[,4]
}

pvals2



### check for confounding from gestational age

# basem2 = base model basem, plus gestational.age
basem2 <- glm(neur_bin ~copep24.imp + copep48.imp + copep72.imp + 
                nse6.imp + nse24.imp + apgar5.bin + gender..1.male. + gestational.age ,family=binomial,data= d2)



# loop through and calculate delta beta hat percent for each coefficient in the two models
base.coeff <- basem$coefficients
base2.coeff <- basem2$coefficients


dbhp <- list()
i <- 2
while (i <= length(base.coeff)) {
  value <- 100*((base.coeff[i] - base2.coeff[i])/base2.coeff[i])
  dbhp[i] <- value
  i <- i+1
}

dbhp

# All values of dbhp are reasonably low, so no adjustments needed




# check where interactions might occur
# repeat this code moving the get(i) term onto each variable
# only find significant interaction term for nse6 with nse24

d3<-cbind(d2[,7], d2[,10:13], d2[,15], d2[,4], d2[,24])
colnames(d3)<-c("neur_bin", "copep24.imp" , "copep48.imp" , "copep72.imp" , 
                "nse6.imp" , "nse24.imp", "gender..1.male.", "apgar5.bin" )

pvals3 <- list()
for(i in names(d3)[-1]){
  m<- glm(neur_bin ~copep24.imp  + copep48.imp + copep72.imp + 
            nse6.imp + nse24.imp* get(i) + apgar5.bin + gender..1.male.,family=binomial,data= d2)
  pvals3[[i]] <-summary(m)$coefficients[,4]
}

pvals3

# model with interaction term --> name it intm

intm<- glm(neur_bin ~copep24.imp+ copep48.imp  +  copep72.imp+
             nse6.imp * nse24.imp + gender..1.male. + apgar5.bin ,family=binomial,data= d2)
summary(intm)$coefficients




# check to see if any of the individual terms are no longer significant once interaction term is included
# remove apgar5 and copep72

intm2<- glm(neur_bin ~copep24.imp+ copep48.imp +   
              nse6.imp * nse24.imp + gender..1.male. ,family=binomial,data= d2)
summary(intm2)$coefficients

# check AIC for these final two model options --> they are comparable

AIC(intm)
AIC(intm2)




### hoslem test for goodness of fit

hoslem.test(intm2$y, intm2$fitted.values, g = 10)

# large p-value indicating no evidence of lack of fit

par(mfrow = c(1,1))



#############################################
##        Individual model building        ##
#############################################



## Question 1 ##



# Outcome = neurological outcome, specifically healthy (0) or any other outcome (1)


# Hosmer-Lemeshow test & calibration curve

n.grp <- 10                     # number of groups
modelName <- intm2              # name of the model to be tested
responseVar <- d2$neur_bin      # response variable for the model being tested


n.cut <- n.grp + 1
pred.prob <- predict(modelName, type="response")
quants <- quantile(pred.prob, probs = seq(0,1, length.out = n.cut))
group.id <- cut(pred.prob, breaks = quants, include.lowest = TRUE)
HL.curve <- cbind(pred.prob = pred.prob, neur_bin = responseVar, group = group.id)
ooo <- order(pred.prob)
HL.curve <- HL.curve[ooo,]

pred.mean <- by(HL.curve[,1], HL.curve[,3], mean)
obs.mean <- by(HL.curve[,2]==2, HL.curve[,3], mean)
obs.sum.cases <- by(HL.curve[,2]==2, HL.curve[,3], sum)
exp.sum.cases <- by(HL.curve[,1], HL.curve[,3], sum)
obs.sum.ctrls <- by(HL.curve[,2]==1, HL.curve[,3], sum)
exp.sum.crtls <- by(1-HL.curve[,1], HL.curve[,3], sum)

HL.test.stat <- sum((obs.sum.cases - exp.sum.cases)^2/exp.sum.cases) + 
  sum((obs.sum.ctrls - exp.sum.crtls)^2/exp.sum.crtls)
HL.test.stat
HL.pval <- pchisq(HL.test.stat, df = n.grp-2, lower.tail = F)
HL.pval

# no evidence of inadequate fit



n.cases <- as.numeric(as.data.frame(HL.curve[HL.curve[,2]==2,1]) %>% summarise(n()))
n.ctrls <- as.numeric(as.data.frame(HL.curve[HL.curve[,2]==1,1]) %>% summarise(n()))

pred.unique <- unique(pred.prob)
pred.unique <- pred.unique[order(pred.unique)]
TP <- rep(NA, NROW(pred.unique))
FP <- rep(NA, NROW(pred.unique))
for(ii in 1:NROW(pred.unique)) {
  TP[ii] <- mean(HL.curve[HL.curve[,2]==2, 1] > pred.unique[ii])
  FP[ii] <- mean(HL.curve[HL.curve[,2]==1, 1] > pred.unique[ii])
}
case.pred <- HL.curve[HL.curve[,2]==2, 1]
ctrl.pred <- HL.curve[HL.curve[,2]==1, 1]
AUC.sum <- 0
for(ii in 1:NROW(ctrl.pred)) {
  AUC.sum <- AUC.sum + sum(case.pred > ctrl.pred[ii])
}
C.index <- AUC.sum/(n.cases*n.ctrls)
C.index

# Fairly good AUC/C index value


# Brier score
BS <- mean((HL.curve[,1]-(HL.curve[,2]==2))^2)
BS

# Also a good Brier score, matching other results




# Calibration curve
plot(pred.mean, obs.mean, xlim=c(0,1),
     ylim = c(0,1), type = "p", lwd = 2,
     xlab = 'Predictions', ylab = 'Outcomes', main = 'Calibration Curve')
abline(a = 0, b=1, lwd = 2, lty = 2)


# ROC curve
plot(c(1, FP, 0), c(1, TP, 0), lty = 1, lwd = 2, col = 'red', type = 'l',
     xlab = 'False Positive', ylab = 'True Positive', main = 'ROC Curve')
abline(a = 0, b = 1, lwd = 2, lty = 2)






## Question 2 ##




# Start with EDA
# Using same EDA from group project; the following are the relevant vars:

# Gender
round(prop.table(table(d2$neur_num,d2$gender), margin=2),3)
round(prop.table(table(d2$neur_bin,d2$gender), margin=2),3)


# apgar5
round(prop.table(table(d2$neur_num,d2$apgar5.bin), margin=2),3)
round(prop.table(table(d2$neur_bin,d2$apgar5.bin), margin=2),3)


# copep168 
# We chose not to use this in model analysis since it has so many missing values (62/81 are missing)


# NSE 168
# As with copep168, we chose not to use this in model analysis since it has so many missing values (also 62)




# Weaker relationships to outcome variable, but still could be relevant:

# Delivery mode
round(prop.table(table(d$neur_num,d$delivery.mode), margin=2),3)
round(prop.table(table(d$neur_bin,d$delivery.mode), margin=2),3)


# pH6
boxplot(d2$pH6.imp ~ d2$neur_bin, main = 'Distribution of pH6 by Neurological Outcome',
        xlab = 'Neurological Outcome (0 = healthy, 1 = not healthy)', ylab = 'pH')


# gestational age
par(mfrow = c(1,2))
boxplot(d2$gestational.age~d2$neur_bin, main="Gestational Age by Neur. Outcome", 
        xlab = 'Neurological Outcome', ylab = 'Gestational Age (weeks)')

boxplot(d2$gestational.age~d2$neur_num, main="Gestational Age by Neur. Outcome", 
        xlab = 'Neurological Outcome', ylab = 'Gestational Age (weeks)')
par(mfrow = c(1,1))

# Healthier babies seem to generally have higher gestational age; however, there is a lot of
# spread in the data, so this appears to only be a weak/moderate relationship.


# apgar10
round(prop.table(table(d2$neur_num,d2$apgar10.bin), margin=2),3)
round(prop.table(table(d2$neur_bin,d2$apgar10.bin), margin=2),3)



# Drop apgar10, since related to apgar5, which will be in the model regardless:
round(prop.table(table(d2$apgar10.bin,d2$apgar5.bin), margin=2),3)
round(prop.table(table(d2$apgar5.bin,d2$apgar10.bin), margin=2),3)

# apgar5 has a closer tie to the outcome, and if we know apgar5, we have a pretty good idea of 
# apgar10, as shown by the above tables


# The removal of apgar10 from consideration leaves us with two definite variables, included
# in every model:    gender, apgar5

# We also have three possible variables, which will be chosen to be included in the coming models
# as they are built:    delivery.mode, pH6.imp, gestational.age




# Next step: Look at how these vars relate to each other and the two "definites," and determine what the
# best three models to test will be!





### delivery.mode and pH6 ###

# new data set of the three variables (we will reuse this in the next part), with missing values filtered out
plotData <- d2 %>% 
  select(delivery.mode, pH6.imp, gestational.age) %>% 
  filter(delivery.mode != '')

# remove missing value as a level of the variable delivery.mode
plotData$delivery.mode <- as.character(plotData$delivery.mode)
plotData$delivery.mode <- as.factor(plotData$delivery.mode)

# plot the data to look at the relationship
plot(plotData$delivery.mode, plotData$pH6.imp, main = 'Distribution of pH6 Values by Delivery Mode',
     xlab = 'Delivery Mode', ylab = 'pH6')

# density plots for each level
ECS <- plotData %>% filter(delivery.mode == 'ECS')
vag <- plotData %>% filter(delivery.mode == 'vaginal')

plot(density(vag$pH6.imp), col = 'blue',
     main = 'Density Plots of pH6, by Delivery Mode\n(Vaginal = blue, ECS = red)',
     xlab = 'pH6'); lines(density(ECS$pH6.imp), col = 'red')

summary(vag$pH6.imp); summary(ECS$pH6.imp)
sd(vag$pH6.imp); sd(ECS$pH6.imp)

# It would seem that vaginal deliveries have slightly higher pH6, but this isn't an incredibly
# strong relationship; looks even weaker in density plots than boxplots





### delivery.mode and gestational.age ###
plot(plotData$delivery.mode, plotData$gestational.age, main = 'Gestational Age by Delivery Mode',
     xlab = 'Delivery Mode', ylab = 'Gestational Age (weeks)')

# density plots for each level
plot(density(vag$gestational.age), col = 'blue',
     main = 'Density Plots of Gestational Age, by Delivery Mode\n(Vaginal = blue, ECS = red)',
     xlab = 'Gestational Age (weeks)'); lines(density(ECS$gestational.age), col = 'red')

summary(vag$gestational.age); summary(ECS$gestational.age)
sd(vag$gestational.age); sd(ECS$gestational.age)

# While the range is the same and spread in general is nearly identical, there does seem to be a
# relationship between gestational age & delivery mode; namely, higher gestational ages are more
# likely to be vaginal deliveries, while lower gestational ages are more likely to be ECS.




### pH6 and gestational age ###
plot(d2$gestational.age, d2$pH6.imp, main = 'pH6 by Gestational Age', xlab = 'Gestational Age (weeks)', 
     ylab = 'pH6') ;abline(lm(d2$pH6.imp~d2$gestational.age), col = 'red')

summary(lm(d2$pH6.imp~d2$gestational.age))

# Very weak correlation; no apparent relationship between the two variables

# This means having both of these variables in the model will likely offer us more information than one
# of the other pairs analyzed above, since there is little to no relationship/dependency/redundancy





### Test significance of all interaction terms, like we did in model construction for group project

d4 <- d2 %>% select(neur_bin, gender..1.male., apgar5.bin, gestational.age, pH6.imp, delivery.mode)

pvals <- list()
for(i in names(d4)[-1]){
  m<- glm(neur_bin ~ gender..1.male.  + apgar5.bin + gestational.age + pH6.imp + delivery.mode* get(i),
          family="binomial", data= d4)
  pvals[[i]] <-summary(m)$coefficients[,4]
}

pvals

# Do not appear to be any significant interaction terms; therefore, none will be included in final models
# All p-values are rather large.




### Creating the models ###
# Using above EDA and analyses as a guide:

# uses two vars (gestational.age and pH6.imp) that are not at all related
model1 <- glm(neur_bin ~ gender..1.male. + apgar5.bin + gestational.age + pH6.imp, data = d2, family = 'binomial')
summary(model1)



# uses two vars (gestational.age and delivery.mode) with strongest relationship to outcome (from EDA)
model2 <- glm(neur_bin ~ gender..1.male. + apgar5.bin + gestational.age + delivery.mode, data = d2, family = 'binomial')
summary(model2)



# Model with two definites, plus strongest from EDA (gestational.age); uses three params for better ratio to events
model3 <- glm(neur_bin ~ gender..1.male. + apgar5.bin + gestational.age, data = d2, family = 'binomial')
summary(model3)




## Analysis of fit for model 1 ##



# Hosmer-Lemeshow

n.grp <- 10                     # number of groups
modelName <- model1             # name of the model to be tested
responseVar <- d2$neur_bin      # response variable for the model being tested



n.cut <- n.grp + 1
pred.prob <- predict(modelName, type="response")
quants <- quantile(pred.prob, probs = seq(0,1, length.out = n.cut))
group.id <- cut(pred.prob, breaks = quants, include.lowest = TRUE)
HL.curve <- cbind(pred.prob = pred.prob, neur_bin = responseVar, group = group.id)
ooo <- order(pred.prob)
HL.curve <- HL.curve[ooo,]

pred.mean <- by(HL.curve[,1], HL.curve[,3], mean)
obs.mean <- by(HL.curve[,2]==2, HL.curve[,3], mean)
obs.sum.cases <- by(HL.curve[,2]==2, HL.curve[,3], sum)
exp.sum.cases <- by(HL.curve[,1], HL.curve[,3], sum)
obs.sum.ctrls <- by(HL.curve[,2]==1, HL.curve[,3], sum)
exp.sum.crtls <- by(1-HL.curve[,1], HL.curve[,3], sum)

HL.test.stat <- sum((obs.sum.cases - exp.sum.cases)^2/exp.sum.cases) + 
  sum((obs.sum.ctrls - exp.sum.crtls)^2/exp.sum.crtls)
HL.test.stat
HL.pval <- pchisq(HL.test.stat, df = n.grp-2, lower.tail = F)
HL.pval

# no evidence of inadequate fit



# AUC of ROC curve
n.cases <- as.numeric(as.data.frame(HL.curve[HL.curve[,2]==2,1]) %>% summarise(n()))
n.ctrls <- as.numeric(as.data.frame(HL.curve[HL.curve[,2]==1,1]) %>% summarise(n()))

pred.unique <- unique(pred.prob)
pred.unique <- pred.unique[order(pred.unique)]
TP <- rep(NA, NROW(pred.unique))
FP <- rep(NA, NROW(pred.unique))
for(ii in 1:NROW(pred.unique)) {
  TP[ii] <- mean(HL.curve[HL.curve[,2]==2, 1] > pred.unique[ii])
  FP[ii] <- mean(HL.curve[HL.curve[,2]==1, 1] > pred.unique[ii])
}
case.pred <- HL.curve[HL.curve[,2]==2, 1]
ctrl.pred <- HL.curve[HL.curve[,2]==1, 1]
AUC.sum <- 0
for(ii in 1:NROW(ctrl.pred)) {
  AUC.sum <- AUC.sum + sum(case.pred > ctrl.pred[ii])
}
C.index <- AUC.sum/(n.cases*n.ctrls)
C.index

# okay AUC; not much better than random


# Brier score
BS <- mean((HL.curve[,1]-(HL.curve[,2]==2))^2)
BS

# BS matches above; okay, could be better



# Calibration curve, ROC curve, TP/FP over threshold
plot(pred.mean, obs.mean, xlim=c(0,1),
     ylim = c(0,1), type = "p", lwd = 2,
     xlab = 'Predictions', ylab = 'Outcomes', main = 'Calibration Curve')
abline(a = 0, b=1, lwd = 2, lty = 2)


plot(c(1, FP, 0), c(1, TP, 0), lty = 1, lwd = 2, col = 'red', type = 'l',
     xlab = 'False Positive', ylab = 'True Positive', main = 'ROC Curve')
abline(a = 0, b = 1, lwd = 2, lty = 2)


plot(c(0, pred.unique, 1), c(1, TP, 0), lwd = 2, lty = 1, xlab = 'Predictions',
     ylab = 'Sensitivity and Specificity', xlim = c(0,1),
     main = 'Precision-Recall Curve', type = 'l') # Sensitivity
lines(c(0, pred.unique, 1), 1-c(1, FP,0), lwd = 2, lty = 2, col = 'red') # Specificity
legend(0.35, 0.2, lty = c(1,2), lwd = 2, col = c('black', 'red'), legend = c('Recall (Sensitivity)',
                                                                             'Precision (Specificity)'), bty = 'n')





## Analysis of fit for model 2 ##



# Hosmer-Lemeshow

n.grp <- 10                     # number of groups
modelName <- model2             # name of the model to be tested
responseVar <- d2$neur_bin      # response variable for the model being tested



n.cut <- n.grp + 1
pred.prob <- predict(modelName, type="response")
quants <- quantile(pred.prob, probs = seq(0,1, length.out = n.cut))
group.id <- cut(pred.prob, breaks = quants, include.lowest = TRUE)
HL.curve <- cbind(pred.prob = pred.prob, neur_bin = responseVar, group = group.id)
ooo <- order(pred.prob)
HL.curve <- HL.curve[ooo,]

pred.mean <- by(HL.curve[,1], HL.curve[,3], mean)
obs.mean <- by(HL.curve[,2]==2, HL.curve[,3], mean)
obs.sum.cases <- by(HL.curve[,2]==2, HL.curve[,3], sum)
exp.sum.cases <- by(HL.curve[,1], HL.curve[,3], sum)
obs.sum.ctrls <- by(HL.curve[,2]==1, HL.curve[,3], sum)
exp.sum.crtls <- by(1-HL.curve[,1], HL.curve[,3], sum)

HL.test.stat <- sum((obs.sum.cases - exp.sum.cases)^2/exp.sum.cases) + 
  sum((obs.sum.ctrls - exp.sum.crtls)^2/exp.sum.crtls)
HL.test.stat
HL.pval <- pchisq(HL.test.stat, df = n.grp-2, lower.tail = F)
HL.pval

# no evidence of inadequate fit



# AUC of ROC curve
n.cases <- as.numeric(as.data.frame(HL.curve[HL.curve[,2]==2,1]) %>% summarise(n()))
n.ctrls <- as.numeric(as.data.frame(HL.curve[HL.curve[,2]==1,1]) %>% summarise(n()))

pred.unique <- unique(pred.prob)
pred.unique <- pred.unique[order(pred.unique)]
TP <- rep(NA, NROW(pred.unique))
FP <- rep(NA, NROW(pred.unique))
for(ii in 1:NROW(pred.unique)) {
  TP[ii] <- mean(HL.curve[HL.curve[,2]==2, 1] > pred.unique[ii])
  FP[ii] <- mean(HL.curve[HL.curve[,2]==1, 1] > pred.unique[ii])
}
case.pred <- HL.curve[HL.curve[,2]==2, 1]
ctrl.pred <- HL.curve[HL.curve[,2]==1, 1]
AUC.sum <- 0
for(ii in 1:NROW(ctrl.pred)) {
  AUC.sum <- AUC.sum + sum(case.pred > ctrl.pred[ii])
}
C.index <- AUC.sum/(n.cases*n.ctrls)
C.index

# okay AUC; not much better than random


# Brier score
BS <- mean((HL.curve[,1]-(HL.curve[,2]==2))^2)
BS

# BS matches above; okay, could be better



# Calibration curve, ROC curve, TP/FP over threshold
plot(pred.mean, obs.mean, xlim=c(0,1),
     ylim = c(0,1), type = "p", lwd = 2,
     xlab = 'Predictions', ylab = 'Outcomes', main = 'Calibration Curve')
abline(a = 0, b=1, lwd = 2, lty = 2)


plot(c(1, FP, 0), c(1, TP, 0), lty = 1, lwd = 2, col = 'red', type = 'l',
     xlab = 'False Positive', ylab = 'True Positive', main = 'ROC Curve')
abline(a = 0, b = 1, lwd = 2, lty = 2)


plot(c(0, pred.unique, 1), c(1, TP, 0), lwd = 2, lty = 1, xlab = 'Predictions',
     ylab = 'Sensitivity and Specificity', xlim = c(0,1),
     main = 'Precision-Recall Curve', type = 'l') # Sensitivity
lines(c(0, pred.unique, 1), 1-c(1, FP,0), lwd = 2, lty = 2, col = 'red') # Specificity
legend(0.35, 0.2, lty = c(1,2), lwd = 2, col = c('black', 'red'), legend = c('Recall (Sensitivity)',
                                                                             'Precision (Specificity)'), bty = 'n')







## Analysis of fit for model 3 ##



# Hosmer-Lemeshow

n.grp <- 9                      # number of groups
modelName <- model3             # name of the model to be tested
responseVar <- d2$neur_bin      # response variable for the model being tested



n.cut <- n.grp + 1
pred.prob <- predict(modelName, type="response")
quants <- quantile(pred.prob, probs = seq(0,1, length.out = n.cut))
group.id <- cut(pred.prob, breaks = quants, include.lowest = TRUE)
HL.curve <- cbind(pred.prob = pred.prob, neur_bin = responseVar, group = group.id)
ooo <- order(pred.prob)
HL.curve <- HL.curve[ooo,]

pred.mean <- by(HL.curve[,1], HL.curve[,3], mean)
obs.mean <- by(HL.curve[,2]==2, HL.curve[,3], mean)
obs.sum.cases <- by(HL.curve[,2]==2, HL.curve[,3], sum)
exp.sum.cases <- by(HL.curve[,1], HL.curve[,3], sum)
obs.sum.ctrls <- by(HL.curve[,2]==1, HL.curve[,3], sum)
exp.sum.crtls <- by(1-HL.curve[,1], HL.curve[,3], sum)

HL.test.stat <- sum((obs.sum.cases - exp.sum.cases)^2/exp.sum.cases) + 
  sum((obs.sum.ctrls - exp.sum.crtls)^2/exp.sum.crtls)
HL.test.stat
HL.pval <- pchisq(HL.test.stat, df = n.grp-2, lower.tail = F)
HL.pval

# no evidence of inadequate fit



# AUC of ROC curve
n.cases <- as.numeric(as.data.frame(HL.curve[HL.curve[,2]==2,1]) %>% summarise(n()))
n.ctrls <- as.numeric(as.data.frame(HL.curve[HL.curve[,2]==1,1]) %>% summarise(n()))

pred.unique <- unique(pred.prob)
pred.unique <- pred.unique[order(pred.unique)]
TP <- rep(NA, NROW(pred.unique))
FP <- rep(NA, NROW(pred.unique))
for(ii in 1:NROW(pred.unique)) {
  TP[ii] <- mean(HL.curve[HL.curve[,2]==2, 1] > pred.unique[ii])
  FP[ii] <- mean(HL.curve[HL.curve[,2]==1, 1] > pred.unique[ii])
}
case.pred <- HL.curve[HL.curve[,2]==2, 1]
ctrl.pred <- HL.curve[HL.curve[,2]==1, 1]
AUC.sum <- 0
for(ii in 1:NROW(ctrl.pred)) {
  AUC.sum <- AUC.sum + sum(case.pred > ctrl.pred[ii])
}
C.index <- AUC.sum/(n.cases*n.ctrls)
C.index

# okay AUC; not much better than random


# Brier score
BS <- mean((HL.curve[,1]-(HL.curve[,2]==2))^2)
BS

# BS matches above; okay, could be better



# Calibration curve, ROC curve, TP/FP over threshold
plot(pred.mean, obs.mean, xlim=c(0,1),
     ylim = c(0,1), type = "p", lwd = 2,
     xlab = 'Predictions', ylab = 'Outcomes', main = 'Calibration Curve')
abline(a = 0, b=1, lwd = 2, lty = 2)


plot(c(1, FP, 0), c(1, TP, 0), lty = 1, lwd = 2, col = 'red', type = 'l',
     xlab = 'False Positive', ylab = 'True Positive', main = 'ROC Curve')
abline(a = 0, b = 1, lwd = 2, lty = 2)


plot(c(0, pred.unique, 1), c(1, TP, 0), lwd = 2, lty = 1, xlab = 'Predictions',
     ylab = 'Sensitivity and Specificity', xlim = c(0,1),
     main = 'Precision-Recall Curve', type = 'l') # Sensitivity
lines(c(0, pred.unique, 1), 1-c(1, FP,0), lwd = 2, lty = 2, col = 'red') # Specificity
legend(0.35, 0.2, lty = c(1,2), lwd = 2, col = c('black', 'red'), legend = c('Recall (Sensitivity)',
                                                                             'Precision (Specificity)'), bty = 'n')




