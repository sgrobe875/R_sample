# Sarah Grobe
# STAT 229
# Group 7
# Logistic Regression (Individual)

library(modEvA)
library(pROC)
library(dplyr)
library(dynpred)
library(rms)
library(ResourceSelection)

###############################################################
##   RUN data_cleaning.R AND FinalModelBuilding.R FIRST!!!   ##
###############################################################

## Question 1 ##

# Outcome = neurological outcome, specifically healthy (0) or any other outcome (1)

# Hosmer-Lemeshow
modEvAmethods('getBins')     # bin.method options
HLfit(intm2, bin.method = "n.bins")
HLfit(intm2, bin.method = "round.prob")
HLfit(intm2, bin.method = "prob.bins")
HLfit(intm2, bin.method = "size.bins")
HLfit(intm2, bin.method = "quantiles")
# Not sure which one we want, but all have the same result

# p-value (n.bins) = 0.763, so there is not significant evidence of poor fit

hoslem.test(intm2$y, intm2$fitted.values, g = 10)





# ROC curve
roc <- roc(intm2$y, intm2$fitted.values, plot = TRUE, print.auc = TRUE, main = "ROC Plot for Group Model")
# AUC = 0.799 which I think is good? Figure out how to interpret that tho
# From book, 0.7-0.8 is acceptable and 0.8-0.9 is excellent; our result right on the cusp of these,
# so it's definitely good!

# Summarize these together, but I'm pretty sure they both show that this is a pretty good fit.



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




#### Interactions?????

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


# Since no interactions are necessary, use the last remaining combination of relevant vars
model4 <- glm(neur_bin ~ gender..1.male. + apgar5.bin + pH6.imp + delivery.mode, data = d2, family = 'binomial')
summary(model4)



# Model with two definites, plus strongest from EDA (gestational.age)
model3 <- glm(neur_bin ~ gender..1.male. + apgar5.bin + gestational.age, data = d2, family = 'binomial')
summary(model3)


### Evaluate the models ###

roc1 <- roc(model1$y, model1$fitted.values, plot = TRUE, print.auc = TRUE, main = "ROC Plot for Model 1")
roc2 <- roc(model2$y, model2$fitted.values, plot = TRUE, print.auc = TRUE, main = "ROC Plot for Model 2")
roc3 <- roc(model3$y, model3$fitted.values, plot = TRUE, print.auc = TRUE, main = "ROC Plot for Model 3")

# Should be noted that the group model is better than any of these three

# Out of these three, model2 is the best, with AUC = 0.659

# This AUC is considered to be poor, according to the text; so this model kind of works, but it isn't
# very helpful or significant

# Also, look into that bit at the end that dips underneath; what does that mean?





####################################################################################################################
# Describe all the steps of your analysis including any assumptions you make. Answer the questions, and discuss    #
# the choices you made in your analysis, any limitations, and how other choices may have altered your conclusions. #
####################################################################################################################







## Paramita's way of doing the ROC curve and such

# pred.prob <- predict(intm2, type="response")
# HL.curve <- cbind(pred.prob = pred.prob, neur_bin = d2$neur_bin)
# str(HL.curve)
# head(HL.curve, 20)
# ooo <- order(pred.prob)
# HL.curve <- HL.curve[ooo,]
# str(HL.curve)
# head(HL.curve, 20)
# n.cut <- 5
# HL.curve <- cbind(HL.curve, group = rep(1:n.cut, each = NROW(HL.curve)/n.cut))
# str(HL.curve)
# head(HL.curve, 20)
# pred.mean <- by(HL.curve[,1], HL.curve[,3], mean)
# obs.mean <- by(HL.curve[,2]==2, HL.curve[,3], mean)
# plot(pred.mean, obs.mean, xlim = c(0,1), ylim = c(0,1), type = 'b')

n.grp <- 10
n.cut <- n.grp + 1
pred.prob <- predict(intm2, type="response")
quants <- quantile(pred.prob, probs = seq(0,1, length.out = n.cut))
group.id <- cut(pred.prob, breaks = quants, include.lowest = TRUE)
table(group.id, useNA = 'a')
HL.curve <- cbind(pred.prob = pred.prob, neur_bin = d2$neur_bin, group = group.id)
ooo <- order(pred.prob)
HL.curve <- HL.curve[ooo,]
head(HL.curve, 20)

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
# indicates adequate fit

n.cases <- 33
n.ctrls <- 38
pred.unique <- unique(pred.prob)
pred.unique <- pred.unique[order(pred.unique)]
TP <- rep(NA, NROW(pred.unique))
FP <- rep(NA, NROW(pred.unique))
for(ii in 1:NROW(pred.unique)) {
  TP[ii] <- mean(HL.curve[HL.curve[,2]==2, 1] > pred.unique[ii])
  FP[ii] <- mean(HL.curve[HL.curve[,2]==1, 1] > pred.unique[ii])
}
head(cbind(FP, TP))
case.pred <- HL.curve[HL.curve[,2]==2, 1]
ctrl.pred <- HL.curve[HL.curve[,2]==1, 1]
AUC.sum <- 0
for(ii in 1:NROW(ctrl.pred)) {
  AUC.sum <- AUC.sum + sum(case.pred > ctrl.pred[ii])
}
C.index <- AUC.sum/(n.cases*n.ctrls)
C.index
# Pretty good! 0.799
# Matches what I had before!!!





# Brier score
BS <- mean((HL.curve[,1]-(HL.curve[,2]==2))^2)
BS
# Decently low



# Calibration curve, ROC curve, TP/FP over threshold
plot(pred.mean, obs.mean, xlim=c(0,1),
     ylim = c(0,1), type = "p", lwd = 2,
     xlab = 'Predictions', ylab = 'Outcomes', main = 'Calibration Curve')
abline(a = 0, b=1, lwd = 2, lty = 2)
#
plot(c(1, FP, 0), c(1, TP, 0), lty = 1, lwd = 2, col = 'red', type = 'l',
     xlab = 'False Positive', ylab = 'True Positive', main = 'ROC Curve')
abline(a = 0, b = 1, lwd = 2, lty = 2)
#
plot(c(0, pred.unique, 1), c(1, TP, 0), lwd = 2, lty = 1, xlab = 'Predictions',
     ylab = 'Sensitivity and Specificity', xlim = c(0,1),
     main = 'Precision-Recall Curve', type = 'l') # Sensitivity
lines(c(0, pred.unique, 1), 1-c(1, FP,0), lwd = 2, lty = 2, col = 'red') # Specificity
legend(0.35, 0.2, lty = c(1,2), lwd = 2, col = c('black', 'red'), legend = c('Recall (Sensitivity',
               'Precision (Specificity)'), bty = 'n')
#