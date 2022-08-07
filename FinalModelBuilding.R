#########################################
##     Run data_cleaning.R first!      ##
#########################################


# Import packages
library(MASS)
library(ResourceSelection)





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






