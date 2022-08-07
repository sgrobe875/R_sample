# Import data
d<-read.csv("/Users/sarahgrobe/Desktop/STAT229/Data/Hypoxia.csv")

### Courtesy of Maria:
d$neur_num<-0
d$neur_num[d$neurological.outcome=="moderate"]<-1
d$neur_num[d$neurological.outcome=="severe"]<-2
d$neur_num[d$neurological.outcome=="death"]<-3
# bin neurological outcome into yes=1/no=0
d$neur_bin<-1
d$neur_bin[d$neurological.outcome=="healthy"]<-0
###


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




m1<-glm(neur_bin~copep6.imp+copep12.imp+copep24.imp+copep48.imp+copep72.imp+nse6.imp+nse12.imp+nse24.imp+nse48.imp+nse72.imp+apgar5.bin+apgar10.bin+lac6.imp+lac12.imp+pH6.imp+pH12.imp+target.T.imp +gestational.age+birthweight..g.+gender..1.male.+delivery.mode, data=d2, family="binomial")
summary(m1)

