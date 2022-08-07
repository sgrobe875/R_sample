d<-read.csv("/Users/sarahgrobe/Desktop/STAT229/Data/Hypoxia.csv")



#assign numerical values to neurological outcome
d$neur_num<-0
d$neur_num[d$neurological.outcome=="moderate"]<-1
d$neur_num[d$neurological.outcome=="severe"]<-2
d$neur_num[d$neurological.outcome=="death"]<-3
#bin neurological outcome into yes=1/no=0
d$neur_bin<-1
d$neur_bin[d$neurological.outcome=="healthy"]<-0

#############################Descriptive Statistics##################################################

######################categorical/ordinal covariates##################################################

#gender is a confounder. Females are more likely to have a neurological outcome
round(prop.table(table(d$neur_num,d$gender), margin=2),3)

#delivery mode does not have a huge effect on overall probability to have a neurological outcome
round(prop.table(table(d$neur_num,d$delivery.mode), margin=2),3)

#bin apgar score into >6=healthy=0, <7=needs extra medical help=1
d$ap5_bin<-1
d$ap5_bin[d$Apgar.5min=="7"]<-0
d$ap5_bin[d$Apgar.5min=="8"]<-0
d$ap5_bin[d$Apgar.5min=="9"]<-0
d$ap5_bin[d$Apgar.5min=="10"]<-0
d$ap10_bin<-1
d$ap10_bin[d$Apgar.10min=="7"]<-0
d$ap10_bin[d$Apgar.10min=="8"]<-0
d$ap10_bin[d$Apgar.10min=="9"]<-0
d$ap10_bin[d$Apgar.10min=="10"]<-0
#apgar at 5 and at 10 min<7 has higher probability of a good neurological outcome
round(prop.table(table(d$neur_num,d$ap5_bin), margin=2),3)
round(prop.table(table(d$neur_num,d$ap10_bin), margin=2),3)


########################assess continuous covariates######################################

#copeptin
par(mfrow=c(2,3))
tvec <- c("6h", "12h","24h", "48h","72h", "168h")

for (k in (3:8)){
  boxplot(d[,k]~d$neur_bin, ylab="copeptin", 
          main=substitute(tvec[k-2]))
}
for (k in (3:8)){
  boxplot(d[,k]~d$neur_num, ylab="copeptin", 
          main=substitute(tvec[k-2]))
}

#nse
par(mfrow=c(2,3))
for (k in (9:14)){
  boxplot(d[,k]~d$neur_bin, ylab="nse", 
          main=substitute(tvec[k-8]))
}
for (k in (9:14)){
  boxplot(d[,k]~d$neur_num, ylab="nse", 
          main=substitute(tvec[k-8]))
}

#pH
tvec2<-c("6h","12h")
par(mfrow=c(1,2))
for (k in (15:16)){
  boxplot(d[,k]~d$neur_bin, ylab="pH", 
          main=substitute(tvec2[k-14]))
}
for (k in (15:16)){
  boxplot(d[,k]~d$neur_num, ylab="pH", 
          main=substitute(tvec2[k-14]))
}

#base excess
par(mfrow=c(1,2))
for (k in (17:18)){
  boxplot(d[,k]~d$neur_bin, ylab="base excess", 
          main=substitute(tvec2[k-16]))
}
for (k in (17:18)){
  boxplot(d[,k]~d$neur_num, ylab="base excess", 
          main=substitute(tvec2[k-16]))
}

#lactation
par(mfrow=c(1,2))
for (k in (19:20)){
  boxplot(d[,k]~d$neur_bin, ylab="lactation", 
          main=substitute(tvec2[k-18]))
}

for (k in (19:20)){
  boxplot(d[,k]~d$neur_num, ylab="lactation", 
          main=substitute(tvec2[k-18]))
}

#gestational age, birthweigh, target temp reached
d$gestational.age<-as.numeric(d$gestational.age)
d$birthweight..g.<-as.numeric(d$birthweight..g.)
d$target.T.reached..h.<-as.numeric(d$target.T.reached..h.)
par(mfrow=c(1,3))
boxplot(d$target.T.reached..h.~d$neur_bin, 
        main="target temp reached")
boxplot(d$birthweight~d$neur_bin,
        main="birthweight(g)")
boxplot(d$gestational.age~d$neur_bin,
        main="gestational age (weeks)")

boxplot(d$target.T.reached..h.~d$neur_num, 
        main="target temp reached")
boxplot(d$birthweight~d$neur_num,
        main="birthweight (g)")
boxplot(d$gestational.age~d$neur_num,
        main="gestational age (weeks)")




########## Copeptin correlation over time ####################

par(mfrow=c(2,3))

time<-c(6,12,24,48,72,168)
mean_cop<-apply(d[,3:8], 2, mean, na.rm=TRUE)
plot(time,mean_cop, pch=16, main="mean copeptin over time")
plot(time,log(mean_cop), pch=16, main="mean log(copeptin) over time")
plot(time,1/(mean_cop), pch=16, main="mean 1/copeptin over time")


################################# NSE correlation over time ####################################

mean_nse<-apply(d[,9:14], 2, mean, na.rm=TRUE)
plot(time, mean_nse, pch=16, main="mean NSE over time")













