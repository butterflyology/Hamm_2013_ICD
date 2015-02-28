#R script for Hamm (201X) paper on hierarchical distance sampling accepted for publication in the journal Insect Conservation and Diversity

#R code by Christopher Hamm (UC Davis), the code for figure 4 was created with help by James Fordyce

require(unmarked)
set.seed(8)

#Import the data (you will have to alter the path to suit you computer)
ds4 <- read.csv('/User/Desktop/ds_Lib_4a_Jul.csv', header=T)
head(ds4)
str(ds4) #Want to confirm that Habitat (which is the management unit the transect occurred in) is treated as a factor
summary(ds4)
sum(ds4[,2:6])#Sum observations across transects


#6th of July 2011
ds6 <- read.csv('/User/Desktop/ds_Lib_6a_Jul.csv', header=T)
head(ds6)
str(ds6)
summary(ds6)
sum(ds6[,2:6])


#7th of July 2011
ds7 <- read.csv('/User/Desktop/ds_Lib_7b_Jul.csv', header=T)
head(ds7)
str(ds7)
summary(ds7)
sum(ds7[,2:6])

#8th of July 2011
ds8 <- read.csv('/User/Desktop/ds_Lib_8b_Jul.csv', header=T)
head(ds8)
str(ds8)
summary(ds8)
sum(ds8[,2:6])


#12th of July 2011
ds12 <- read.csv('/User/Desktop/ds_Lib_12b_Jul.csv', header=T)
head(ds12)
str(ds12)
summary(ds12)
sum(ds12[,2:6])



#Create the unmarked frames
#4 July
J4 <- unmarkedFrameDS(y=cbind(ds4$D0, ds4$D1, ds4$D2, ds4$D3, ds4$D4), siteCovs=data.frame(Unit=ds4$Unit), dist.breaks=c(0,1,2,3,4,5), tlength=ds4$Length, survey='line', unitsIn='m')
summary(J4)

#6th of July 2011
J6 <- unmarkedFrameDS(y=cbind(ds6$D0, ds6$D1, ds6$D2, ds6$D3, ds6$D4), siteCovs=data.frame(Unit=ds6$Unit), dist.breaks=c(0,1,2,3,4,5), tlength=ds6$Length, survey='line', unitsIn='m')
summary(J6)

#7th of July 2011
J7 <- unmarkedFrameDS(y=cbind(ds7$D0, ds7$D1, ds7$D2, ds7$D3, ds7$D4), siteCovs=data.frame(Unit=ds7$Unit), dist.breaks=c(0,1,2,3,4,5), tlength=ds7$Length, survey='line', unitsIn='m')
summary(J7)

#8th of July 2011
J8 <- unmarkedFrameDS(y=cbind(ds8$D0, ds8$D1, ds8$D2, ds8$D3, ds8$D4), siteCovs=data.frame(Unit=ds8$Unit), dist.breaks=c(0,1,2,3,4,5), tlength=ds8$Length, survey='line', unitsIn='m')
summary(J8)

#12th of July 2011
J12 <- unmarkedFrameDS(y=cbind(ds12$D0, ds12$D1, ds12$D2, ds12$D3, ds12$D4), siteCovs=data.frame(Unit=ds12$Unit), dist.breaks=c(0,1,2,3,4,5), tlength=ds12$Length, survey='line', unitsIn='m')
summary(J12)


# Sum the total transect lengths
TLen <- matrix(c(sum(ds4$Length[1:24]), sum(ds6$Length[1:21]), sum(ds7$Length[1:22]), sum(ds8$Length[1:30]), sum(ds12$Length)), nrow=5, ncol=1, byrow=T, dimnames = list(c('4 July', '6 July', '7 July', '8 July', '12 July'), 'Transect Length'))
TLen
sum(TLen)

#Sum the number of butterflies observed
Nbfly <- matrix(c(sum(ds4[1:24, 2:6]), sum(ds6[1:21, 2:6]), sum(ds7[1:22, 2:6]), sum(ds8[1:30, 2:6]), sum(ds12[1:16, 2:6])), nrow=5, ncol=1, byrow=T, dimnames= list(c('4 July', '6 July', '7 July', '8 July', '12 July'), 'N Butterflies'))
Nbfly
sum(Nbfly)


#Histograms of distance distributions by survey date.
#Figure 3
pdf(file='Dates3.pdf', bg='white')
par(mfrow=c(1,5), mar=c(4.5,4.2,2,1.5))
hist(J4, col='grey', main='', xlab='Distance (m)', ylab='# Butterflies Observed', ylim=c(0,50), las=1, cex.lab=1.5, cex.axis=1.5)#cex to make text larger?
text(1, 50, 'A', cex=1.5)
hist(J6, col='grey', main='', xlab='Distance (m)', ylab='', ylim=c(0, 30), las=1, cex.lab=1.5, cex.axis=1.5)
text(1, 30, 'B', cex=1.5)
hist(J7, col='grey', main='', xlab='Distance (m)', ylab='', ylim=c(0, 25), las=1, cex.lab=1.5, cex.axis=1.5) #Odd distribution, but these are the data
text(1, 25, 'C', cex=1.5)
hist(J8, col='grey', main='', xlab='Distance (m)', ylab='', ylim=c(0,20), las=1, cex.lab=1.5, cex.axis=1.5)
text(1, 20, 'D', cex=1.5)
hist(J12, col='grey', main='', xlab='Distance (m)', ylab='', ylim=c(0, 7), las=1, cex.lab=1.5, cex.axis=1.5)
text(1, 7, 'E', cex=1.5)
dev.off()


#With the data imported and visualized (the data do look appropriate for HDS analysis) I fit the halfnormal and negative exponential functions. Both the hn and nexp seem biologically plausible for these data. Temperature and Relative Humidity did not vary during each survey, but the transect did occur on different management units. So I will use "Unit" as a covariate on detection and abundance. 
#Half-normal models
j4.h1 <- distsamp(~1 ~1, J4, keyfun='halfnorm', output='density', method='BFGS', unitsOut='ha', se=T)
j4.h2 <- update(j4.h1, formula = ~Unit ~1)
j4.h3 <- update(j4.h1, formula = ~1 ~Unit) 
j4.h4 <- update(j4.h1, formula = ~Unit ~ Unit)
#Negative exponential models
j4.e1 <- distsamp(~1 ~1, J4, keyfun='exp', output='density', method='BFGS', unitsOut='ha', se=T)
j4.e2 <- update(j4.e1, formula= ~Unit ~ 1)
j4.e3 <- update(j4.e1, formula= ~1 ~Unit)
j4.e4 <- update(j4.e1, formula= ~Unit ~Unit)

#figure 1
pdf(file='Fig1.pdf', bg='white')
hist(j4.h1, col='grey', xlab='Perpendicular Distance (m)', main='', ylim=c(0, 0.4), lwd=3, las=1, cex.lab=1.5, cex.axis=1.5)
abline(h=0.361, lwd=3.5, lty=2)
dev.off()

#Now I create a fitList to conduct AIC model selection and prediction
j4.list <- fitList('H ~1 ~1' = j4.h1, 'H ~Unit ~1' = j4.h2, 'H ~1 ~Unit'=j4.h3, 'H ~Unit ~Unit'=j4.h4, 'E ~1 ~1' = j4.e1, 'E ~Unit ~1' = j4.e2, 'E ~1 ~Unit' = j4.e3, 'E ~Unit ~Unit' = j4.e4)
modSel(j4.list)#j4.e2 is AIC best and the hn sister model is very close in dAIC units
backTransform(j4.e2, type='state')#Parameter estimates from the AIC best model

#I want to make sure that the AIC best model reasonably approximates the data, so I will conduct a chi-square goodness of fit test
chisq <- function(fm){
	observed <- getY(fm@data)
	expected <- fitted(fm)
	sum((observed - expected)^2/expected)
}
#1000 bootstrap pseudo replicates, this will take a while
(pb.j4.e2 <- parboot(j4.e2, statistic=chisq, nsim=1000))#The model passes the test, df=102
plot(pb.j4.e2, col='grey')

#I want to generate model weighted averages for the density estimates
newdataj4 <- data.frame(Unit=c('A', 'B', 'C', 'D'))
(pj4 <- predict(j4.list, type='state', newdata=newdataj4, appendData=T))

#Bootstrap the AIC best model to extrapolate density estimate to abundance for the entire area
est <- function(fm, A){
	D <- exp(coef(fm, type="state"))
	A <- 8.321 #total area of the JCC site 
	N <- D*A
	return(N)
}

parboot(j4.e2, stat=est, nsim=1000)#1305 (944-1798) Plotting the pdf would be nice here. Df = 102 (n-1) 



#6th of July, only Units A, B were visited
j6.h1 <- distsamp(~1 ~1, J6, keyfun='halfnorm', output='density', method='BFGS', unitsOut='ha', se=T)
j6.h2 <- update(j6.h1, formula = ~Unit ~1)
j6.h3 <- update(j6.h1, formula = ~1 ~Unit) 
j6.h4 <- update(j6.h1, formula = ~Unit ~Unit)

j6.e1 <- distsamp(~1 ~1, J6, keyfun='exp', output='density', method='BFGS', unitsOut='ha', se=T)
j6.e2 <- update(j6.e1, formula= ~Unit ~ 1)
j6.e3 <- update(j6.e1, formula= ~1 ~Unit)
j6.e4 <- update(j6.e1, formula= ~Unit ~Unit)

j6.list <- fitList('H ~1 ~1' = j6.h1, 'H ~Unit ~1' = j6.h2, 'H ~1 ~Unit'=j6.h3, 'H ~Unit ~Unit'=j6.h4, 'E ~1 ~1' = j6.e1, 'E ~Unit ~1' = j6.e2, 'E ~1 ~Unit' = j6.e3, 'E ~Unit ~Unit' = j6.e4)
modSel(j6.list)#j6.e3 is AIC best

(pb.j6.e3 <- parboot(j6.e3, statistic=chisq, nsim=1000))#Model passes
plot(pb.j6.e3, col='grey')#df = 57


newdataj6 <- data.frame(Unit=c('A', 'B', 'D'))
(pj6 <- predict(j6.list, type='state', newdata=newdataj6, appendData=T))

parboot(j6.e3, stat=est, nsim=1000)#469 (181-723), df = 57




#7 July 2011
j7.h1 <- distsamp(~1 ~1, J7, keyfun='halfnorm', output='density', method='BFGS', unitsOut='ha', se=T)
j7.h2 <- update(j7.h1, formula = ~Unit ~1)
j7.h3 <- update(j7.h1, formula = ~1 ~Unit) 
j7.h4 <- update(j7.h1, formula = ~Unit ~Unit)

j7.e1 <- distsamp(~1 ~1, J7, keyfun='exp', output='density', method='BFGS', unitsOut='ha', se=T)
j7.e2 <- update(j7.e1, formula= ~Unit ~ 1)
j7.e3 <- update(j7.e1, formula= ~1 ~Unit)
j7.e4 <- update(j7.e1, formula= ~Unit ~Unit)

j7.list <- fitList('H ~1 ~1' = j7.h1, 'H ~Unit ~1' = j7.h2, 'H ~1 ~Unit'=j7.h3, 'H ~Unit ~Unit'=j7.h4, 'E ~1 ~1' = j7.e1, 'E ~Unit ~1' = j7.e2, 'E ~1 ~Unit' = j7.e3, 'E ~Unit ~Unit' = j7.e4)
modSel(j7.list) #j7.e3 is AIC best, 
(pb.j7.e3 <- parboot(j7.e3, statistic=chisq, nsim=1000))#passed the test

newdataj7 <- data.frame(Unit=c('A', 'B', 'D'))
(j7.pred <- predict(j7.list, type='state', newdata=newdataj7, appendData=T))

(parboot(j7.e3, statistic=est, nsim=1000))#439 (143 - 699), df=44

#8 July 2011
j8.h1 <- distsamp(~1 ~1, J8, keyfun='halfnorm', output='density', method='BFGS', unitsOut='ha', se=T)
j8.h2 <- update(j8.h1, formula = ~Unit ~1)
j8.h3 <- update(j8.h1, formula = ~1 ~Unit) 
j8.h4 <- update(j8.h1, formula = ~Unit ~ Unit)

j8.e1 <- distsamp(~1 ~1, J8, keyfun='exp', output='density', method='BFGS', unitsOut='ha', se=T)
j8.e2 <- update(j8.e1, formula= ~Unit ~ 1)
j8.e3 <- update(j8.e1, formula= ~1 ~Unit)
j8.e4 <- update(j8.e1, formula= ~Unit ~Unit)

j8.list <- fitList('H ~1 ~1' = j8.h1, 'H ~Unit ~1' = j8.h2, 'H ~1 ~Unit'=j8.h3, 'H ~Unit ~Unit'=j8.h4, 'E ~1 ~1' = j8.e1, 'E ~Unit ~1' = j8.e2, 'E ~1 ~Unit' = j8.e3, 'E ~Unit ~Unit' = j8.e4)
modSel(j8.list)#j8.e2 is AIC best
(pb.j8.e2 <- parboot(j8.e2 , statistic= chisq, nsim=1000))#passed test

newdataj8 <- data.frame(Unit=c('A', 'B', 'C', 'D'))
(predict(j8.list, type='state', newdata=newdataj8, appendData=T))

(parboot(j8.e2, statistic=est, nsim=1000))#531 (396 - 798), df=53

#12 July 2011 (B,C)
j12.h1 <- distsamp(~1 ~1, J12, keyfun='halfnorm', output='density', method='BFGS', unitsOut='ha', se=T)
j12.h2 <- update(j12.h1, formula = ~Unit ~1)
j12.h3 <- update(j12.h1, formula = ~1 ~Unit) 
j12.h4 <- update(j12.h1, formula = ~Unit ~ Unit)

j12.e1 <- distsamp(~1 ~1, J12, keyfun='exp', output='density', method='BFGS', unitsOut='ha', se=T)
j12.e2 <- update(j12.e1, formula= ~Unit ~ 1)
j12.e3 <- update(j12.e1, formula= ~1 ~Unit)
j12.e4 <- update(j12.e1, formula= ~Unit ~Unit)

j12.list <- fitList('H ~1 ~1' = j12.h1, 'H ~Unit ~1' = j12.h2, 'H ~1 ~Unit'=j12.h3, 'H ~Unit ~Unit'=j12.h4, 'E ~1 ~1' = j12.e1, 'E ~Unit ~1' = j12.e2, 'E ~1 ~Unit' = j12.e3, 'E ~Unit ~Unit' = j12.e4)
modSel(j12.list)#j12.h4 is AIC best
(pb.j12.e4 <- parboot(j12.e4 , statistic= chisq, nsim=1000))#passed test

newdataj12 <- data.frame(Unit=c('B', 'C'))
(predict(j12.list, type='state', newdata=newdataj12, appendData=T))

(parboot(j12.h4, statistic=est, nsim=1000))#523 (176 - 1655), df = 14


#Bootstrap to extrapolate predictions for area (because error may not scale linearly)
est <- function(fm, A){
	D <- exp(coef(fm, type="state"))
	A <- 3.659 #Area
	N <- D*A
	return(N)
}



#Create figure 4
moneyShot<-function(x){
	where<-x[1]
	dot<-x[2]
	uc<-x[3]
	lc<-x[4]	
		segments(where,uc,where,lc, lwd=3.5)
		points(where,dot, pch=15, col='grey', cex=1.5)
}

moneyShot2<-function(x){
	where<-x[1]
	dot<-x[2]
	uc<-x[3]
	lc<-x[4]	
		segments(where,uc,where,lc, lwd=3.5)
		points(where,dot, pch=17, col='black ', cex=1.5)
}

moneyShot3<-function(x){
	where<-x[1]
	dot<-x[2]
	uc<-x[3]
	lc<-x[4]	
		segments(where,uc,where,lc, lwd=3.5)
		points(where,dot, pch=19, col='light grey', cex=1.5)
}



#Plotting AIC model averaged predictions for density
pdf(file='MSB_Dens3.pdf', bg='white')
#jpeg(file='MSB_Dens3.jpeg')
plot(3:10,seq(1,210,length.out=8),type="n", ylab='Butterflies / ha', xaxt='n', xlab='July 2011', las=1, cex.axis=1.2, cex.lab=1.2)#will need to add axis eventually
axis(1, at=c(3, 4, 5, 6, 7, 8, 9, 10), lab=expression(3, 4, 5, 6, 7, 8, 12, 13), cex.axis=1.2, cex.lab=1.2)

#4 July, 2012 (ABCD)
moneyShot(c(3.8, 142, 79, 205)) #Unit A
moneyShot2(c(4, 147, 93, 201)) #Unit B
moneyShot3(c(4.2, 147, 93, 202)) #Unit C

#6 July, 2012 (ABD)
moneyShot(c(5.9, 57, 5.6, 108)) #Unit A
moneyShot2(c(6.1, 106, 55, 158)) #Unit B


#7 July, 2012 (ABD)
moneyShot(c(6.9, 42, 4.8, 79.6)) #Unit A
moneyShot2(c(7.1, 74, 32, 115)) #Unit B

#8 July, 2012
moneyShot(c(7.8, 50, 14, 87))#Unit A
moneyShot2(c(8, 58, 32, 84)) #Unit B
moneyShot3(c(8.2, 71, 28, 113)) #Unit C

#12 July, 2012 #Note the difference in date and where it is plotted
moneyShot2(c(8.9, 64, 0, 137))
moneyShot3(c(9.1, 32, 0, 65))

#Sample sizes excluding unit D
text(4, 213, 'n = 103', cex=1.2)#n for J4
text(6, 166, 'n = 58', cex=1.2)#n for J6
text(7, 123, 'n = 45', cex=1.2)#n for J7
text(8, 121, 'n = 54', cex=1.2)#n for J8
text(9, 144, 'n = 15', cex=1.2)#n for J12

legend('topright', legend=c('Unit A', 'Unit B', 'Unit C', '95% CI'), col=c('dark grey', 'black', 'light grey', 'black', 'black'), pch=c(15, 17, 19, NA), lwd=c(NA, NA, NA, 3.5), pt.cex=c(1.4, 1.4, 1.4, NA), cex=1.2)
dev.off()




#Relationship between CI (range) and sample size, AIC best model range, (regress them?)
#If I want to test these idea right, I need to do this Ian style
CIrange <- c((1798-944), (723-181), (699-143), (798-396), (1655-176))
Bfly <- c(103, 58, 45, 54, 15)

plot(Bfly, CIrange, type='p', xlab='Buterflies Observed', ylab= '95% CI range')
lm1 <- lm(CIrange ~ Bfly) #P=0.09
summary(lm1)
par(mfrow=c(2,2))
plot(lm1)


lCI <- log(CIrange)
lBfly <- log(Bfly)
lm1a <- lm(lCI ~ Bfly)
summary(lm1a)
plot(Bfly, lCI)
abline(lm1a)


lm2 <- lm(CIrange ~ lBfly)
summary(lm2)
lm2a <- lm(lCI ~ lBfly)
summary(lm2a)

plot(lBfly, lCI) #more butterflies = lower CI, but few data points
abline(lm2a)

#Cost calculations
#MRR costs
#8 persons * 5 days * $15/hour *10 hours/day
8*5*15*10

#TM
#2 persons * 5 days * $15/hour * 5 hours/day
2*5*15*3

#Distance 
#five hours for set up (once) 
#1 person * 5 days * 5 hours/day * $15/hour + 5*15 (set up cost)
1*5*5*15+(5*15)

Cost <- matrix(c(8*5*15*10, 2*5*15*3, 1*5*5*15+(5*15)), nrow=3, ncol=1, byrow=T, dimnames=list(c('MRR', 'TM', 'Distance'), 'Dollars'))
Cost

450/6000 #cost is ~7.5%

Days <- seq(1:21)
MRR <- 4*15*8
TM <- 2*15*3
Dist <- 3*15

Mday <- MRR*Days
Tdays <- TM*Days
Ddays <- Dist*Days
#need to add cost of setting up transect 5*15
Ddays <- Ddays + 75
Ddays

#Figure 5
pdf(file='Cost1.pdf')
plot(Days, Mday, type='l', lwd=3.5, lty=1, ylab='Cost USD', ylim=c(0, 10000), yaxt='n', cex.axis=1.5, cex.lab=1.5)
axis(2, at=c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000), lab=expression('0', '1K', '2K', '3K', '4K', '5K', '6K', '7K', '8K', '9K', '10K'), las=1, cex.axis=1.5, cex.lab=1.5)
lines(Tdays, lwd=3.5, lty=3)
lines(Ddays, lwd=3.5, lty=6)
legend('topleft', legend=c('MRR', 'TM', 'HDS'), lwd=3, lty=c(1,3,6), cex=1.5)
dev.off()


