#required packages
library(readr)
library(stats)
library(tseries)
library(ggplot2)
library(forecast)
library(TTR)
library(ggthemes)
library(gganimate)
library(magick)
library(astsa)
library(lmtest)
library(TSA)
library(MTS)
library(vars)
library(urca)
library(fUnitRoots)
library(fpp2)
library(dplyr)
#-----------------------------------

#loading the data
data=Arctic_Death_Spiral_data
ice.extent=na.omit((data$ice.extent))
length(ice.extent)
tail(ice.extent)
#--------------------------

#making a ts obejct
ice.series=ts(ice.extent,start=c(2006,1),frequency = 12)
class(ice.series)
autoplot(ts(ice.extent,start=c(2006,1),frequency = 12)) +geom_line(color="blue4", size=0.7)+
  ggtitle("Time Series Plot")+xlab("Year")+ylab("ice_extent (in sq km)")+theme(plot.title = element_text(hjust = 0.5))+ 
 theme_economist(base_size = 14) + scale_color_economist()
animated_plot <- ggplot(data, aes(x=date, y=extent1)) +
  geom_line(color="blue4", size=0.7) + geom_point(size=2) + transition_reveal(date)+ theme(plot.title = element_text(hjust = 0.9)) + theme_economist(base_size = 14) + scale_color_economist() +xlab("Year")+ylab("Ice Cover (sq km)")
animation <- animate(animated_plot, nframes=50, renderer=magick_renderer())
animation
#----------------------------------------

#checking if additive/multiplicative and (if) transformation 
ets(ice.series) #ANA reflects addtive model
BoxCox.lambda(ice.series)#lamda= 1 (no transformation)
#----------------------------------------

#components of ts 
decompose=decompose(ice.series)
plot(decompose)
seasonal=decompose$seasonal
plot(seasonal)
trend=decompose$trend
plot(trend) 
#----------------------------

#stationary test 
adf.test(ice.series) #stationary
ndiffs(ice.series,test="adf") #non seasonal differences
nsdiffs(ice.series) #seasonal differences
kpss.test(ice.series,null="Level") #not mean stationary 
#-------------------------------------

#checking for 1st difference 
df=diff(ice.series)
adf.test(df)
var(ice.series)>var(df)
#-------------------------------------

#test and training data 
training=window(ice.series,end=c(2017,12))
test=window(ice.series,start=c(2018,1),end=c(2020,8))
#-------------------------------------

#Holtwinters
holtwinter=HoltWinters(training,seasonal="additive")
checkresiduals(holtwinter)

#for Box.test generate residuals by org-fitted - not required 
residuals=training- holtwinter$fitted[,1]
Box.test(residuals,type="Ljung-Box")
#---------------------------------------------------

#Box Jenkins - SARMA Model
acf=acf(training,lag.max = 36,plot=F)
acf$lag <- acf$lag * 12# transform lag from years to months
plot(acf)
pacf=pacf(training,lag.max = 36,plot=F)
pacf$lag=pacf$lag*12 
plot(pacf)

sarma=auto.arima(training,ic="aic",test="adf",seasonal=T)
checkresiduals(sarma)
forecast2=forecast(sarma,h=32)

#check
autoplot(forecast2,lwd=0.73,fcol="dodgerblue4",flwd=0.73,showgap=F)+autolayer(test,series = "Test Data",color="indianred3",lwd=0.72)+theme(panel.background = element_rect(fill = 'aliceblue')) +xlab("Year")+ylab("Training series (in sq km)")
#----------------------------------------------

#Manual SARIMA
s.training=diff(training,lag=12)
ns.training=diff(s.training)

acf=acf(ns.training,lag.max = 36,plot=F)
acf$lag <- acf$lag * 12# transform lag from years to months
autoplot(acf)+ggtitle("ACF Plot")+theme(plot.title = element_text(hjust = 0.5)) #q=1,Q=1 oR q=2
pacf=pacf(ns.training,lag.max = 36,plot=F)
pacf$lag=pacf$lag*12 
autoplot(pacf)+ggtitle("PACF Plot")+theme(plot.title = element_text(hjust = 0.5)) #P=1,p=2

sarima=sarima(training, p = 2, d = 1, q = 1, P = 1, D = 1, Q = 1, S = 12)
checkresiduals(sarima2$fit)
forecast6=sarima.for(training,n.ahead=32, p = 2, d = 1, q = 1, P = 1, D = 1, Q = 1, S = 12)

#check
plot <- rbind(cbind(fortify(training), sd = 0), cbind(fortify(forecast6$pred), sd = as.numeric(forecast6$se)))
plot$upper <- plot$y + plot$sd * 1.96
plot$lower <- plot$y - plot$sd * 1.96
plot$upper2<- plot$y + plot$sd * 1.282
plot$lower2<- plot$y - plot$sd * 1.282
ggplot(plot, aes(x = x ,y = y)) + 
  geom_line() + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2,fill="deepskyblue3") + geom_ribbon(aes(ymin = lower2, ymax = upper2), alpha = 0.2,fill="deepskyblue3") +
  ylab("Training series (in sq km)") + xlab("Year")+theme(panel.background = element_rect(fill = 'aliceblue'))+
  ggtitle("Forecast from ARIMA(2,1,1)(1,1,1)[12]")+theme(plot.title = element_text(hjust = 0.5)) +
  autolayer(test,series = "Test Data",color="indianred3",lwd=0.73)
#-----------------------------------------------

#ets 
ets=ets(training)
checkresiduals(ets)
#----------------------------------------------

#BATS
bats=bats(training,use.box.cox = NULL,use.trend = NULL,use.damped.trend = NULL,seasonal.periods = 12,use.arma.errors = TRUE)
checkresiduals(bats)

#for Ljung Box Test 
resids=bats$errors
Box.test(resids,type="Ljung-Box") 

forecast8=forecast(bats,h=32)

#check
autoplot(forecast8,lwd=0.73,fcol="dodgerblue4",flwd=0.73,showgap=F)+autolayer(test,series = "Test Data",color="indianred3",lwd=0.72)+theme(panel.background = element_rect(fill = 'aliceblue')) +xlab("Year")+ylab("Training series (in sq km)")
#-------------------------------------

#Univariate checking accuracy- use MAPE to check test value 
accuracy(forecast2,ice.series) #sarma
accuracy(forecast6$pred,ice.series) #sarima - best
accuracy(forecast8,ice.series) #bats
#-----------------------------------------

#SARIMAX 
#check cross correlation
ccf=ccf(temperature,ice.series,plot=F,lag.max=24)
ccf$lag <- ccf$lag * 12
plot(ccf) #lag at -4 

#fitting 
temp=data$temperature
temperature=ts(temp,start=c(2006,1),frequency = 12)
autoplot(temperature,col="red4")+theme(panel.background = element_rect(fill = 'gold1')) +xlab("Year")+ylab("Temperature")

t.training=head(temp,length(training))
t.test=tail(temp,length(test))

t.sarimax=Arima(training,order=c(2L,1L,1L),seasonal = list(order=c(1L,1L,1L),frequency=12),xreg = t.training)
checkresiduals(t.sarimax)
forecast7=forecast::forecast(t.sarimax,h=length(t.test),xreg=t.test)

#check
autoplot(forecast7,lwd=0.73,fcol="dodgerblue4",flwd=0.73,showgap=F)+autolayer(test,series = "Test Data",color="indianred3",lwd=0.72)+theme(panel.background = element_rect(fill = 'aliceblue')) +xlab("Year")+ylab("Training series (in sq km)")

#MAPE check 
accuracy(forecast7,ice.series)#best accuracy
#------------------------------------------------

#Final Forecast 
sarimax=Arima(ice.series,order=c(2L,1L,1L),seasonal = list(order=c(1L,1L,1L),freq=12),xreg= temperature)
checkresiduals(sarimax)

#forecasting temperature first 
t.arima=Arima(temperature,order=c(4L,1L,1L),seasonal=c(2L,0L,0L))#on checking various models, best MAPE
t.forecast=forecast(t.arima,h=12)
t.temp=(t.forecast$mean)

forecast=forecast(sarimax,h=length(t.temp),xreg=t.temp)
autoplot(forecast7,lwd=0.73,fcol="dodgerblue4",flwd=0.73,showgap=F)+theme(panel.background = element_rect(fill = 'aliceblue')) +xlab("Year")+ylab("Training series (in sq km)")
#-----------------------------------------

#temperature forecast 
#checking if additive/multiplicative and (if) transformation 
ets(ice.series) #ANA reflects addtive model
BoxCox.lambda(ice.series)#lamda= 1 (no transformation)
#----------------------------------------

#components of ts 
t.decompose=decompose(temperature)
plot(t.decompose)
t.seasonal=decompose$seasonal
plot(t.seasonal)
t.trend=decompose$trend
plot(t.trend) 
#----------------------------

#stationary test 
adf.test(temperature) #not stationary
ndiffs(temperature,test="adf") #non seasonal differences
nsdiffs(temperature) #seasonal differences
kpss.test(temperature,null="Level") #not mean stationary 
#-------------------------------------

#test and training data 
t.training=window(temperature,end=c(2017,12))
t.test=window(temperature,start=c(2018,1),end=c(2020,8))
#-------------------------------------

#Holtwinters
t.holtwinter=HoltWinters(t.training,seasonal="additive")
checkresiduals(t.holtwinter)

#for Box.test generate residuals by org-fitted - not required 
t.residuals=t.training- t.holtwinter$fitted[,1]
Box.test(t.residuals,type="Ljung-Box")
#---------------------------------------------------

#Box Jenkins - SARMA Model
acf(t.training)
pacf(t.training)

t.sarma=auto.arima(t.training,ic="aic",test="adf",seasonal=F)
checkresiduals(t.sarma)
forecast2t=forecast(t.sarma,h=32)
#----------------------------------------------

#ets 
t.ets=ets(t.training)
checkresiduals(t.ets)
forecast3t=forecast(t.ets,h=32)
#----------------------------------------------

#temperature checking accuracy- use MAPE to check test value 
accuracy(forecast2t,temperature) #ets 
accuracy(forecast3t,temperature) #sarma - best MAPE 
#----------------------------------------------

