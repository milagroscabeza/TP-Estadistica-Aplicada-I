#-----------TRABAJO PRACTICO E.A.I-----------

#----------PARTE A----------

datos<-as.matrix(datos[,1])

x<-as.call(datos)

rosa<-"lightsalmon"
azul<-"steelblue4"

n<-100

hist(datos,freq=FALSE,breaks=15,col=rosa,main='Densidad de X estimada',ylab='Densidad estimada',
     xlab = 'X: Registros de máxima presión [Bar]',xlim=c(170,245))
lines(density(datos),lwd=3,col=azul)

#calculo la media y la varianza muestral
media_muestral<-mean(datos, na.rm=FALSE)
varianza_muestral<-var(datos)

#calculo el coeficiente de asimetría de fisher
library(moments)
coeficiente_asimetria_fisher<-skewness(datos)

#calculo los cuartiles
cuartiles<-quantile(datos, prob=c(0,0.25,0.5,0.75,1))

#realizo el boxplot
boxplot(datos, horizontal=TRUE,ylim=c(165,250),main='Boxplot de losregistros de máxima presión [Bar]', 
        col=rosa)

datos_ordenados<-sort(c(datos))

maximo<-max(c(datos))
minimo<-min(c(datos))

#realizo el boxplot sobre el histograma

layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))

par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(datos , main="Boxplot e histograma",horizontal=TRUE , ylim=c(165,250), xaxt="n" , col=azul , frame=F)
par(mar=c(4, 3.1, 1.1, 2.1))
hist(datos , breaks=15 , col=rosa  , main="" , xlab="Registros de máxima presión [Bar]",xlim=c(165,250))
abline(v=quantile(datos,0),lwd=3,col=azul)
abline(v=quantile(datos,0.25),lwd=3,col=azul)
abline(v=quantile(datos,0.5),lwd=3,col=azul)
abline(v=quantile(datos,0.75),lwd=3,col=azul)
abline(v=quantile(datos, 1),lwd=3,col=azul)

#estimo los parámetros de la lognormal

parametro_m<-((1/n)*sum(log(datos)))

parametro_D2<-((1/n)*sum((log(datos)-parametro_m)^(2)))

parametro_D<-sqrt(parametro_D2)

#calculo la logverosimilitud de la lognormal

logverosimilitud_ln<-sum(dlnorm(datos, meanlog = parametro_m, sdlog = parametro_D, log = TRUE))

#ploteo la densidad estimada y la densidad lognormal

hist(datos,freq=FALSE,breaks=15,col=rosa,main='Densidad de X estimada + Densidad Lognormal',
     ylab='Densidad estimada',xlab = 'X: Registros de máxima presión [Bar]',xlim=c(170,245),ylim=c(0,0.04))

curve(dlnorm(x, meanlog = parametro_m, sdlog = parametro_D, log = FALSE),lwd=3,add= TRUE, col=azul)

lines(density(datos),lwd=1.8,col='black')

#estimo los parámetros de la gumbel del máximo

parametro_beta<-((sqrt(6)/pi)*sqrt(varianza_muestral))

parametro_theta<-(media_muestral-0.5772157*parametro_beta)


#calculo la logverosimilitud de la gumbel del máximo

library("ordinal")

logverosimilitud_gumbel<-sum(dgumbel(datos, location = parametro_theta, scale = parametro_beta, log = TRUE, max = TRUE))

#ploteo la densidad estimada y la densidad gumbel del máximo

hist(datos,freq=FALSE,breaks=15,col=rosa, main='Densidad de X estimada + densidad Gumbel del Máximo',
     ylab='Densidad estimada', xlab = 'X: Registros de máxima presión',xlim=c(170,245),ylim=c(0,0.04))

curve(dgumbel(x, location = parametro_theta, scale = parametro_beta, log = FALSE, 
              max = TRUE),lwd=3,add= TRUE, col=azul)

lines(density(datos),lwd=1.8,col='black')

#estimo los parámetros de la gamma

parametro_alpha<-(((media_muestral)^(2))/varianza_muestral)

parametro_gammabeta<-(varianza_muestral/media_muestral)

#calculo la logverosimilitud de la gamma

logverosimilitud_gamma<-sum(dgamma(datos, shape=parametro_alpha, 
                                   scale = parametro_gammabeta, log = TRUE))

#ploteo la densidad estimada y la densidad de la gamma
hist(datos,freq=FALSE,breaks=15,col=rosa,main='Densidad de X estimada + densidad Gamma',ylab='Densidad estimada',
     xlab = 'X: Registros de máxima presión',xlim=c(170,245),ylim=c(0,0.04))

curve(dgamma(x, shape=parametro_alpha, scale = parametro_gammabeta, log=FALSE),lwd=3,add= TRUE,col=azul)

lines(density(datos),lwd=1.8,col='black')


#----------PARTE B----------
set.seed(190)
#MUESTRA DE 50
#calculo los intervalos de confianza de cada muestra y genero un data frame

t50<-qt(0.95,df=49)

simulaciones50<-list()
medias_simulaciones50<-c()
s_simulaciones50<-c()
error_simulaciones50<-c()
li_simulaciones50<-c()
ls_simulaciones50<-c()

for(i in 1:100){
  simulacion50<-sample(c(datos), 50, replace = TRUE, prob = NULL)
  simulaciones50[[i]]<-simulacion50
  medias_simulaciones50[i]<-mean(simulacion50)
  s_simulaciones50[i]<-sd(simulacion50)
  error_simulaciones50[i]<-t50*(s_simulaciones50[i]/sqrt(50))
  li_simulaciones50[i]<-medias_simulaciones50[i]-error_simulaciones50[i] 
  ls_simulaciones50[i]<-medias_simulaciones50[i]+error_simulaciones50[i]
}

int_conf_simulaciones50<-data.frame(muestra=c(1:100), 
                                    media=medias_simulaciones50,
                                    lower=li_simulaciones50, 
                                    upper=ls_simulaciones50)

#ploteo los intervalos de confianza

library(ggplot2)
ggplot() +
  geom_errorbar(data=int_conf_simulaciones50,mapping=aes(x=muestra,ymin=li_simulaciones50,ymax=ls_simulaciones50), width=1, size=1, color=rosa) + 
  geom_point(data=int_conf_simulaciones50, mapping=aes(x=muestra, y=media), size=1, shape=21, fill=azul)+
  geom_hline(yintercept = media_muestral) 

#calculo la proporción de veces que aparece la media en las muestras

contador50<-0
for(i in 1:100){
  if((media_muestral>=li_simulaciones50[i]&media_muestral<=ls_simulaciones50[i])){
    contador50=contador50+1
  }
}


#genero un histograma con las medias
hist(medias_simulaciones50,freq=FALSE,col=rosa, main='Medias para 100 muestras de 50',
     ylab='Densidad estimada',xlab='Media estimada',ylim=c(0,0.20),xlim=c(185,210))
lines(density(medias_simulaciones50),lwd=2)


#MUESTRA DE 15
#calculo los intervalos de confianza de cada muestra y genero un data frame

t15<-qt(0.95,df=14)

simulaciones15<-list()
medias_simulaciones15<-c()
s_simulaciones15<-c()
error_simulaciones15<-c()
li_simulaciones15<-c()
ls_simulaciones15<-c()

for(i in 1:100){
  simulacion15<-sample(c(datos), 15, replace = TRUE, prob = NULL)
  simulaciones15[[i]]<-simulacion15
  medias_simulaciones15[i]<-mean(simulacion15)
  s_simulaciones15[i]<-sd(simulacion15)
  error_simulaciones15[i]<-t15*(s_simulaciones50[i]/sqrt(15))
  li_simulaciones15[i]<-medias_simulaciones15[i]-error_simulaciones15[i] 
  ls_simulaciones15[i]<-medias_simulaciones15[i]+error_simulaciones15[i]
}

int_conf_simulaciones15<-data.frame(muestra=c(1:100), 
                                    media=medias_simulaciones15,
                                    lower=li_simulaciones15, 
                                    upper=ls_simulaciones15)

#ploteo los intervalos de confianza

library(ggplot2)
ggplot() + 
  geom_errorbar(data=int_conf_simulaciones15,mapping=aes(x=muestra,ymin=li_simulaciones15, ymax=ls_simulaciones15), width=1, size=1, color=rosa) + 
  geom_point(data=int_conf_simulaciones15, mapping=aes(x=muestra, y=media), size=1, shape=21, fill=azul)+
  geom_hline(yintercept = media_muestral) 

#calculo la proporción de veces que aparece la media en las muestras

contador15<-0
for(i in 1:100){
  if((media_muestral>=li_simulaciones15[i]&media_muestral<=ls_simulaciones15[i])){
    contador15=contador15+1
  }
}



#genero un histograma con las medias
hist(medias_simulaciones15,freq=FALSE,col=rosa, main='Medias para 100 muestras de 15',
     ylab='Densidad estimada',xlab='Media estimada',xlim=c(185,210),ylim=c(0,0.12))
lines(density(medias_simulaciones15),lwd=2)


