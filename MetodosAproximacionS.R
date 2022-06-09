###MetodosAproximacionS
library(actuar)
library(fitdistrplus)
library(vcd)
library(ADGofTest)
library(mixtools)
library(mclust)
library(MASS)

###Datos 
###Base de datos para pérdidas agregadas

base<-read.csv("C:/Users/Salvador/Desktop/Cursos2022-II/PerdidasAgregadas.csv")
head(base)
attach(base)

###La variable frecuencia guarda la frecuencia de reclamaciones de estas 200 pólizas

barplot(table(frecuencia),col=rainbow(15),main="Tabla de frecuencias de reclamación",col.main="darkblue")

###Ajuste modelo discreto: Poisson

disc.fit<-goodfit(frecuencia,type = "poisson", method = "ML")
plot(disc.fit)
summary(disc.fit)
disc.fit$par

###Entonces, el número medio de reclamaciones es: 4.925

###Ajuste distribución continua: montos de reclamación

###Guardando los montos de reclamación en un vector y eliminando las "NA"

severidad<-c(base[,2],base[,3],base[,4],base[,5],base[,6],base[,7],base[,8],base[,9],base[,10],base[,11],base[,12],base[,13])
severidad1<-severidad[!is.na(severidad)]
length(severidad1)
severidad1[1:50]

plot(density(severidad1))


### Posibles densidades asociadas (no todas)

descdist(severidad1,discrete=FALSE,boot=5000,obs.pch = 19, boot.col="darkviolet")

### Mi función (Recuerden, ustedes no la tienen)

source("C:/Users/Salvador/Desktop/Cursos2021-I/Riesgo2021/autodistfit.R")

fitData(severidad1,c("gamma","weibull","lognormal"),sample=1)

fG<-fitdist(severidad1,"gamma")

ks.test(severidad1,"pgamma",shape=fG$estimate[1],rate=fG$estimate[2])
ad.test(severidad1, pgamma,shape=fG$estimate[1],rate=fG$estimate[2])

plot(fG)

### Entonces: N~Poisson(4.925) y X~Gamma(forma=50.98746956,rate=0.02906342). rate=1/scale, scale=34.40751

########################################################################################################################################
########################################### APROXIMACIÓN NORMAL ########################################################################

muN<-4.925; muX<-mgamma(1,shape=fG$estimate[1],rate=fG$estimate[2]);varN<-muN;varX<-mgamma(2,shape=fG$estimate[1],rate=fG$estimate[2])-muX^2
muN;muX;varN;varX

ES<-muN*muX;VS<-muN*mgamma(2,shape=fG$estimate[1],rate=fG$estimate[2])

ES;VS

Z<-function(s){pnorm(s,ES,sqrt(VS))}

curve(Z,from=0,to=25000,main="Aproximación: Distribución normal",lwd=2,col.main="darkblue",col="magenta")

1-Z(20000)

Z1<-function(s){dnorm(s,ES,sqrt(VS))}
curve(Z1,from=-5000,to=25000,main="Aproximación: Densidad normal",lwd=2,col.main="darkred",col="darkmagenta")

Z(ES)

qnorm(0.95,ES,sqrt(VS));qnorm(0.99,ES,sqrt(VS)),qnorm(0.995,ES,sqrt(VS))

###Si no hubieramos logrado ajustar modelos a nuestros datos

meanX<-mean(severidad1);varX<-var(severidad1);meanN<-mean(frecuencia);varN<-var(frecuencia)

E.S<-meanX*meanN;V.S<-meanN*varX+varN*meanX^2

E.S;V.S

Z2<-function(s){pnorm(s,E.S,sqrt(V.S))}

curve(Z2,from=0,to=25000,main="Aproximación: Distribución normal",lwd=2,col.main="darkblue",col="magenta")

Z3<-function(s){dnorm(s,E.S,sqrt(V.S))}
curve(Z3,from=-5000,to=25000,main="Aproximación: Densidad normal",lwd=2,col.main="darkred",col="darkmagenta")

###Comparación Aproximación Normal vs. Aproximación Normal Empírica

curve(Z,from=-5000,to=25000,main="Comparaciones: Aproximación normal vs. Aproximación normal empírica",lwd=2,col.main="darkred",col="darkmagenta")
curve(Z2,from=0,to=25000,col="magenta",xlab="",ylab="",add=TRUE,lwd=2)
legend(locator(1),legend=c("Aproximación Normal","Aproximación Normal Empírica"),col=c("darkmagenta","magenta"),bty="n",lwd=2)

curve(Z1,from=-5000,to=25000,main="Comparaciones: Aproximación normal vs. Aproximación normal empírica",lwd=2,col.main="darkred",col="darkmagenta")
curve(Z3,from=0,to=25000,col="magenta",xlab="",ylab="",add=TRUE,lwd=2)
legend(locator(1),legend=c("Aproximación Normal","Aproximación Normal Empírica"),col=c("darkmagenta","magenta"),bty="n",lwd=2)

#############################################################################################################################################
############################################# APROXIMACIÓN LOG-NORMAL #######################################################################
###Mismos datos anteriores. N~Poisson(4.925) y X~Gamma(forma=50.98746956,rate=0.02906342). rate=1/scale, scale=34.40751
###El sistema de ecuaciones es 2*log(E(S))=2mu+sigma^2   y  log(E(S^2))=2mu+2sigma^2
### Ya tenemos E(S) y E(S^2)=Var(S)+E(S)^2
### Resolvamos el sistema apoyándonos en R

A=matrix(c(2,1,2,2),2,byrow=TRUE)
A
b<-c(2*log(ES),log(VS+ES^2))
sol<-solve(A,b)

Z4<-function(s){pnorm(log(s),sol[1],sqrt(sol[2]))}
curve(Z4,from=0.01,to=25000,main="Aproximación: Distribución log-normal",lwd=2,col.main="darkred",col="darkmagenta")

Z5<-function(s){dnorm(log(s),sol[1],sqrt(sol[2]))}
curve(Z5,from=0.01,to=25000,main="Aproximación: Densidad log-normal",lwd=2,col.main="darkred",col="darkmagenta")

###Aproximación log normal empírica

B=matrix(c(2,1,2,2),2,byrow=TRUE)
B
c<-c(2*log(E.S),log(V.S+E.S^2))
sol1<-solve(B,c)

Z41<-function(s){pnorm(log(s),sol1[1],sqrt(sol1[2]))}
curve(Z41,from=0.01,to=25000,main="Aproximación: Distribución log-normal",lwd=2,col.main="darkred",col="darkmagenta")

Z51<-function(s){dnorm(log(s),sol1[1],sqrt(sol1[2]))}
curve(Z51,from=0.01,to=25000,main="Aproximación: Densidad log-normal",lwd=2,col.main="darkred",col="darkmagenta")


curve(Z4,from=0.01,to=25000,main="Aproximación: Distribución log-normal",lwd=2,col.main="darkred",col="darkmagenta")
curve(Z41,from=0.01,to=25000,lwd=2,col="darkblue",add=TRUE,xlab="",ylab="")
legend(locator(1),legend=c("Aproximación Log-Normal","Aproximación Log-normal empírica"),col=c("darkmagenta","darkblue"),bty="n",lwd=2)


###Comparando con las aproximaciones normales anteriores

curve(Z,from=-5000,to=25000,main="Comparaciones: Aproximación normal vs. Aproximación normal empírica vs. Aproximación Log-normal",lwd=2,col.main="darkred",col="darkmagenta")
curve(Z2,from=0,to=25000,col="magenta",xlab="",ylab="",add=TRUE,lwd=2)
curve(Z4,from=0.01,to=25000,col="darkblue",xlab="",ylab="",add=TRUE,lwd=2)

legend(locator(1),legend=c("Aproximación Normal","Aproximación Normal Empírica","Aproximación Log-normal"),col=c("darkmagenta","magenta","darkblue"),bty="n",lwd=2)


###Aproximación gamma transladada. Observemos que la distribución de nuestros montos !!!es una gamma!!!
###Ya que la distribución discreta es una Poisson, entonces la distribución de S es Poisson compuesta
###y el sesgo de S, se calcula como lambda*E(X^3)/VS^(3/2)

mu<-ES;sigma2<-VS;tau<-meanN*mgamma(3,shape=fG$estimate[1],rate=fG$estimate[2])/VS^(3/2)
mu;sigma2;tau
alpha<-4/tau^2;theta<-sqrt(sigma2)*tau/2;k<-mu-2*sqrt(sigma2)/tau
alpha;theta;k

Z6<-function(s){pgamma(s-k,shape=alpha,scale=theta)}
curve(Z6,from=0,to=25000,main="Aproximación: Distribución gamma transladada",lwd=2,col.main="darkred",col="darkmagenta")

Z7<-function(s){dgamma(s-k,shape=alpha,scale=theta)}
curve(Z7,from=0,to=25000,main="Aproximación: Densidad gamma transladada",lwd=2,col.main="darkred",col="darkmagenta")


###Comparación de todas las aproximaciones

curve(Z,from=-5000,to=25000,main="Comparaciones: Aproximación normal vs. Aproximación normal empírica vs. Aproximación Log-normal",lwd=2,col.main="darkred",col="darkmagenta")
curve(Z2,from=0,to=25000,col="magenta",xlab="",ylab="",add=TRUE,lwd=2)
curve(Z4,from=0.01,to=25000,col="darkblue",xlab="",ylab="",add=TRUE,lwd=2)
curve(Z6,from=0,to=25000,col="darkorange",add=TRUE,lwd=2)
legend(locator(1),legend=c("Aproximación Normal","Aproximación Normal Empírica","Aproximación Log-normal","Gamma transladada"),col=c("darkmagenta","magenta","darkblue","darkorange"),bty="n",lwd=2)

###################################################################################################################################
###################################################################################################################################
### Aproximación Poisson compuesta. Mismo ejemplo que en las notas de clase originales
### Considerando los datos de la compañía aseguradora en donde se cubren 3 diferentes grupos de asegurados, 
### Utilizar la Aproximación Poisson Compuesta para el modelo individual y posteriormente con la aproximación normal encontrar FS(1900) para los 3 valores que puede tomar ?.
### La tabla correspondiente a este seguro es

#### La tabla
##          
## i   # de pólizas |Probabilidad de reclamación |Monto de Reclamación
## 1      1000      |0.05                        |10
## 2      2000      |0.10                        |5
## 3       500      |0.02                        |20

### Aproximación poisson

### El objetivo básico es calcular, para cada valor de lambda_{i} asociado al valor de q_{i} la función de densidad
### f(x)=sum (lambda_{i}/lambda*f_{B_{i}} y con ella calcula media y varianza de X, y después media y varianza de S

q1<-0.05
q2<-0.1
q3<-0.02
x<-c(10,5,20)
fB5<-c(0,1,0)
fB10<-c(1,0,0)
fB20<-c(0,0,1)
n1<-1000
n2<-2000
n3<-500
n<-c(n1,n2,n3)

#############¿Porqué funciona esta aproximación?

ln<-function(x){log(x,exp(1))}
dpois(c(0,0,0,1,1,1),lambda=c(q1,q2,q3));dbinom(c(0,0,0,1,1,1),1,c(q1,q2,q3))
dpois(c(0,0,0,1,1,1),lambda=-ln(1-c(q1,q2,q3)));dbinom(c(0,0,0,1,1,1),1,-ln(1-c(q1,q2,q3)))
dpois(c(0,0,0,1,1,1),lambda=c(q1,q2,q3)/(1-c(q1,q2,q3)));dbinom(c(0,0,0,1,1,1),1,c(q1,q2,q3)/(1-c(q1,q2,q3)))

###################################################################################################

#Metodo 1 con lambdai=qi

lambdai1<-c(q1, q2, q3)

lambda1<-sum(n*lambdai1)

f1<- c(sum(n*(lambdai1/lambda1)*fB10),sum(n*(lambdai1/lambda1)*fB5),sum(n*(lambdai1/lambda1)*fB20))

Ex1<-x*f1

Ex2_1<-(x^2)*f1

ES1<-sum(lambda1*Ex1)

VarS1<- sum(lambda1*Ex2_1)

pnorm(1900, mean=ES1, sd=sqrt(VarS1),lower.tail=TRUE)

#Metodo 2 con lambdai=-ln(1-qi)

lambdai2<- c(-ln(1-q1), -ln(1-q2), -ln(1-q3))
lambda2<- sum(n*lambdai2)
f2<- c(sum(n*(lambdai2/lambda2)*fB10),sum(n*(lambdai2/lambda2)*fB5),sum(n*(lambdai2/lambda2)*fB20))
Ex2<- x*f2
Ex2_2<- (x^2)*f2
ES2<- sum(lambda2*Ex2)
VarS2<- sum(lambda2*Ex2_2)
pnorm(1900, mean=ES2, sd=sqrt(VarS2),lower.tail=TRUE)

#Metodo 3 con lambdai=qi/(1-qi)

lambdai3<- c(q1/(1-q1), q2/(1-q2), q3/(1-q3))
lambda3<- sum(n*lambdai3)
f3<- c(sum(n*(lambdai3/lambda3)*fB10),sum(n*(lambdai3/lambda3)*fB5),sum(n*(lambdai3/lambda3)*fB20))
Ex3<- x*f3
Ex2_3<- (x^2)*f3
ES3<- sum(lambda3*Ex3)
VarS3<- sum(lambda3*Ex2_3)
pnorm(1900, mean=ES3, sd=sqrt(VarS3),lower.tail=TRUE)

#Grafica de las funciones de densidad

curve(dnorm(x,mean=ES1, sd=sqrt(VarS1)),from=0,to=3000,xlab="s", ylab=expression(f[S](s)), ylim=c(0,0.004),type="l",col="darkgreen", 
main="Comparación aproximación Poisson",col.main="darkorange",col.lab="darkorange")
curve(dnorm(x,mean=ES2, sd=sqrt(VarS2)),from=0,to=3000,xlab="s",col="darkblue",add=T)
curve(dnorm(x,mean=ES3, sd=sqrt(VarS3)),from=0,to=3000,xlab="s",col="darkred",add=T)
t<-c("Método 1", "Método 2", "Método 3")
legend(100, 0.003, paste(t), col=c("darkgreen", "darkblue", "darkred"),lty=1)

curve(pnorm(x,mean=ES1, sd=sqrt(VarS1)),from=0,to=3000,xlab="s", ylab=expression(F[S](s)), ylim=c(0,1),type="l",col="darkgreen", 
main="Comparación aproximación Poisson",col.main="darkorange",col.lab="darkorange")
curve(pnorm(x,mean=ES2, sd=sqrt(VarS2)),from=0,to=3000,xlab="s",col="darkblue",add=T)
curve(pnorm(x,mean=ES3, sd=sqrt(VarS3)),from=0,to=3000,xlab="s",col="darkred",add=T)
t<-c("Método 1", "Método 2", "Método 3")
legend(100, 0.4, paste(t), col=c("darkgreen", "darkblue", "darkred"),lty=1)


### Convoluciones
### Convolución de dos distribuciones Poisson
### Sabemos que la suma (convolución) de dos v.s.a.s Poisson (lambda) es otra Poisson con parámetro lambda1+lambda2

x<-dpois(seq(0,20,1),2)

y<-dpois(seq(0,20,1),2)

z<-convolve(x,rev(y),type="o")  ###La convolución de las dos Poisson

w<-dpois(seq(0,40,1),4)  ###La distribución teórica

v<-as.matrix(cbind(z[1:20],w[1:20]))

barplot(t(v),beside=TRUE,col=c("blue","red"),legend.text=c("Convolución","Distribución teórica"))


### Convolución de dos exponenciales 0.1 (media=10)
### Sabemos que la suma (convolución) de dos v.s.a.s exponenciales (alpha, lambda) es otra una gamma con parámetro de
### forma alpha1+alpha2 y parámetro de escala lambda
### Para usar la función convolve necesitamos discretizar estas distribuciones continuas

x1<-discretize(pexp(x,0.1), method = "rounding", from = 0, to = 20,step=0.1)

x2<-discretize(pexp(x,0.1), method = "rounding", from = 0, to = 20,step=0.1)

x3<-convolve(x1,rev(x2),type="o")

x4<-discretize(pgamma(x,2,0.1), method = "rounding", from = 0, to = 20,step=0.1)

x5<-as.matrix(cbind(x3[1:20],x4[1:20]))

barplot(t(x5),beside=TRUE,col=c("green","violet"),legend.text=c("Convolución","Distribución teórica"))

### Nuestro ejemplo de notas
################### convolucion

fr<-c(0.1,0.2,0.3,0.4)  ###Igualamos los vectores con ceros donde sea necesario
fs<-c(0,0.4,0.6,0)      ###Utilizando la definición de densidad: "cero en cualquier otro lado"

Fs<-aggregateDist("convolution", model.freq = fr, model.sev = fs)

quantile(Fs)

plot(Fs,col="blue",main="Distribución de S: Método convolución",col.main="red",col.lab="green",sub="",xlab="s",ylab=expression(F[S](s)))

win.graph()

plot(Fs, do.points = FALSE, verticals = TRUE, xlim = c(0, 5),col="darkblue",lwd=2)

CDFs<-Fs(c(0,1,2,3,4,5,6))

dfs<-diff(c(0,CDFs))

barplot(dfs,col="darkred",main="Densidad de S: Método convolución",col.main="darkblue",col.lab="darkgreen",sub="",)

par(mfrow=c(1,2))

plot(Fs,col="blue",main="Distribución de S: Método convolución",col.main="red",col.lab="green",sub="")

barplot(dfs,col="darkred",main="Densidad de S: Método convolución",col.main="darkblue",col.lab="darkgreen",sub="")


VaR(Fs)

VaR(Fs,0.995)

CTE(Fs)

CTE(Fs,0.995)


### Datos reales (Seguro dental)

### Ejemplo convolucion de distribuciones discretas (tablas)
### Obsérvese que los recorridos de X y de N no son iguales
### X no toma valor en cero y N sí. X toma valores en 9 y 10, y N, no
### Debemos de completar con CEROS en cada caso, para igualar los recorridos

### Los originales

x1<-seq(1,10,1)

x2<-c(0.150,0.200,0.250,0.125,0.075,0.050,0.050,0.050,0.025,0.025)

x<-data.frame(cbind(x2,x1))

n1<-seq(0,8,1)

n2<-c(0.05,0.10,0.15,0.20,0.25,0.15,0.06,0.03,0.01)

n<-data.frame(cbind(n2,n1))

###Los "completados"

fs<-c(0,x2)   # vector de probabilidades asociado a la severidad (obsérvese que aumentamos el valor cero)

fr<-c(n2,0,0)#vector de probabilidades asociado a la frecuencia (obsérvese que aumentamos dos ceros para que los vectores sean de igual tamaño)

Fs<-aggregateDist("convolution", model.freq = fr, model.sev = fs)

quantile(Fs)

plot(Fs,col="darkblue",main="Función de distribución seguro dental",col.main="darkred",col.lab="darkgreen",sub="",xlab="s",ylab=expression(F[S](s)))

win.graph()

plot(Fs, do.points = FALSE, verticals = TRUE,col="darkblue",main="Función de distribución seguro dental",col.main="darkred",
     col.lab="darkgreen",xlab="s",ylab=expression(F[S](s)))

CDFs<-Fs(c(seq(0,40,1)))

dfs<-diff(c(0,CDFs))

barplot(dfs,col="orange",main="Función de distribución de S via convolución",xlab="s",ylab=expression(f[S](s)),col.main="darkred",
        col.lab="darkgreen")

VaR(Fs)

VaR(Fs,0.995)

CTE(Fs)

CTE(Fs,0.995)

### Tienen sentido estos valores?. "La neta, NO".
### Cómo debería de hacerse?

Fs<-aggregateDist("convolution", model.freq = fr, model.sev = fs,x.scale=25)

quantile(Fs)

plot(Fs,col="darkblue",main="Función de distribución seguro dental",col.main="darkred",col.lab="darkgreen",sub="")

win.graph()

plot(Fs, do.points = FALSE, verticals = TRUE,col="darkblue",main="Función de distribución seguro dental",col.main="darkred",

     col.lab="darkgreen")

CDFs<-Fs(c(seq(0,1000,25)))

dfs<-diff(c(0,CDFs))

barplot(dfs,col="orange",main="Función de distribución de S via convolución",xlab="s",ylab=expression(f[S](s)),col.main="darkred",
        col.lab="darkgreen")

VaR(Fs)

VaR(Fs,0.995)

CTE(Fs)

CTE(Fs,0.995)

###Con colvolve (Pendiente)

barplot(convolve(fr,rev(fs),type="o"),col="orange",main="Función de distribución de S via convolución",xlab="s",ylab=expression(f[S](s)),col.main="darkred",
        col.lab="darkgreen")

###PANJER notas
### Ejemplo Panjer 

###X~BN(5,0.4) p=1/(1+beta)
###N~Poisson(2)

x1<-seq(0,100)
fx1<-dnbinom(x1,5,0.4)
Fs1<- aggregateDist("recursive",model.freq="poisson", model.sev=fx1,lambda=2)
1-Fs1(3)

S1<-Fs1(knots(Fs1))
fS1<-diff(c(0,S1))
plot(Fs1, main="Distribución de Reclamaciones Agregadas: Panjer",col="darkblue",col.main="red",sub="",
xlab= "Aproximación Método Recursivo (Panjer)",col.lab="darkgreen")

barplot(fS1,main="Densidad de Reclamaciones Agregadas: Panjer",col="darkgreen",col.main="red",sub="",
xlab= "Aproximación Método Recursivo (Panjer)",col.lab="darkred",ylab=expression(f[S](x)))

par(mfrow=c(1,2))

plot(Fs1, main="Distribución de Reclamaciones Agregadas: Panjer",col="darkblue",col.main="red",sub="",
xlab= "Aproximación Método Recursivo (Panjer)",col.lab="darkgreen")

barplot(fS1,main="Densidad de Reclamaciones Agregadas: Panjer",col="darkgreen",col.main="red",sub="",
xlab= "Aproximación Método Recursivo (Panjer)",col.lab="darkred",ylab=expression(f[S](x)))

VaR(Fs1)

VaR(Fs1,0.995)

CTE(Fs1)

CTE(Fs1,0.995)

###La tabla

cbind(x=seq(0,10),fS=fS1[1:11],FS=S1[1:11])

#########################################################################
###Pajer: Seguro dental

### Datos reales

x1<-seq(1,10,1)

x2<-c(0.150,0.200,0.250,0.125,0.075,0.050,0.050,0.050,0.025,0.025)

x<-data.frame(cbind(x2,x1))

n1<-seq(0,8,1)

n2<-c(0.05,0.10,0.15,0.20,0.25,0.15,0.06,0.03,0.01)

n<-data.frame(cbind(n2,n1))

fs<-c(0,x2)   # vector de probabilidades asociado a la severidad (obsérvese que aumentamos el valor cero)

fr<-c(n2,0,0)#vector de probabilidades asociado a la frecuencia (obsérvese que aumentamos dos ceros para que los vectores sean de igual tamaño)

poisfit1<-goodfit(n, type = "poisson", method = "ML")

summary(poisfit1)
plot(poisfit1)
poisfit1$par


###Curiosida'

poisfit2<-goodfit(x, type = "poisson",method = "ML")
summary(poisfit2)
plot(poisfit2)
poisfit2$par


### Modelo de perdidas agregadas utilizando Panjer

Fs2<-aggregateDist("recursive", model.freq = "poisson", model.sev=fs, lambda=3.4, x.scale=25)

plot(Fs2,main="Distribución de Reclamaciones Agregadas",col="darkblue",col.main="red",sub="",
xlab= "Aproximación Método Recursivo (Panjer)",col.lab="darkgreen")

FCDFs2<-Fs2(c(seq(0,1000,1)))
		      
fCDFs2<-diff(c(0,FCDFs2))     

barplot(fCDFs2,main="Densidad de Reclamaciones Agregadas",col="darkgreen",col.main="red",sub="",
xlab= "Aproximación Método Recursivo (Panjer)",col.lab="darkred",col.lab="darkblue",ylab=expression(f[S](x)))

VaR(Fs2)

VaR(Fs2,0.995)

CTE(Fs2)

CTE(Fs2,0.995)

###Mismas cantidades, pero con el método de convolución

VaR(Fs)

VaR(Fs,0.995)

CTE(Fs)

CTE(Fs,0.995)


###Método de redondeo
###Ejemplo notas: pareto(4,50)

fx1<-discretize(ppareto(x,shape=4,scale=50), method = "rounding", from = 0, to = 10000,step=0.9)
sum(fx1)
cbind(x<-seq(0,10),f<-fx1[1:11])

###Comparación

par(oma = c( 0, 0, 5, 0 ),col="darkblue")
par(mfrow=c(1,2))

fx<- discretize(ppareto(x, shape=4, scale=50), method="rounding",from=0,to=100,step=0.9)
v<- c(0:110)
plot(stepfun(v, diffinv(fx)),xlab="x",ylab=expression(F[X](x)),col="darkred",main="Comparación: Funciones de distribución",
     col.main="red",col.lab="darkgreen")
curve( ppareto(x, shape=4, scale=50),from=0,to=120,col="blue",add=TRUE)
t<- c("Real", "Discretizacion")
legend(20, 0.4, paste(t), lty=1,col=c("blue", "darkred"))


barplot(fx,xlab="x",ylab=expression(f[X](x)),col="green",main="Comparación: Funciones de densidad",
        col.main="red",col.lab="darkgreen")
curve(dpareto(x, shape=4, scale=50),from=0,to=120,col="blue",add=TRUE)
t<- c("Real", "Discretizacion")
legend(20, 0.04, paste(t), lty=1, col=c("blue", "green"))

title("Comparación: Pareto(4,50) continua vs. discretizada",,outer = TRUE)


### Método de simulación
###Ejemplo: Poisson-BinomialNegativa

model.freq<-expression(data=rpois(2))
model.sev<-expression(data=rnbinom(size=5,prob=0.4))
SimF<-aggregateDist("simulation",nb.simul=10000,model.freq,model.sev)

plot(SimF, main="Distribución de Reclamaciones Agregadas (simulación)",col="darkblue",col.main="red",sub="",
     xlab= "Método de simulación",col.lab="darkgreen")

plot(SimF,do.points=FALSE,verticals=TRUE,main="Distribución de Reclamaciones Agregadas (simulación)",col="darkblue",col.main="red",sub="",
     xlab="Método de simulación",col.lab="darkgreen")

val<-SimF(knots(SimF))

Simf<-diff(c(0,val))

par(mfrow=c(1,2))

barplot(Simf, main="Densidad pérdidas agregadas: Método de simulación",col.lab="darkviolet",
        sub= "Aproximación vía simulación",col="red",col.main="darkgreen",col.sub="skyblue",ylab=expression(f[S](s)))

plot(SimF, main="Distribución de Reclamaciones Agregadas (simulación)",col="darkgreen",col.main="red",sub="",xlab= "Método de simulación",col.lab="darkgreen",ylab=expression(F[S](s)))

VaR(SimF)

CTE(SimF)

###################################################################################################
###Comparando (que es gerundio) los métodos vistos utilizando los datos de seguro dental

x1<-seq(1,10,1)

x2<-c(0.150,0.200,0.250,0.125,0.075,0.050,0.050,0.050,0.025,0.025)

x<-data.frame(cbind(x2,x1))

n1<-seq(0,8,1)

n2<-c(0.05,0.10,0.15,0.20,0.25,0.15,0.06,0.03,0.01)

n<-data.frame(cbind(n2,n1))

###Los "completados"

fs<-c(0,x2) 

fr<-c(n2,0,0)

ConSD<-aggregateDist("convolution", model.freq = fr, model.sev = fs)

PanSD<-aggregateDist("recursive", model.freq = "poisson", lambda=3.4,model.sev = fs)

model.freq<-expression(data=rpois(3.4))
model.sev<-expression(data=rpois(3.7))

SimSD<-aggregateDist("simulation",nb.simul=10000,model.freq,model.sev)
SimSD1<-SimSD(25*seq(0:100))

plot(ConSD,do.points=FALSE,verticals=TRUE,main="Comparación: Convolución, Panjer y Simulación",col="darkred",col.main="red",sub="",xlab="x",col.lab="darkgreen")

plot(PanSD,do.points=FALSE,verticals=TRUE,main="",col="darkgreen", xlab=" ",add=T,sub="")

plot(SimSD,do.points=FALSE,verticals=TRUE,main="",col="darkorange",xlab="",sub="",add=T)

legend(locator(1),legend=c("Convolución","Panjer","Simulación"),bty="n",col=c("darkred","darkgreen","darkorange"),lty=1)

res<-matrix(c(25*VaR(ConSD),25*VaR(ConSD,0.995),25*CTE(ConSD),25*CTE(ConSD,0.995),25*VaR(PanSD),25*VaR(PanSD,0.995),25*CTE(PanSD),25*CTE(PanSD,0.995),25*VaR(SimSD),25*VaR(SimSD,0.995),
              25*CTE(SimSD),25*CTE(SimSD,0.995)),nrow=3,byrow=TRUE)
res

rownames(res)<-c("Convolución","Panjer","Simulación");colnames(res)<-c("VaR90","VaR95","VaR99","VaR99.5","TVaR90","TVaR95","TVaR99","TVaR99.5")
t(res)


###########################################################################
#######################Método de De Pril ##################################
### Ejemplo de las notas

r <- 4 # Valor máximo del índice de los montos de reclamaciones
m <- 4 # Valor máximo del índice de probabilidades de reclamación
X <- 30 # Valor máximo que toma x en fs(x); X=3000,6000,...,90000

###Tabla

n<-diag(c(20,14,8,24))
n

q<-c(0.02,0.012,0.05,0.013)

###Función h(i,k)

h <- function(i,k) {
hold<-0
for (j in 1:m) {
hold <- hold+n[i,j]*(q[j]/(1-q[j]))^k
}
hold <- i*((-1)^(k-1))*hold
return(hold)
}

###Cálculo de fs

fs <-seq(1,30)
f <- function(x) {
hold<-1
if (x==0) {
for (i in 1:r) {
for (j in 1:m) {
hold <- hold*((1-q[j])^n[i,j])
}
}
return(hold)
}
else
{
hold <- 0
for (i in 1:min(x,r)) {
for (k in 1:floor(x/i)) {
if (x-i*k==0) {hold <- hold + f(0)*h(i,k) }
else {hold <- hold + fs[x-i*k]*h(i,k)}
}
}
hold <- hold/x
fs[r] <- hold
return(hold)
}
}

fs0<-f(0)


###La función de densidad de S

for (i in 1:X) {fs[i] <- f(i)}

Fs<-cumsum(c(fs0,fs))

options(scipen = 999)

round(cbind(fs=c(fs0,fs),Fs),7)

par(mfrow=c(1,2))

barplot(c(fs0,fs),main="Densidad de Reclamaciones Agregadas",sub="Método de De Pril", xlab="s", 
ylab=expression(f[S](s)), xlim=c(0,30), ylim=c(0,0.3),col="red",col.main="darkblue",col.lab="darkgreen",col.sub="darkred")


#Grafica de la funcion de distribucion
plot(Fs, main="Distribucion de Reclamaciones Agregadas",sub="Método de De Pril", pch=19, xlab="x", ylab=expression(F[S](x)), xlim=c(0,30), ylim=c(0,1), 
     col="darkmagenta",col.main="darkred",col.lab="darkblue",col.sub="darkred")


########################################################################################################################################################
###MEZCLAS
###MEZCLA DE DISTRIBUCIONES

### Famoso ejemplo de la literatura

###measurements give time in minutes between eruptions of the Old Faithful geyser in Yellowstone National Park

data(faithful)
attach(faithful)

faithful

hist(waiting, main = "Tiempo entre las erupciones del Old Faithful", xlab = "Minutos", ylab = "", cex.main = 1.5, 
cex.lab = 1.5, cex.axis = 1.4,col="darkgreen",border="white",prob=T)
lines(density(waiting),col="red",lwd=2)

###Algoritmo "a mano"

### Valores para iniciar el algoritmo

s = c(0.5,50,80,5,5)

EM = function(x,s) {
pi.hat= s[1]*dnorm(x, s[2], sqrt(s[4]))/(s[1]*dnorm(x, s[2], sqrt(s[4])) +(1-s[1])*dnorm(x, s[3], sqrt(s[5])))
s[1] = mean(pi.hat)
s[2] = sum(pi.hat*x) / sum(pi.hat)
s[3] = sum((1-pi.hat)*x) / sum(1-pi.hat)
s[4] = sqrt(sum(pi.hat*(x-s[2])^2) / sum(pi.hat))
s[5] = sqrt(sum((1-pi.hat)*(x-s[3])^2) / sum((1-pi.hat)))
s
}

iter = function(x, s) {
s1 = EM(x,s)
for (i in 1:5) {
if (abs(s[i]-s1[i]) > 0.00001) {
s=s1
iter(x,s)
}
else s1
}
s1
}

x<-waiting

iter(x,s)

###Probabilidades de clasificación aposteriori

pi.hat<-function(s,x){s[1]*dnorm(x, s[2], sqrt(s[4]))/(s[1]*dnorm(x, s[2], sqrt(s[4])) +(1-s[1])*dnorm(x, s[3], sqrt(s[5])))}

prob.apos<-pi.hat(c(0.3510847, 54.2247274, 79.9173437,  5.4662678^2,  5.9875619^2),waiting)

head(prob.apos,20)

options(scipen=999)
cbind(waiting[1:20],head(prob.apos,20),1-head(prob.apos,20))

f<-function(x,s){s[1]*dnorm(x, s[2], sqrt(s[4])) +(1-s[1])*dnorm(x, s[3], sqrt(s[5]))}

x<-seq(min(waiting),max(waiting),0.01)
plot(x,f(x,c(0.3510847, 54.2247274, 79.9173437,  5.4662678^2,  5.9875619^2)),col="darkmagenta",ylab="",xlab="",lwd=2)
lines(density(waiting),col="darkblue",main="Mezcla de distribuciones normales",col.main="darkviolet",ylab="",xlab="",lwd=2)


######################################################################################################################
###Con mixtools

mezajuste<-normalmixEM(waiting, lambda = 0.5, mu = c(55, 80), sigma = 5) #Se supone que las normales en la mezcla tienen varianzas iguales

plot(mezajuste, density = TRUE, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.8, main2 = "Ajuste mezcla: tiempo entre las erupciones del Old Faithful", 
xlab2 = "Minutos")


mezajuste[c("lambda", "mu", "sigma")]


### Y ajusta bien?

x<-seq(min(waiting),max(waiting),length.out=272)

p<-mezajuste$lambda
m<-mezajuste$mu
s<-mezajuste$sigma


### Ajuste de densidad

fest<-(p[1]*dnorm(x,m[1],s[1])+p[2]*dnorm(x,m[2],s[2]))

hist(waiting,col="darkgreen",border="white",prob=T,ylim=c(0,0.045))

lines(x,fest,col="blue",lwd=2)

### Ajuste de distribución

edest<-(p[1]*pnorm(x,m[1],s[1])+p[2]*pnorm(x,m[2],s[2]))

plot(ecdf(waiting),verticals=TRUE,do.points=FALSE,col.hor="red", col.vert="bisque",ylab="")
lines(x,edest,col="darkblue")

###Prueba formal de bondad de ajuste

mix.norm<-function(x, mu, sigma, pmix) {
  pmix[1]*pnorm(x,mu[1],sigma[1])+pmix[2]*pnorm(x,mu[2],sigma[2])
}

ks.test(waiting,mix.norm,mu=m,sigma=s,pmix=p)
ad.test(waiting,mix.norm,mu=m,sigma=s,pmix=p)



### OTRA FORMA

mix.norm <- function(x, mu, sigma, pmix) {
  pmix[1]*pnorm(x,mu[1],sigma[1])+pmix[2]*pnorm(x,mu[2],sigma[2])
}

ks.test(sort(waiting),mix.norm,mu=m,sigma=s,pmix=p)
ad.test(sort(waiting),mix.norm,mu=m,sigma=s,pmix=p)


### Y si no suponemos varianzas iguales???

mezajuste1<-normalmixEM(waiting, lambda = 0.5, mu = c(55, 80), sigma =c(5,5))

plot(mezajuste1, density = TRUE, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.8, main2 = "Tiempo entre las erupciones del Old Faithful", 
xlab2 = "Minutos")


mezajuste1[c("lambda", "mu", "sigma")]

###Sin dar valores iniciales

mezajuste2<-normalmixEM(waiting, lambda = NULL, mu = NULL, sigma =NULL)

plot(mezajuste2, density = TRUE, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.8, main2 = "Tiempo entre las erupciones del Old Faithful", 
xlab2 = "Minutos")


mezajuste2[c("lambda", "mu", "sigma")]


###
### Y si utilizamos alguna distribucion mas acorde al hecho que los datos son positivos???

ajuste3<-gammamixEM(waiting,lambda = c(0.5,0.5), alpha = NULL, beta = NULL, k = 2, epsilon = 1e-08, maxit = 1000, maxrestarts=20,verb = FALSE)

plot(ajuste3, , whichplots = 1, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.8, main2 = "Tiempo entre las erupciones del Old Faithful", 
xlab2 = "Minutos")

### No hay gráfica ajustada

peso<-ajuste3$lambda
alfa<-ajuste3$gamma.pars[1,]
beta<-ajuste3$gamma.pars[2,]


summary(ajuste3)# No hay resumen ni graficacion automatica

### Mixture gamma

gmixture<-function(x,alpha,beta,lambda) {

  lambda<-lambda/sum(lambda)

  K<-length(lambda)

  tmp<-NULL

  for(k in 1:K){

    tmp<-cbind(tmp,dgamma(x,shape=alpha[k],scale=beta[k]))

  }

  return(tmp%*%lambda)
  
}

x <- seq(min(waiting),max(waiting),by=0.01)


hist(waiting, main = "Tiempo entre las erupciones del Old Faithful", xlab = "Minutos", ylab = "", cex.main = 1.5, 
cex.lab = 1.5, cex.axis = 1.4,probability=TRUE,col="green",ylim=c(0,0.05),border="white")

lines(x,gmixture(x,c(alfa[1],alfa[2]),c(beta[1],beta[2]),c(peso[1],peso[2])),type="l",ylab="density",col="darkblue")


pgmixture<-function(x,alpha,beta,lambda) {

  lambda<-lambda/sum(lambda)

  K<-length(lambda)

  tmp<-NULL

  for(k in 1:K){

    tmp<-cbind(tmp,pgamma(x,shape=alpha[k],scale=beta[k]))

  }

  return(tmp%*%lambda)
  
}

x <- seq(min(waiting),max(waiting),by=0.01)

plot(ecdf(waiting),verticals=TRUE,do.points=FALSE,col.hor="red", col.vert="bisque",ylab="")
lines(x,pgmixture(x,c(alfa[1],alfa[2]),c(beta[1],beta[2]),c(peso[1],peso[2])),type="l",ylab="density",col="darkblue")

gedest<-pgmixture(x,c(alfa[1],alfa[2]),c(beta[1],beta[2]),c(peso[1],peso[2]))


### PRUEBAS FORMALES DE BONDAD DE AJUSTE

mix.gamma<- function(x,alpha,beta, pmix) {
  pmix[1]*pgamma(x,alpha[1],beta[1])+pmix[2]*pgamma(x,alpha[2],beta[2])
}

ks.test(waiting,mix.gamma,alpha=alfa,beta=beta,pmix=peso)

ad.test(waiting,mix.gamma,alpha=alfa,beta=beta,pmix=peso)

### Cuál es la bronca, si se veía que ajustaba muy bien???


### Mezclas con el paquete: mclust

densWaiting <- densityMclust(faithful$waiting)
summary(densWaiting, parameters = TRUE)

plot(densWaiting)

plot(densWaiting, data = faithful$waiting)

1plot(densWaiting,what="diagnostic")

plot(densWaiting,what="BIC")

faithfulDens<-densityMclust(faithful)

plot(faithfulDens)

plot(faithfulDens, faithful, col = "blue", nlevels = 10, drawlabels = FALSE)

plot(faithfulDens, type = "image", col = topo.colors(50))

plot(faithfulDens, type = "persp", col = topo.colors(50))



#################################################################################################################

############################################ PRIMAS   ##########################################################

###Ejercicio (notas) de cálculo de primas
### Recuerde que cada reclamación es Bernoulli con qk prob. de reclamación y (1-qk) de no reclamación
### El modelo es el modelo individual: E(S)=sum(qj*E(Bj)); V(S)=qjV(Bj)+qj(1-qj)E(Bj)^2 (pags. 80-81 notas grandes)
### Con nuestros datos X=B: Montos B1=20, B1~U(0,20), B2=30, B2~U(0,30), B3=30, B3~U(0,50), B4=100, B4~U(0,100)

### Cálculo de E(N) y V(N)
### Cálculos para Nj~Bernoulli(qj) (Nj:# de estructuras)
### Este cálculo es un cálculo simple de una Bernoulli o Binomial para todos

q<-c(0.04,0.04,0.02,0.02)
est<-c(500,300,500,200) 

EN<-sum(est*q)

EN

VN<-sum(est*q*(1-q))
VN

###Cálculos E(S) y V(S)

#### Ahora hay que hacer los cálculos para S
### X~U(0,A) entonces E(X)=A/2; V(X)=A^2/12 -> Bk~U(0,Ak) -> E(Bk)=Ak/2 V(Bk)=Ak^2/12; k=1,2,3,4
### A1=20, A2=30, A3=50, A4=100
### S=S1+S2+S3+S4  con Sk total de reclamaciones para la k-ésima estructura

### El modelo es el modelo individual: E(S)=sum(qj*E(Bj)); V(S)=qjV(Bj)+qj(1-qj)E(Bj)^2

Ak<-c(20,30,50,100) ###Bj=Ak: j=1,2,3,4, k=1,2,3,4

ESk<-est*(Ak/2)*q
ESk

ES<-sum(ESk)
ES

### Var(Sk)=qkV(Sk)+qk(1-qk)E^2(Sk) (Modelo individual) ##En notas: qj*V(Bj)+qj*pj*E^2(Bj) y
### en nuestro caso, Bk~U(0,Ak)

VSk<-est*(q*Ak^2/12+q*(1-q)*(Ak/2)^2)
VSk

VS<-sum(VSk)
VS

### Valor de theta en la primera asignación del factor de recargo

### PI= (1+2*theta)[E(S1)+E(S2)]+(1+theta)[E(S3)+E(S4)]
###   =(2*theta)*380+(theta)*450+E(S1)+E(S2)+E(S3)+E(S4)
###   =1210*theta+ES
###  Entonces 0.99=P(S<1210*theta+ES) (Asumiendo una aproximación normal para S)
###  Entonces 0.99=P(S-ES<1210*theta)
###  Entonces 0.99=P((S-ES)/sd(S)<1210*theta/sd(S))
###  Entonces 0.99=P(Z<1210*theta/sd(S)) 

theta<-qnorm(0.99)*sqrt(VS)/1210  ###Utilizando aprox. Normal
theta

### Segunda asignación

### PII= (1+2*theta)*2*[E(S1) E(S2)]+(1+theta)[E(S3)+E(S4)]
###   =(1+2*theta)*2*380+(1+theta)*450
###   =1970*theta+ES
### Var(S)=2*Var(S1)+2*Var(S2)+Var(S3)+Var(S4)
###       =Var(S1)+Var(S2)+Var(S)
###  Entonces 0.99=P(S<1970*theta+ES)
###  Entonces 0.99=P(S-ES<-1970*theta)
###  Entonces 0.99=P((S-ES)/sd(S)<1970*theta/sd(S))


theta1<-qnorm(0.99)*sqrt(VS+sum(VSk[1:2]))/1970
theta1

### En términos generales, el factor de recargo disminuirá cuando el volumen del negocio (asegurados) aumenta si la frecuencia relativa y la gravedad (o tipo) de las reclamaciones
### sigue siendo el mismo (i.e. si las condiciones iniciales se mantienen intactas)

###Datos

base2<-read.csv("C:/Users/Salvador/Desktop/Cursos2020-II/SolvenciaII 2020/peragregas2.csv")

head(base2)

table(base2$frecuencia)

barplot(table(base2$frecuencia),col=rainbow(15))

c(mean(base2$frecuencia),var(base2$frecuencia))

poisfit2<-goodfit(base2$frecuencia, type = "poisson",method = "ML")
summary(poisfit2)
plot(poisfit2)
poisfit2$par


### Para la distribución de montos de reclamación

kk<-c(base2[,2],base2[,3],base2[,4],base2[,5],base2[,6],base2[,7],base2[,8],base2[,9],base2[,10],base2[,11],base2[,12],base2[,13],base2[,14])
kk1<-kk[!is.na(kk)]

plot(density(kk1),col="darkblue",lwd=2,main="Densidad subyacente",col.main="darkmagenta",xlab="")

###También puede hacerse 

vv<-as.vector(base2[,-1])
vv<-vv[ !is.na(vv)]

plot(density(vv),col="darkblue",lwd=2,main="Densidad subyacente",col.main="darkmagenta",xlab="")

### Posibles densidades asociadas (no todas)

descdist(kk1,discrete=FALSE,,boot=5000,obs.pch = 19, boot.col="darkviolet")

### Mi función (Recuerden, ustedes no la tienen)

source("C:/Users/Salvador/Desktop/Cursos2021-I/Riesgo2021/autodistfit.R")

fitData(kk1,c("gamma","weibull","lognormal","pareto","burr"),sample=1)

fB<-fitdist(kk1,"burr",method="mle",start=list(shape1=1,shape2=1,scale=1))

ks.test(kk1,"pburr",shape1=fB$estimate[1],shape2=fB$estimate[2],scale=fB$estimate[3])
ad.test(kk1, pburr,shape1=fB$estimate[1],shape2=fB$estimate[2],scale=fB$estimate[3])

plot(fB)

###Mismo argumento que antes: Valor de theta, factor de recargo, que nos de 99% de probabilidad
###que las primas exceden los reclamos
###Entonces, S tiene distribución Poisson compuesta con N~Poisson(5.225) y X~Burr(2.836407,1.056454,259.309589)

lambda<-5.225

ES<-lambda*mburr(1,shape1=fB$estimate[1],shape2=fB$estimate[2],scale=fB$estimate[3])
VS<-lambda*mburr(2,shape1=fB$estimate[1],shape2=fB$estimate[2],scale=fB$estimate[3])

ES;VS

###0.99=(1+theta)*ES
###0.99=P(S<ES+theta*ES)
###0.99=P(S-ES<theta*ES)
###0.99=P((S-ES)/sd(S)<theta*ES/sd(S))
###0.99=P(Z<theta*ES/sd(S))

theta<-qnorm(0.99)*sqrt(VS)/ES
theta


###CREDIBILIDA' (Como diría "El peje")

####################################################################################################################################################################################
########################################################################### CREDIBILIDAD ###########################################################################################
### Ejemplos de las notas
### Ejemplo 2: Tomamos la "marca de clase" como el valor representativo de todos los individuos en el intervalo

mr<-seq(200,3400,400)
mr 
N<-c(2,24,32,21,10,6,3,1,1)

cbind(mr,N)

EN<-mean(N)
VN<-var(N)
EN;VN

EX<-sum(mr*N)/(sum(N))
EX

EX2<-sum(N*mr^2)/(sum(N))
EX2

ES<-EN*EX
ES

VS<-EN*(EX2-EX^2)+VN*EX^2
VS

z<-qnorm(0.95)
k<-0.05

m<-z^2*VS/(k^2*ES^2)
m

###Credibilidad Bayesiana
### FRAUDE
###DistribuciónBetaInicial

alpha<-3.648 ; beta<-91.2

curve(dbeta(x,alpha,beta),col="darkorange",from=0.01,to=0.25,main="Distribución inicial: Fraude",col.main="darkblue",ylab=expression(p(theta)),xlab=expression(theta),col.lab="darkred",lwd=2)

y<-22
n<-240
a<-3.648
b<-91.2
x<-seq(0,1,length=5000)

beta.prior<-dbeta(x,a,b)

curve(dbeta(x,a,b),col="darkorange",from=0.01,to=0.25,main="Distribución inicial: Fraude",col.main="darkblue",ylab=expression(p(theta)),xlab=expression(theta),col.lab="darkred",lwd=2)

beta.lik<-dbeta(x,y+1,n-y+1)

curve(beta.lik(x,y+1,n-y+1),col="darkorange",from=0.01,to=0.25,main="Verosimilitud: Fraude",col.main="darkblue",ylab=expression(p(x|theta)),xlab=expression(theta),col.lab="darkred",lwd=2)

beta.posterior<-dbeta(x, y+a, n-y+b)

plot(x,beta.prior,type="l",main="Comparación: Inicial, Final, Verosimilitud",col="darkred",xlab=expression(theta),ylab="",col.main="darkorange",col.lab="darkblue",ylim=c(0,30),lwd=2,xlim=c(0,0.2))
lines(x,beta.posterior,col="darkgreen",lwd=2)
lines(x,beta.lik,col="darkblue",lwd=2,lty=2)
legend("topright",c("Distribución Inicial","Distribución Final","Verosimilitud"),col=c("darkred","darkgreen","darkblue"), pch=15, bty="n")


###Credibilidad con datos

D<-read.csv("C:/Users/Salvador/Desktop/Cursos2020-II/Riesgo 2020-II/Credibilidad.csv")

attach(D)

D1<-c(D[,2],D[,3],D[,4],D[,5],D[,6],D[,7],D[,8],D[,9],D[,10],D[,11],D[,12],D[,13])
D2<-D1[!is.na(D1)]

plot(density(D2),col="darkblue",main="Densidad subyacente",xlab="x",col.main="darkred",col.lab="darkorange")

descdist(D2,discrete=FALSE,,boot=5000,obs.pch = 19, boot.col="darkviolet")

source("C:/Users/Salvador/Desktop/Cursos2020-II/Riesgo 2020-II/autodistfit.R")

fitData(D2,c("weibull","lognormal","pareto"),sample=1)

###Distribución Weibull ajusta MUY BIEN

fW<-fitdist(D2,"weibull")

ks.test(D2,"pweibull",shape=fW$estimate[1],scale=fW$estimate[2])
ad.test(D2, pweibull,shape=fW$estimate[1],scale=fW$estimate[2])

plot(fW)

###Parte discreta

barplot(table(Freq),col=rainbow(10),main="Tabla de frecuencias de reclamación",col.main="darkblue")

###Ajuste modelo discreto: Poisson

disc.fit<-goodfit(Freq,type = "poisson", method = "ML")
plot(disc.fit)

summary(disc.fit)
disc.fit$par

###Los parámetros de las distribuciones ajustadas
###Weibull(0.530762112327404,4038.1)
###Poisson(4.97)
### Número total de reclamaciones

N<-sum(Freq)
N

###Credibilidad completa (0.05,0.9)

k<-0.05
p<-0.9

lambda<-4.97
alpha<-0.53
beta<-4038.10

###S es Poisson compuesto

ES<-lambda*mweibull(1,shape=alpha,scale=beta)
VS<-lambda*mweibull(2,shape=alpha,scale=beta)

ES;VS

m=qnorm((1+p)/2)^2*VS/(k^2*ES^2)
m

###Chin..."no nos alcanza N para credibilida' completa)

###¿Pa' qué sí alcanza?...credibilida' parcial

Z<-k*ES*sqrt(N)/(qnorm((1+p)/2)*sqrt(VS))
Z

###¿Cuánto nos falta para alcanzar credibilida' completa?

m-N

























