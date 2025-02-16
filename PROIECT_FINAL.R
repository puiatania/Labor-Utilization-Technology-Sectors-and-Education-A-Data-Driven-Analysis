library(corrplot)
library(moments)
library(ggplot2)
date<-read.table(file="p.txt",sep="\t",header=TRUE,dec=".",row.names=1)
date
summary(date)

#SEMINARII 1-2
#Vectorul mediilor
medii <- apply(date,2,mean)
medii

# Vectorul coeficientilor de variatie
cv <- function (x) {
  c <- sd(x)/mean(x)
  return(c)}
cv(date[,1])##aplic coef de var pentru prima coloana

cvv <- apply(date,2,cv)##aplic pe tot setul de date
cvv

# Matrice de corelatie- Identificati corelatiile puternice
C <- cor(date)
corrplot(C,method="number",type="upper") 

## reprezentare grafica
library(corrplot)
corrplot(C,method="number",type="upper")

#SEMINAR 3
#Import în R pentru datele extrase de pe Eurostat. Observații la nivel regional (Xnxp). 
X<- read.table(file="p.txt",sep="\t",header=TRUE,dec=".",row.names=1) 
X
summary(X) 
dim(X) 
names(X)

#Extrageți vectorul coloană pentru care se înregistrează cea mai mare abatere standard. 
abateri_std <-apply(X,2,sd)
which.max(abateri_std)
Vc1<-X[,which.max(abateri_std)]

#Calculați cosinusul unghiului dintre vectorul extras la punctul 3 și vectorul dat de valorile 
#variabilei cu cea mai mică abatere standard. Centrați cei doi vectori și comparați cu valoarea 
#coeficientului de corelație dintre cele 2 variabile. 
Vc2 <-X[,which.min(abateri_std)]
#scale = standardizare

Vc1_centrat <-scale(Vc1, center=T,scale=F)
Vc2_centrat <-scale(Vc2, center=T,scale=F)

n1<-norm(Vc1_centrat,type="2")
n2<-norm(Vc2_centrat,type="2")

cos <-(t(Vc1_centrat)%*%Vc2_centrat)/(n1*n2)
cos

#Calculați matricea ce corelație pentru  cele p variabile. Interpretați. 
R<-cor(X)
round(R,3) #zecimale

corrplot(R, method="number", type='upper')
corrplot(R,  method="number")

#Calculați matricea produselor încrucișate pentru variabilele centrate. 
X_centrat <-as.matrix(scale(X,center=T, scale=F))
#matricea produselor incrucisate

prod_inc <-1/169*t(X_centrat)%*%X_centrat
round(prod_inc)
round(cov(X_centrat),2)
#obs ; matricea prod incucisate = matricea de corelatie pt elem centrate

#Calculați matricea de covarianța pentru variabilele standardizate.
X_std <-scale(X)
round(apply(X_std,2,mean),5)
round(apply(X_std,2,sd),5)

hist(X_std[,1])
round(cov(X_std),3)
round(cor(X_std),3)
round(cor(X),3)

#Extrageți valorile proprii ale matricii de covarianță.
S <- cov(X_std)
e <- eigen(S)
e
e$values
e$vectors

#SEMINAR 4
#Ce reprezinta matricea produselor incrucisate (pe date centrate)?
R11<-cor(X)
corrplot(R11,method="number",type="upper")

X_final<-X
R22<-cor(X_final)

corrplot(R22,method="number",type="upper")

X_centrat<-scale(X_final,center=T,scale=F)
n<-dim(X_centrat)[1]
n

#Matricea produsului incrucisat
Pi<-t(X_centrat)%*%X_centrat/n
C<-cov(X_centrat)
round(Pi-C,3)

#Descompuneți folosind valorile proprii și vectorii proprii
#RECAP: de ce valoriile proprii sunt reale și pozitive ?
E<-eigen(C) 
E
Lambda<-E$values # valorile proprii
Lambda
V<-E$vectors # matricea vectorilor proprii
V
#un vect propr inmultit cu o constanta este tot vector propriu

#Calculați norma primului vector propriu.  
#ptr vector 1
vector1<-V[,1]# vect prop
norm(vector1,type="2")#imp sa fie 1
sum(vector1^2)

#ptr vector 6
vector6<-V[,6]#vector propriu
norm(vector6,type="2")#imp sa fie 1
sum(vector6^2)

#OBS: Vectorii prop ai matricei de covarianta au norma 1

#Verificați descompunerea 
descompunereJ <- V%%diag(Lambda)%%t(V) 
round(Pi-descompunereJ,3) #Pi=C

#Folosiți funcția princomp() in R 
#comp care extrage comp princ/nu aplicam pe matricea de corelatie 
acp <- princomp(X_centrat,cor=F,scores=TRUE) 
summary(acp)
#coloanele loadings=vectorii proprii ai matricei de covarianta !!!!

#Care este legătura dintre valorile cu abaterile standard de mai sus și valorile Lambda? 
sqrt(Lambda) 
#varianta lui Z este Lambda 

#Ce sunt loadings? Comparati cu V? Ce observati? 
acp$loadings 
V 
#OBS:  Înmulțirea unui vector propriu cu o constantă nu schimbă calitatea de vector propriu al unei matrici.
#ce facem cu vect prop 
#cu ajutorul lor scriem forma principala a comp principale

#Ce reprezintă scores? Ce dimensiune are matricea scoruri? 
Coeficienti<-acp$loadings

scoruri <- acp$scores 
Ponderare <- X_centrat %*% Coeficienti
X_centrat[1:3,]#datele initiale centrate pe cele 3 regiuni
scoruri[1:3,]
Ponderare[1:3,]

dim(scoruri)
dim(X_centrat)
#obt aceeasi dimensiune

C_L<-cov(scoruri)
round(C_L,3)
Lambda
#Lambda mare=cov ptr Z
#0 in rest, pe diag principale val

#DPDV al comp principale=aceste val lambda sunt in ordine descrescatoare

#corelatie componente principale
cor_cp<-cor(scoruri)
round(cor_cp)

#5)Folosiți funcția princomp() in R 
acp <- princomp(X_centrat,cor=T,scores=TRUE) 
summary(acp) 

#6)Care este legătura dintre valorile cu abaterile standard de mai sus și valorile Lambda? 
sqrt(Lambda) 

#7)Ce sunt loadings? Comparati cu V? Ce observati? 
acp$loadings 
V 

#8)Ce reprezintă scores? Ce dimensiune are matricea scoruri? 
Coeficienti<-acp$loadings

scoruri <- acp$scores 
Ponderare <- X_centrat %*% Coeficienti
X_centrat[1:3,]
scoruri[1:3,]
Ponderare[1:3,]


dim(scoruri)
dim(X_centrat)

C_L<-cov(scoruri)
round(C_L,3)
Lambda

#corelatie componente principale
cor_cp<-cor(scoruri)
round(cor_cp)

#SEMINAR 5
acp <- princomp(X_final,cor=T,scores=T)
var <- acp$sdev
#ab std ale componentelor principale
coef <- acp$loadings
scoruri <- acp$scores

variante<-acp$sdev*acp$sdev
#varianta

summary(acp)

#Propr comp princ
sum(variante)

plot(acp,type="l")
#val lambda in raport cu fiec comp
#screeplot

print(acp)

R <- cor(X_final)
E <- eigen(R)
V <- E$vectors # vectorul coeficienților. Uneori sunt coeficienții “loadings” cu semn schimbat.
lambda <- E$values
sqrt(lambda) # abaterile standard ale componentelor principale

comblin <- scale(X_final,center=T,scale=T)%*%coef # sunt scorurile
coef

scale(X_final,center=T,scale=T)[1:3,]
comblin[1:3,]
acp$scores[1:3,]

View(comblin)

#Matricea factor, folosită pentru a da noilor caracteristici construite o interpretare concretă,
#evidențiază corelațiile dintre variabilele X și componentele principale Z.
MF=cor(X_final,scoruri)
corrplot(MF,method="number")

biplot(acp)

#TEMA:
#1. Instalați cele două biblioteci: 
library(factoextra) 
library(FactoMineR)

#2. Folosiți funcția PCA() pentru a realiza analiza componentelor principale. 
acp2 <- PCA(X_final) 

#3. Extrageți outputurile folosind următoarele funcții:
summary(acp2)
fviz_eig(acp2)
fviz_pca_ind(acp2) 
fviz_pca_var(acp2,col.var="contrib")

#SEMINAR 6
X_final2<-X_final
dim(X_final2)

acp2<-PCA(X_final2,scale.unit=T,graph=T)
acp2
summary(acp2)

#Pentru VARIABLES
#1
sc<-acp2$ind$coord
round(acp2$var$coord,3)
round(cor(X_final2,sc),3)

#2
acp2$var$cos2
acp2$var$coord*acp2$var$coord

#3
acp2$eig
round(acp2$var$cos2,3)
sum(acp2$var$cos2[,1])
sum(acp2$var$cos2[,2])
sum(acp2$var$cos2[,3])
sum(acp2$var$cos2[,4])

#4
acp2$var$contrib
acp2$var$cos2[,1]/sum(acp2$var$cos2[,1])*100
acp2$var$cos2[,1]/sum(acp2$eig[1,1])*100

#grafice
fviz_eig(acp2)
acp2$eig
fviz_pca_ind(acp2)# repr grafica a observatiilor in planul principal 
#punctele=coordonatele de la indiv
acp2$ind$coord[1:5,]

fviz_pca_var(acp2,col.var="contrib")#cercul corelațiilor
# cu cat sageata este mai mare, cu atat este mai mare corelatia
summary(acp2)

sc<-data.frame(acp2$ind$coord)
ggplot(sc,aes(x=sc[,1],y=sc[,2]))+
  geom_point(shape=16,alpha=.4)+
  geom_text(label=row.names(sc),vjust=0,hjust=0,size=5)+
  labs(x="Z1",y="Z2")+
  theme_minimal()+
  scale_color_gradient(low="#0091ff",high="#f0650e")+
  theme(axis.text.y=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        legend.position="botton")

#SEMINAR 9
library(psych)
R<-cor(X_final)
KMO(R)

cortest.bartlett(R)#ipoteze test Bartlett ???
#stat B=-(n-1-(2p+5)/6)*ln|R|~chisq cu p(p-1)/2 DF
#n=nr de regiuni
#matricea de corelatie=R
#p=nr de variabile
#p-value<0.5 => se respinge ipoteza H0, conform careia variabilele sunt ortogonale
#45 de grade de libertate

#MAXIMUM LIKELIHOOD=ML =>metoda verosimilitatii maxime
#sau PC,PA=principal axes 

##Utilizați funcția factanal () pentru a efectua analiza factorială. Soluția se bazează pe metoda 
##verosimilității maxime.  
af1<-factanal(X_final,factors=2,scores="regression",rotation="none")
af1
summary(af1)
af1$uniqueness#unicitati
af1$scores[1:5,]#scoruri
af1$loadings#intensitatile

##Utilizați funcția fa() din biblioteca psych. 
#metoda ML
dim(X_final)
af<-fa(R,nfactors=2,n.obs=154,rotate="none",fm="ml")
fa.diagram(af)# nu va afisa corelatii mai mici de 0.4

#metoda PA
af<-fa(R,nfactors=2,n.obs=154,rotate="none",fm="pa")
fa.diagram(af)
summary(af)
af
#col PA1,PA2 inseamna loadings
#col h2 este comunalitate=varianta explicata de factori
#val sunt suma patratelor coef loadings
# h2=PA1^2+PA2^2
# pentru primul este 1.164=(-0.75)^2+(0.78)^2
# cu cat este mai mare unicitatea, cu atat factorii surprind mai putin din varianta variabilei
# cand avem u2 negativ sau >1, datele nu sunt ok, nu sunt calculate bine =>caz ultra-Heywood
#coloana com este indicele de complexitate Hoffman
##SS loadings sunt valorile proprii, suma patratelor loadings
Ld<-af$loadings
ld1<-Ld[,1]
sum(ld1^2)#gal cu SS loadings
ld2<-Ld[,2]
sum(ld2^2)
#Proportion Var: Proporția preluată din varianța totală 
#Cumulative var: cumulat 
#Proportion explained: relativ =Proportion var/sum(proportion var) 
#0.60=0.37/(0.37+0.25)
#Cumulative proportion: relativ cumulat 

#scoruri factoriale
fsc<- factor.scores(X_final, af)     
scoruri <- fsc$scores
scoruri[1:10,]

#rotire varimax
af<-fa(R,nfactors=2,n.obs=154,rotate="none",fm="ml")
fa.diagram(af)
af<-fa(R,nfactors=2,n.obs=154,rotate="varimax",fm="ml")
fa.diagram(af)
#varimax vine de la variance maximization
#rotatia ortogonala a a axelor prin max variantei unei variabile
# sea ccentueaza semnificatia informationala a unor variabile in raport cu altele
af<-fa(R,nfactors=2,n.obs=154,rotate="varimax",fm="pa")
fa.diagram(af)

#numarul de factori
library(nFactors)
library(psych)
paralel<-fa.parallel(X_final)
paralel$nfact
#fiec punct de pe liniile albastre(actual data) care este deasupra liniilor rosii 
#este un factor care trebuie pastrat

##SEMINAR 13
X<-read.table(file = "p.txt",sep="\t",header=TRUE,dec=".",row.names=1)
dim(X)
names(X)
acp<-princomp(X,cor=T,scores=T)
#sau fctia PCA - $ind$coord
scoruri<-data.frame(acp$scores[,1:2])
names(scoruri)<-c("Z1","Z2")
scoruri[1:5,]

library(cluster)
###1.Matricea distantelor
d<-dist(scoruri)
d[1]
##distanta euclidiana dintre primele 2 regiuni, BE10 și BE21, 
##fiind radical din (-0.07+0.06)^2+(-1.49+0.31)^2

ierarhie<-hclust(d,method="single")
#"single","complete","average","centroid"
#"ward.D2"
#?hclust
plot(ierarhie)
##dendrograma pentru metoda celor mai apropiati vecini 
##Criterii de alegere a numarului de clase
##1) Criteriul dendogramei
##se realizeaza o taietura in dendograma a.i. dist dintre 2 pasi consecutivi sa fie cea mai mare
##nr de intersectii orizontale ale taieturii cu dendograma=nr de clase
##realizam o statistica descriptiva a datelor
##extragem cp, realizam metodele ierarhica ptr a determina numarul de clase
##algoritmul camils ptr a determ apart tarilor la clase
##ES61 pare a fi un outlier, height=dist de agregare, in cazul meu creste
##2) Criteriul -cum vrea autorul
##3) Criteriul NbCluster
##intoarce un rezultat (votul majoritatii) 
ierarhie<-hclust(d,method="ward.D2")
plot(ierarhie)

ierarhie$height##distantele de agregare
ierarhie$labels##etichetele regiunilor noastre

windows()
plot(ierarhie)
rect.hclust(ierarhie,k=5,border=2:5)
?rect.hclust()

solutie5<-cutree(ierarhie,k=5)##apartenenta la clase
table(solutie5)##ptr fiec clasa avem nr de obiecte care fac parte din ea
aggregate(scoruri,list(solutie5),mean)##putem analiza clasele si sa le determinam caracteristicile cu ajutorul centroizilor 

library(factoextra)
sol<-hcut(scoruri,k=5,hc_method="ward.D2")
#reprezentare solutie
windows()
fviz_cluster(list(data=scoruri,cluster=solutie5))
##pentru solutiile extrase cu cutree
fviz_cluster(sol)
#?fviz_cluster()

##repr grafic silhouette
?silhouette()
s<-silhouette(solutie5,d)
plot(s)
s[1:5,]
##daca si tinde catre 1, atunci observatia este bine incadrata in clasa
##daca si se apropie de 0 inseamna ca observatia respectiva este intre 2 clase
##daca avem si<0, atunci observatia este gresit incadrata 

##Algoritm kmeans
km<-kmeans(scoruri,5)
clase<-km$cluster
#apartenenta la clase
table(clase)
fviz_cluster(km,scoruri)

s<-silhouette(clase,d)
plot(s)

#?kmeans()
km$centers#centroizii
km$totss#variabilitatea totala
km$withinss#variabilitatea intra clasa
km$tot.withinss#variabilitatea intra clasa totala
km$betweenss#variabilitate inter clasa
