date<-read.table(file="p.txt",row.names=1,sep="\t",header=TRUE,dec=".")
date
summary(date)

##1. Identificați regiunile care obțin scoruri peste Q75% pentru prima componentă principală. 

X_f<-scale(date,center=T,scale=T)
X_f
library(corrplot)
R<-cor(date)
corrplot(R,method="number",type="upper")
R_f<-cor(X_f)
corrplot(R_f,method="number",type="upper")
C<-cov(X_f)
C
library(moments)
library(FactoMineR)
library(factoextra)
library(nFactors)
library(dplyr)
acp<-PCA(X_f,scale.unit=TRUE)
acp
summary(acp)
fviz_eig(acp)
fviz_pca_var(acp,col.var="contrib")
acp2<-princomp(X_f,cor=TRUE,scores=TRUE)
acp2$scores
scoruri<-acp2$scores
q<-quantile(scoruri[,1],probs=c(0.25,0.5,0.75))
q
r<-which(acp2$scores[,1]>1.5165200)
r
View(acp2$scores)

##2. Care sunt variabilele X cu care această componentă este puternic corelată.
MF<-cor(X_f,scoruri)
corrplot(MF,method="number")

##3. Reprezentați cercul corelațiilor și extrageți două concluzii. 
fviz_pca_var(acp,col.var="contrib")

##4. Alegeți 2 variabile X corelate negativ. Transformați-le in variabile categoriale cu 4 nivele și aplicați 
##analiza corespondențelor. Identificați categoriile cu cele mai mari contribuții la inerția totală.
q1<-quantile(date$X6,probs=c(0.25,0.5,0.75))
q1
dim(date)
X6_cat<-rep(0,224)
X6_cat[date$X6<=q1[1]]<-"L1"
X6_cat[date$X6>q1[1]& date$X6<=q1[2]]<-"ML1"
X6_cat[date$X6>q1[2]& date$X6<=q1[3]]<-"MH1"
X6_cat[date$X6>q1[3]]<-"H1"
table(X6_cat)
q2<-quantile(date$X10,probs=c(0.25,0.5,0.75))
q2
X10_cat<-rep(0,224)
X10_cat[date$X10<=q2[1]]<-"L2"
X10_cat[date$X10>q2[1]& date$X10<=q2[2]]<-"ML2"
X10_cat[date$X10>q2[2]& date$X10<=q2[3]]<-"MH2"
X10_cat[date$X10>q2[3]]<-"H2"
table(X10_cat)
contingenta<-table(X6_cat,X10_cat)
contingenta
library(ca)
ac<-ca(contingenta)
summary(ac)
plot(ac)

##5. Efectuați analiza factorială pentru 50 de regiuni extrase aleator dintre cele din data frame-ul X. 
dim(date)
set.seed(100) 
indici<-sample(1:154,50,replace=FALSE)
date2<-date[indici,]
dim(date2)
library(psych)
d2_f<-scale(date2,center=TRUE,scale=TRUE)
R2<-cor(d2_f)
KMO(R2)
cortest.bartlett(R2)
af1<-factanal(d2_f,factors=2,scores="regression",rotation="none")
af1
summary(af1)
af1$uniqueness
af<-fa(d2_f,nfactors=2,n.obs=50,rotate="none",fm="ml")
fa.diagram(af)
af2<-fa(d2_f,nfactors=2,n.obs=50,rotate="none",fm="pa")
fa.diagram(af2)
summary(af2)
fsc<- factor.scores(d2_f, af)     
scoruri <- fsc$scores
scoruri[1:10,]
af3<-fa(d2_f,nfactors=2,n.obs=50,rotate="varimax",fm="ml")
fa.diagram(af3)
af4<-fa(d2_f,nfactors=2,n.obs=50,rotate="varimax",fm="pa")
fa.diagram(af4)
paralel<-fa.parallel(d2_f)
paralel$nfact
