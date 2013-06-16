#########################################
# classifySubjects.R - Classify subjects
# Author: Hoai Tuong Nguyen
# Created: 30/05/2013
# Modified: 14/06/2013
# CMD: R --no-save --no-restore --slave -f classifySubjects.R "--args input.dir='?' output.dir='?'"
#########################################

#First read and parse the arguments for on the fly run
args=(commandArgs(TRUE))

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}


library(mc)

#default working path for inner run
#setwd("MetaCARDIS/")
#setwd("C:/Users/TUONG/Documents/metacardis/7.scripts/")

#default input and output paths
if (!exists("input.dir"))
  input.dir="../1.data/metacardis"
if (!exists("output.dir"))
  output.dir="../3.results/CLASSIFICATION/Subjects/2013-06-14"

dir.create(output.dir, showWarnings = FALSE)

#Load imputed data
clinic.data<-get(load(sprintf("%s/DonneeCliniques0_6_12Sem_3avril2010.xls.imputed.RData",input.dir)))
#seqcount.data<-get(load(sprintf("%s/DonneeSequences_t0_global_counting_table.csv.gz.opt.RData",input.dir)))


#Annotation
#"diam_clin" = adipocyte diametre
#"tour.de.taille_clin"= Waist
#"MG.kg_clin" = Fat mas en KG (absolu)
#"MG.percent_clin" = Fat mass/body weight %
#"Taille"= Height


########################################
# PLOT CORRELATION AND ADD LOWESS LINE #
########################################

#t0 - "Adipocyte diameter" vs "Fat mass" by LOWESS
reg.plot.mc(clinic.data$MG.kg_clin[1:49],clinic.data$diam_clin[1:49],
         type="lowess",
         pch=ifelse(clinic.data$Sexe[1:49]=="M", 0, 1),         
         subjects=as.vector(clinic.data$NOM[1:49]),
         title="CORRELATION (T0) - LOWESS",
         xlab="Fat mass",ylab="Adipocyte diameter",
         legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
         imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass.pdf",output.dir),
         pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))

#t0 - "Waist" vs "Fat mass/Height*Height" by LOWESS
x=clinic.data$MG.kg_clin[1:49]/(clinic.data$Taille[1:49]*clinic.data$Taille[1:49])
reg.plot.mc(x,clinic.data$tour.de.taille_clin[1:49],
         type="lowess",
         pch=ifelse(clinic.data$Sexe[1:49]=="M", 0, 1),
         subjects=as.vector(clinic.data$NOM[1:49]),
         title="CORRELATION (T0) - LOWESS",
         xlab="Fat mass/Height*Height",ylab="Waist",
         legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
         imgfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2.pdf",output.dir),
         pointsfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))



#t0 - "Adipocyte diameter" vs "Waist/Height*Height" by LOWESS
reg.plot.mc(clinic.data$tour.de.taille_clin[1:49]/(clinic.data$Taille[1:49]*clinic.data$Taille[1:49]),clinic.data$diam_clin[1:49],
         type="lowess",
         pch=ifelse(clinic.data$Sexe[1:49]=="M", 0, 1),
         subjects=as.vector(clinic.data$NOM[1:49]),
         title="CORRELATION (T0) - LOWESS",
         xlab="Waist/Height*Height",ylab="Adipocyte diameter",
         legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
         imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2.pdf",output.dir),
         pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))


#t0 - "Adipocyte diameter" vs "Fat mass/Height*Height" by LOWESS
x=clinic.data$MG.kg_clin[1:49]/(clinic.data$Taille[1:49]*clinic.data$Taille[1:49])
reg.plot.mc(x,clinic.data$diam_clin[1:49],
         type="lowess",
         pch=ifelse(clinic.data$Sexe[1:49]=="M", 0, 1),
         subjects=as.vector(clinic.data$NOM[1:49]),
         title="CORRELATION (T0) - LOWESS",
         xlab="Fat mass/Height*Height",ylab="Adipocyte diameter",
         legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
         imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2.pdf",output.dir),
         pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))






########################################
# PLOT CORRELATION AND ADD LINEAR LINE #
########################################


#t0 - "Adipocyte diameter" vs "Fat mass" by Linear Regression
reg.plot.mc(clinic.data$MG.kg_clin[1:49],clinic.data$diam_clin[1:49],
         type="lm",
         pch=ifelse(clinic.data$Sexe[1:49]=="M", 0, 1),
         subjects=as.vector(clinic.data$NOM[1:49]),
         title="CORRELATION (T0) - Linear Regression",
         xlab="Fat mass",ylab="Adipocyte diameter",
         legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
         imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass.pdf",output.dir),
         pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))

#t0 - "Waist" vs "Fat mass/Height*Height" by Linear Regression
reg.plot.mc(clinic.data$MG.kg_clin[1:49]/(clinic.data$Taille[1:49]*clinic.data$Taille[1:49]),clinic.data$tour.de.taille_clin[1:49],
         type="lm",
         pch=ifelse(clinic.data$Sexe[1:49]=="M", 0, 1),
         subjects=as.vector(clinic.data$NOM[1:49]),
         title="CORRELATION (T0) - Linear Regression",
         xlab="Fat mass/Height*Height",ylab="Waist",
         legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
         imgfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2.pdf",output.dir),
         pointsfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))



#t0 - "Adipocyte diameter" vs "Waist/Height*Height" by Linear Regression
reg.plot.mc(clinic.data$tour.de.taille_clin[1:49]/(clinic.data$Taille[1:49]*clinic.data$Taille[1:49]),clinic.data$diam_clin[1:49],
         type="lm",
         pch=ifelse(clinic.data$Sexe[1:49]=="M", 0, 1),
         subjects=as.vector(clinic.data$NOM[1:49]),
         title="CORRELATION (T0) - Linear Regression",
         xlab="Waist/Height*Height",ylab="Adipocyte diameter",
         legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
         imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2.pdf",output.dir),
         pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))




#t0 - "Adipocyte diameter" vs "Fat mass/Height*Height" by Linear Regression
x=clinic.data$MG.kg_clin[1:49]/(clinic.data$Taille[1:49]*clinic.data$Taille[1:49])
reg.plot.mc(x,clinic.data$diam_clin[1:49],
            type="lm",
            pch=ifelse(clinic.data$Sexe[1:49]=="M", 0, 1),
            subjects=as.vector(clinic.data$NOM[1:49]),
            title="CORRELATION (T0) - Linear Regression",
            xlab="Fat mass/Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))






#################################
# BOXPLOT FOR CLASSIFIED GROUPS #
#################################

classified.wfmh2.t0<-read.csv(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))

pdf(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_boxplot.pdf",output.dir))
outfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(4:6,10:98),function(x) boxplot.class.mc(data=clinic.data[1:49,],x,class=classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()


classified.adfm.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))

pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(4:6,10:98),function(x) boxplot.class.mc(data=clinic.data[1:49,],x,class=classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()


classified.adwh2.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))

pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(4:6,10:98),function(x) boxplot.class.mc(data=clinic.data[1:49,],x,class=classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()


classified.adfmh2.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))

pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(4:6,10:98),function(x) boxplot.class.mc(data=clinic.data[1:49,],x,class=classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()






alim.data<-get(load(sprintf("%s/DonneeCliniques0_6_12Sem_3avril2010_Alimentaire.xls.RData",input.dir)))

classified.wfmh2.t0<-read.csv(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
alim.classified.wfmh2.t0<-classified.adwh2.t0[which(classified.adwh2.t0[,1] %in% alim.data$NOM[1:24]),]
alim.classified.wfmh2.t0<-alim.classified.wfmh2.t0[order(alim.classified.wfmh2.t0$Name),]

pdf(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_Alim_boxplot.pdf",output.dir))
outfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_Alim_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(4,8:44),function(x) boxplot.class(alim.data[1:24,],x,alim.classified.wfmh2.t0[,4],outfile))
dev.off()


classified.adfm.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
alim.classified.adfm.t0<-classified.adfm.t0[which(classified.adfm.t0[,1] %in% alim.data$NOM[1:24]),]
alim.classified.adfm.t0<-alim.classified.adfm.t0[order(alim.classified.adfm.t0$Name),]

pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Alim_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Alim_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(4,8:44),function(x) boxplot.class(alim.data[1:24,],x,alim.classified.adfm.t0[,4],outfile))
dev.off()



classified.adwh2.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
alim.classified.adwh2.t0<-classified.adwh2.t0[which(classified.adwh2.t0[,1] %in% alim.data$NOM[1:24]),]
alim.classified.adwh2.t0<-alim.classified.adwh2.t0[order(alim.classified.adwh2.t0$Name),]

pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_Alim_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_Alim_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(4,8:44),function(x) boxplot.class(alim.data[1:24,],x,alim.classified.adwh2.t0[,4],outfile))
dev.off()








#TUONG

classified.wfmh.t0<-read.csv(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height_points.csv",output.dir))
alim.classified.wfmh.t0<-classified.adwh.t0[which(classified.adwh.t0[,1] %in% alim.data$NOM[1:24]),]
alim.classified.wfmh.t0<-alim.classified.wfmh.t0[order(alim.classified.wfmh.t0$Name),]

pdf(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height_Alim_boxplot.pdf",output.dir))
outfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height_Alim_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(4,8:44),function(x) boxplot.class(alim.data[1:24,],x,alim.classified.wfmh.t0[,4],outfile))
dev.off()


classified.adwh.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height_points.csv",output.dir))
alim.classified.adwh.t0<-classified.adwh.t0[which(classified.adwh.t0[,1] %in% alim.data$NOM[1:24]),]
alim.classified.adwh.t0<-alim.classified.adwh.t0[order(alim.classified.adwh.t0$Name),]

pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height_Alim_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height_Alim_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(4,8:44),function(x) boxplot.class(alim.data[1:24,],x,alim.classified.adwh.t0[,4],outfile))
dev.off()


classified.ahw2h4.t0<-read.csv(sprintf("%s/lm_t0_Height_vs_Waist2Height_4pi_points.csv",output.dir))
alim.classified.ahw2h4.t0<-classified.ahw2h4.t0[which(classified.ahw2h4.t0[,1] %in% alim.data$NOM[1:24]),]
alim.classified.ahw2h4.t0<-alim.classified.ahw2h4.t0[order(alim.classified.ahw2h4.t0$Name),]

pdf(sprintf("%s/lm_t0_Height_vs_Waist2Height_4pi_Alim_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Height_vs_Waist2Height_4pi_Alim_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(4,8:44),function(x) boxplot.class(alim.data[1:24,],x,alim.classified.ahw2h4.t0[,4],outfile))
dev.off()


p.clinic.adfm<-read.table(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_t-test.csv",output.dir),sep=";")
p.clinic.wfmh2<-read.table(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_t-test.csv",output.dir),sep=";")
p.clinic.adwh2<-read.table(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_t-test.csv",output.dir),sep=";")
p.alim.adfm<-read.table(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Alim_t-test.csv",output.dir),sep=";")
p.alim.wfmh2<-read.table(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_Alim_t-test.csv",output.dir),sep=";")
p.alim.adwh2<-read.table(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_Alim_t-test.csv",output.dir),sep=";")

p.clinic.all<-rbind(c("clinic.factors","ad.fm","w.fmh2","ad.wh2"),
                    cbind(as.character(p.clinic.adfm[,1]),
                          p.clinic.adfm[,3],
                          p.clinic.wfmh2[,3],
                          p.clinic.adwh2[,3]))

write.table(p.clinic.all,sprintf("%s/diff_high-low-risk_clinic_p-values.csv",output.dir),col.names=F,row.names=F,quote=F,sep=";")


p.alim.all<-rbind(c("alim.factors","ad.fm","w.fmh2","ad.wh2"),
                    cbind(as.character(p.alim.adfm[,1]),
                          p.alim.adfm[,3],
                          p.alim.wfmh2[,3],
                          p.alim.adwh2[,3]))

write.table(p.alim.all,sprintf("%s/diff_high-low-risk_alim_p-values.csv",output.dir),col.names=F,row.names=F,quote=F,sep=";")

p.star.clinic.all<-rbind(c("clinic.factors","ad.fm","w.fmh2","ad.wh2"),
                    cbind(as.character(p.clinic.adfm[,1]),
                          as.character(p.clinic.adfm[,4]),
                          as.character(p.clinic.wfmh2[,4]),
                          as.character(p.clinic.adwh2[,4])))




write.table(p.star.clinic.all,sprintf("%s/diff_high-low-risk_clinic_p-values-star.csv",output.dir),col.names=F,row.names=F,quote=F,sep=";")

p.star.alim.all<-rbind(c("alim.factors","ad.fm","w.fmh2","ad.wh2"),
                  cbind(as.character(p.alim.adfm[,1]),
                        as.character(p.alim.adfm[,4]),
                        as.character(p.alim.wfmh2[,4]),
                        as.character(p.alim.adwh2[,4])))


write.table(p.star.alim.all,sprintf("%s/diff_high-low-risk_alim_p-values-star.csv",output.dir),col.names=F,row.names=F,quote=F,sep=";")




diff.alim.groups.star<-read.table(sprintf("%s/diff_high-low-risk_alim_p-values-star.csv",output.dir),header=T,sep=";")
library(xtable)
print(xtable(diff.alim.groups.star),include.rownames=FALSE)


diff.clinic.groups.star<-read.csv(sprintf("%s/diff_high-low-risk_clinic_p-values-star.csv",output.dir),header=T,sep=";")
print(xtable(diff.clinic.groups.star),include.rownames=FALSE)


clinic.data.rate<-clinic.data[,10:ncol(clinic.data)]
clinic.data.rate<-sapply(10:ncol(clinic.data), function (x) clinic.data[,x]/max(clinic.data[,x]))
colnames(clinic.data.rate)<-names(clinic.data)[10:ncol(clinic.data)]
colnames(clinic.data.rate)

plot(clinic.data.rate[,1],clinic.data.rate[,4])
cor(clinic.data.rate[,1],clinic.data.rate[,4])

plot(clinic.data[,10],clinic.data[,13])
cor(clinic.data[,10],clinic.data[,13])

plot.mc<-function(x,y){
  x.rate<-x/max(x)
  y.rate<-y/max(y)
  plot(x.rate,y.rate,xlim=c(0,1),ylim=c(0,1))
}

plot.mc(clinic.data[,10],clinic.data[,13])

#'@name normality.mc
#'@aliases normality.mc
#'@docType methods
#'@rdname normality.mc
#'@title Normality Test
#'@description Perform a normality test for variables 
#'@param m variable (list, matrix, data frame...)
#'@param alpha p-value threshold
#'@return logical value indicating whether variable is normally distributed
#'@author Hoai Tuong Nguyen
#'@example
#' attach(mtcars)
#' normality.mc(mtcars)
normality.mc<-function(m,alpha=0.05){
  if (class(m)=="data.frame" || class(m)=="matrix")
    return(sapply(1:ncol(m), function(x) normality.mc(m[,x])$p.value<=alpha))
  else return(shapiro.test(m)$p.value<=alpha)
}

normality.mc(clinic.data[,10])

#'@name class.mc
#'@aliases class.mc
#'@docType methods
#'@rdname class.mc
#'@title Object Classes
#'@description Get class of variable
#'@param variable (list, matrix, data frame...)
#'@return class of variables or of columns of matrix/data frame
#'@author Hoai Tuong Nguyen
#'@example
#' attach(mtcars)
#' class.mc(mtcars)
class.mc<-function(m){
  if (class(m)=="data.frame" || class(m)=="matrix")
    return(sapply(1:ncol(m),function(x) class(m[,x])))
  else return(class(m))
}



data<-clinic.data[,c(4:6,10:55,97,98)]

cor.mat<-cor(data,data)
cor.mat<-cor.mat*cor.mat
diag(cor.mat)<-0
head(cor.mat[(cor.mat[,1]),1],20)

adjac.mat<-cor.mat>0.2
sapply(1:ncol(adjac.mat), function(x) which(adjac.mat[,x])

