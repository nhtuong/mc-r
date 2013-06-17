#########################################
# classifySubjects.R - Classify subjects
# Author: Hoai Tuong Nguyen
# Created: 30/05/2013
# Modified: 17/06/2013
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
  output.dir="../3.results/CLASSIFICATION/Subjects/2013-06-17"

dir.create(output.dir, showWarnings = FALSE)


#Load imputed data
#seqcount.data<-get(load(sprintf("%s/DonneeSequences_t0_global_counting_table.csv.gz.opt.RData",input.dir)))
clinic.data<-get(load(sprintf("%s/T0ClinicalParametersMicrObesja_All.RData",input.dir)))

#Pre-processing data
if(F){
  clinic.data.1<-get(load(sprintf("%s/DonneeCliniques0_6_12Sem_3avril2010.xls.imputed.RData",input.dir)))
  tmp<-clinic.data[sapply(1:49, function(x) which(clinic.data$Names == clinic.data.1$NOM[x])),]
  heigh_clin_m<-clinic.data.1$Taille[1:49]
  names(heigh_clin_m)
  clinic.data<-cbind(tmp[,1:5],heigh_clin_m,tmp[,6:ncol(tmp)])
  save(clinic.data,file=sprintf("%s/T0ClinicalParametersMicrObesja_All.RData",input.dir))
  write.table(clinic.data,sprintf("%s/T0ClinicalParametersMicrObesja_All.csv",input.dir),row.names=F,sep=";")
}


#Annotation
#"diam_clin" = adipocyte diametre
#"tour.de.taille_clin"= Waist
#"MG.kg_clin" = Fat mas en KG (absolu)
#"MG.percent_clin" = Fat mass/body weight %
#"Taille"= Height


########################################
# PLOT CORRELATION AND ADD LOWESS LINE #
########################################

attach(clinic.data)

#t0 - "Adipocyte diameter" vs "Fat mass" by LOWESS
reg.plot.mc(fatmass_clin_kg,diam_clin_cm,
            type="lowess",
            pch=ifelse(Sexe=="M", 0, 1),         
            subjects=as.vector(Names),
            title="CORRELATION (T0) - LOWESS",
            xlab="Fat mass",ylab="Adipocyte diameter",
            legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass.pdf",output.dir),
            pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))

#t0 - "Waist" vs "Fat mass/Height*Height" by LOWESS
x=fatmass_clin_kg/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,waist.circumference_clin_cm,
            type="lowess",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - LOWESS",
            xlab="Fat mass/Height*Height",ylab="Waist",
            legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))



#t0 - "Adipocyte diameter" vs "Waist/Height*Height" by LOWESS
x=waist.circumference_clin_cm/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,diam_clin_cm,
            type="lowess",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - LOWESS",
            xlab="Waist/Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))


#t0 - "Adipocyte diameter" vs "Fat mass/Height*Height" by LOWESS
x=fatmass_clin_kg/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,diam_clin_cm,
            type="lowess",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - LOWESS",
            xlab="Fat mass/Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))






########################################
# PLOT CORRELATION AND ADD LINEAR LINE #
########################################


#t0 - "Adipocyte diameter" vs "Fat mass" by Linear Regression
reg.plot.mc(fatmass_clin_kg,diam_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),         
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="Fat mass",ylab="Adipocyte diameter",
            legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))

#t0 - "Waist" vs "Fat mass/Height*Height" by Linear Regression
x=fatmass_clin_kg/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,waist.circumference_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="Fat mass/Height*Height",ylab="Waist",
            legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))



#t0 - "Adipocyte diameter" vs "Waist/Height*Height" by Linear Regression
x=waist.circumference_clin_cm/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,diam_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="Waist/Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))




#t0 - "Adipocyte diameter" vs "Fat mass/Height*Height" by Linear Regression
x=fatmass_clin_kg/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,diam_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="Fat mass/Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))






########################################################
# BOXPLOT FOR CLASSIFIED GROUPS WITH LINEAR REGRESSION #
########################################################

lm.classified.wfmh2.t0<-read.csv(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
lm.classified.adfm.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
lm.classified.adwh2.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
lm.classified.adfmh2.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))


pdf(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lm.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lm.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lm.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lm.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()











pdf(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lm.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lm.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lm.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lm.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()












pdf(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lm.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lm.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lm.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lm.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



















########################################################
# BOXPLOT FOR CLASSIFIED GROUPS WITH LOWESS REGRESSION #
########################################################


lw.classified.wfmh2.t0<-read.csv(sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
lw.classified.adfm.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
lw.classified.adwh2.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
lw.classified.adfmh2.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))


pdf(sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lw.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lw.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lw.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lw.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()











pdf(sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lw.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lw.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lw.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lw.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()












pdf(sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lw.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lw.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lw.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lw.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()











