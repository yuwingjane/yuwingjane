#Linear mixed model, LMM
library(lme4)#to build a LMM model
library(car)#to use Anova function
library(MuMIn)# to use r.squaredGLMM function
model2<-lmer(estimatedOC_log~biogeoCenterDepth.y_log+clayTotal_log+
               silt_plus_clay_log+Ca_plus_Mg_log+phH2o_log+
               alOxalate_log+feOxalate_log+mnOxalate_log+
               pOxalate_log+fedfeo_log+MAP_mm+MAT_C+
               (1+biogeoCenterDepth.y|plotID)+(1|siteID)+(1|domainID),
             data=m8rnsr,REML = FALSE)#ML was used for comparing models with
            #different fixed effects. REML often estimates the random effects better
summary(model2)
Anova(model2)#please note that anova() included in lmer package
#is not equal to Anova in car package. The former is used to
#compare model performances among models, while the latter is 
#used to get F-value and p-value
r.squaredGLMM(model2)
#Marginal R_GLMM2 represents the variance explained by fixed factors
#Conditional R_GLMM2 is interpreted as variance explained by both fixed and random factors
AIC(model2)# AIC is a proxy for model fit
vif(model2)#vif<3(or<10) or correlation coefficient between -0.7 and 0.7
#suggested low colinearity among variates
library(ggResidpanel)
resid_panel(model2,plots=c("qq","resid","yvp"))

#generalized additive mixed model, GAMM
library(mgcv)
model2<-gamm(estimatedOC_log~s(biogeoCenterDepth.y,bs="cr")+s(clayTotal,bs="cr")
             +s(silt_plus_clay,bs="cr")+s(Ca_plus_Mg,bs="cr")+s(phH2o,bs="cr")+
               s(alOxalate,bs="cr")+s(feOxalate,bs="cr")+s(mnOxalate,bs="cr")+
               s(pOxalate,bs="cr")+s(fedfeo,bs="cr")+s(MAT_C,bs="cr")+s(MAPminusPET,bs="cr"),random =list(plotID=~1,siteID=~1,domainID=~1),data=m8rnsr,method="ML")
summary(model2$gam)
gam.check(model2$gam)#Produce some diagnostic information about the fitting procedure and results
plot(model2$gam,pages=1,resid=T,pch=16,se=T,shade=T)#Plot the component smooth functions that make it up, on the scale of the linear predictor

#random forest model (RFM)
require(randomForest)
require(VSURF)
require(rfUtilities)

#Get the dataset
v14 <- m8rnsr[,c(1,3,18,19,7,10:13,17,45,48,31)]
r<-v14[complete.cases(v14),]
X=r[,-13]#Set up dataframe of predictors (i.e. environmental drivers)
Y = r[,13]#Define response variable 

#Run the RFM 
set.seed(42)#to ensure the outcomes of different iterations are repeatable
rf = randomForest(y = Y,x = X,keep.inbag=TRUE,importance=TRUE,ntree=10000)
rf.regression.fit(rf)#Evaluate fit & overfit of rfm

#Determine order of variable importance based on %IncMSE (increase of the mean squared error) value
windows(8,6)
varImpPlot(rf) 
importance(rf)

#Look at Partial Dependency Plots which visualize the marginal effects of features on the predicted outcome 
par(mfrow=c(6,3),mar=c(4,2.2,0,0))
#windows(4,4)
for(i in 1: ncol(X)){pplot<-partialPlot(rf,X,x.var = names(X)[i],xlab=names(X)[i],main="")}

#Look at predictor interaction
library(iml)#faster than 'pdp' package
library(processx)
set.seed(42)
predictor=Predictor$new(rf,data=X,y=Y)#Create a model object
interact=Interaction$new(predictor)#Measure the interaction strength, for each feature the interactions with all other features are estimated
interact$results
plot(interact)# Plot the resulting leaf nodes
interactMAT=Interaction$new(predictor,feature='MAT_C')#By selecting one feature name, the 2-way interactions of this feature with all other features are estimated
interactMAT$results
efffeoxdep=FeatureEffect$new(predictor,feature=c('feOxalate','biogeoCenterDepth.y'),method='pdp')# FeatureEffect plots support up to two features
efffeoxdep$results