

# Final Project
# Statistical Connectomics
# Alan Juliano


install.packages('e1071',dependencies=TRUE) 
install.packages("MASS")
library(e1071)
library(MASS)
install.packages("class")
install.packages("DAAG")
library(class)
librar("DAAG")
data<-read.csv("~/statconnfinal.csv")
attach(data)
names(data)

#Using gender as a factor
data$gender<-factor(data$gender)

# ADHD and Depression Regions of Brain 
# Biologically, research has indicated that ADHD is a deficiency within the premotor cortex along with
# the superior prefrontal cortex. Depression is displayed through defficiencies in the frontal cortex, 
# hippocampus, amygdala, striatum and thalamus regions of the brain. 

# In our data, V8-V15 and V22-30 are the regions of interest when examining ADHD
# Depression effects activity in the V56-V60 (Hippocampus), V63-65(Amygdala), V32-35(Striatum), V78-82 (Thalamus)  

#ADHD Region of the Brain

# From here we are able to run logistic regression 

 mylogit<-glm(adhd~V8+V9+V10+V11+V12+V13+V14+V15+V22+V23+V24+V25+V26+V27+V28+V29+V30+gender,data=data,family="binomial")
 summary(mylogit)

#Call:
#glm(formula = adhd ~ V8 + V9 + V10 + V11 + V12 + V13 + V14 + 
#    V15 + V22 + V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30 + 
#    gender, family = "binomial", data = data)

#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.8080  -0.5263  -0.4744  -0.4216   2.4223  

#Coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -2.1026690  0.0999466 -21.038   <2e-16 ***
#V8          -0.0115285  0.0103788  -1.111   0.2667    
#V9          -0.0079060  0.0059359  -1.332   0.1829    
#V10         -0.0042993  0.0060658  -0.709   0.4785    
#V11         -0.0054432  0.0071121  -0.765   0.4441    
#V12          0.0016451  0.0076320   0.216   0.8293    
#V13         -0.0065502  0.0084538  -0.775   0.4384    
#V14         -0.0059843  0.0086006  -0.696   0.4865    
#V15          0.0169776  0.0080167   2.118   0.0342 *  
#V22          0.0041620  0.0051870   0.802   0.4223    
#V23          0.0033985  0.0066132   0.514   0.6073    
#V24         -0.0015768  0.0079091  -0.199   0.8420    
#V25          0.0046741  0.0041884   1.116   0.2644    
#V26         -0.0003988  0.0031019  -0.129   0.8977    
#V27          0.0057526  0.0057730   0.996   0.3190    
#V28         -0.0046816  0.0050992  -0.918   0.3586    
#V29          0.0035973  0.0074636   0.482   0.6298    
#V30         -0.0055721  0.0074746  -0.745   0.4560    
#gender1      0.0635755  0.1394597   0.456   0.6485    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#    Null deviance: 1463.6  on 2047  degrees of freedom
#Residual deviance: 1445.7  on 2029  degrees of freedom
#AIC: 1483.7

#Number of Fisher Scoring iterations: 5



#V15 is the most statistically significant region when identifying ADHD patients. 

#Depression Region
 mylogit2<-glm(depression~V56+V57+V58+V59+V60+V63+V64+V65+V32+V33+V34+V35+V26+V78+V79+V80+V81+V82+gender,data=data,family="binomial")
summary(mylogit2)
#Call:
#glm(formula = depression ~ V56 + V57 + V58 + V59 + V60 + V63 + 
#    V64 + V65 + V32 + V33 + V34 + V35 + V26 + V78 + V79 + V80 + 
#    V81 + V82 + gender, family = "binomial", data = data)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.0015  -0.6490  -0.5912  -0.5190   2.1303  

#Coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -1.637e+00  8.435e-02 -19.402   <2e-16 ***
#V56         -1.023e-02  7.642e-03  -1.338   0.1809    
#V57         -9.181e-03  9.203e-03  -0.998   0.3185    
#V58         -1.532e-03  1.040e-02  -0.147   0.8829    
#V59         -5.884e-03  7.474e-03  -0.787   0.4311    
#V60         -7.605e-03  7.093e-03  -1.072   0.2836    
#V63         -1.908e-04  6.096e-03  -0.031   0.9750    
#V64          1.885e-02  7.762e-03   2.428   0.0152 *  
#V65          7.952e-04  3.988e-03   0.199   0.8420    
#V32         -2.652e-05  5.039e-03  -0.005   0.9958    
#V33         -3.578e-03  7.563e-03  -0.473   0.6362    
#V34         -7.759e-03  7.998e-03  -0.970   0.3320    
#V35          2.687e-03  3.318e-03   0.810   0.4180    
#V26          6.395e-04  2.421e-03   0.264   0.7917    
#V78          1.008e-02  4.643e-03   2.172   0.0299 *  
#V79          5.251e-03  3.107e-03   1.690   0.0911 .  
#V80         -2.231e-03  3.635e-03  -0.614   0.5393    
#V81          1.826e-03  6.045e-03   0.302   0.7627    
#V82         -3.966e-03  8.093e-03  -0.490   0.6241    
#gender1      1.169e-01  1.177e-01   0.993   0.3207    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#    Null deviance: 1895.1  on 2047  degrees of freedom
#Residual deviance: 1873.4  on 2028  degrees of freedom
#AIC: 1913.4

#Number of Fisher Scoring iterations: 4


# In our logistic model, we see that V59 is the most statistically significant when identifying depression. 

#Removing gender variable from Depression
> mylogit3<-glm(depression~V56+V57+V58+V59+V60+V63+V64+V65+V32+V33+V34+V35+V26+V78+V79+V80+V81+V82,data=data,family="binomial")
> summary(mylogit3)

#Call:
#glm(formula = depression ~ V56 + V57 + V58 + V59 + V60 + V63 + 
#    V64 + V65 + V32 + V33 + V34 + V35 + V26 + V78 + V79 + V80 + 
#    V81 + V82, family = "binomial", data = data)

#Deviance Residuals: 
 #   Min       1Q   Median       3Q      Max  
#-1.0211  -0.6467  -0.5910  -0.5211   2.1516  

#Coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -1.578e+00  5.939e-02 -26.573   <2e-16 ***
#V56         -1.030e-02  7.636e-03  -1.348   0.1776    
#V57         -9.501e-03  9.196e-03  -1.033   0.3015    
#V58         -1.370e-03  1.040e-02  -0.132   0.8952    
#V59         -5.673e-03  7.469e-03  -0.760   0.4475    
#V60         -7.679e-03  7.091e-03  -1.083   0.2788    
#V63         -6.515e-07  6.088e-03   0.000   0.9999    
#V64          1.870e-02  7.757e-03   2.411   0.0159 *  
#V65          7.183e-04  3.980e-03   0.180   0.8568    
#V32          2.240e-04  5.026e-03   0.045   0.9645    
#V33         -3.576e-03  7.567e-03  -0.473   0.6365    
#V34         -7.908e-03  8.001e-03  -0.988   0.3229    
#V35          2.721e-03  3.311e-03   0.822   0.4113    
#V26          6.096e-04  2.420e-03   0.252   0.8011    
#V78          9.921e-03  4.638e-03   2.139   0.0324 *  
#V79          5.223e-03  3.109e-03   1.680   0.0930 .  
#V80         -2.155e-03  3.634e-03  -0.593   0.5531    
#V81          1.697e-03  6.039e-03   0.281   0.7787    
#V82         -3.988e-03  8.089e-03  -0.493   0.6220    
---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

 #   Null deviance: 1895.1  on 2047  degrees of freedom
#Residual deviance: 1874.4  on 2029  degrees of freedom
#AIC: 1912.4

Number of Fisher Scoring iterations: 4
# V64 and V68 are again the most statistically significant

#Removing gender variable from ADHD
> mylogit4<-glm(formula = adhd ~ V8 + V9 + V10 + V11 + V12 + V13 + V14 + 
+     V15 + V22 + V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30, family = "binomial", data = data)
> summary(mylogit4)
#
#Call:
#glm(formula = adhd ~ V8 + V9 + V10 + V11 + V12 + V13 + V14 + 
#    V15 + V22 + V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30, 
#    family = "binomial", data = data)

#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.7989  -0.5260  -0.4739  -0.4217   2.4143  

#Coefficients:
 #            Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -2.071095   0.071230 -29.076   <2e-16 ***
#V8          -0.011611   0.010371  -1.120    0.263    
#V9          -0.007883   0.005932  -1.329    0.184    
#V10         -0.004260   0.006065  -0.702    0.482    
#V11         -0.005543   0.007104  -0.780    0.435    
#V12          0.001734   0.007630   0.227    0.820    
#V13         -0.006515   0.008448  -0.771    0.441    
#V14         -0.006065   0.008599  -0.705    0.481    
#V15          0.016985   0.008013   2.120    0.034 *  
#V22          0.004175   0.005183   0.806    0.421    
#V23          0.003457   0.006606   0.523    0.601    
#V24         -0.001678   0.007906  -0.212    0.832    
#V25          0.004711   0.004189   1.125    0.261    
#V26         -0.000399   0.003100  -0.129    0.898    
#V27          0.005730   0.005772   0.993    0.321    
#V28         -0.004730   0.005098  -0.928    0.354    
#V29          0.003619   0.007465   0.485    0.628    
#V30         -0.005663   0.007471  -0.758    0.448    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

 #   Null deviance: 1463.6  on 2047  degrees of freedom
#Residual deviance: 1445.9  on 2030  degrees of freedom
#AIC: 1481.9

#Number of Fisher Scoring iterations: 5




################

#Depression 

mylogit2<-glm(depression~V56+V57+V58+V59+V60+V63+V64+V65+V32+V33+V34+V35+V26+V78+V79+V80+V81+V82+gender,data=data,family="binomial")
summary(mylogit2)

### Removing gender from the data
mylogit3<-glm(depression~V56+V57+V58+V59+V60+V63+V64+V65+V32+V33+V34+V35+V26+V78+V79+V80+V81+V82,data=data,family="binomial")
summary(mylogit3)
# Comparing the two models (with and without gender)
anova(mylogit3,mylogit2,test="Chisq")

# ADHD
mylogit<-glm(adhd~V8+V9+V10+V11+V12+V13+V14+V15+V22+V23+V24+V25+V26+V27+V28+V29+V30+gender,data=data,family="binomial")
 summary(mylogit)
 #removing gender
mylogit4<-glm(formula = adhd ~ V8 + V9 + V10 + V11 + V12 + V13 + V14 + 
V15 + V22 + V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30, family = "binomial", data = data)
summary(mylogit4)
anova(mylogit4,mylogit,test="Chisq")

# Refining Depression Interpretation
 mylogit5<-glm(depression~V56+V57+V59+V60+V64+V33+V34+V35+V78+V79+V80+V82,data=data,family="binomial")
summary(mylogit5)
# Displays further connectivity
mylogit6<-glm(depression~V56+V57+V59+V60+V64+V34+V35+V78+V79,data=data,family="binomial")
summary(mylogit6)
# Further exemplifies that the Thalamas plays an integral role in depression
anova(mylogit6,mylogit5,test="Chisq")

#Refining ADHD
mylogit7<-glm(formula = adhd ~ V8 + V9 + V10 + V11 + V13 + V14 + 
V15 + V22 + V27 + V28 + V30, family = "binomial", data = data)
summary(mylogit7)

# Comparison Between Statistically Significant Regions
adhdfinal<-data[,c("V8","V9","V10","V11","V13","V14","V15","V22","V27","V28","V30")]
depressionfinal<-data[,c("V56","V57","V59","V60","V64","V34","V35","V78","V79")]
correlation<-cor(adhdfinal,depressionfinal)
write.table(correlation)
 
# Plotting Statistically Significant Regions vs. Disease Presence
par(mfrow=c(1,3))
plot(depression,V64,xlab="Depression (1=Yes 0=No)",ylab="fMRI of V64") 
plot(depression,V78,xlab="Depression (1=Yes 0=No)",ylab="fMRI of V78") 
plot(adhd,V15,xlab="ADHD (1=Yes 0=No)",ylab="fMRI of V15") 
mtext("ADHD & Depression vs. Statistically Significant Regions of the Brain", side=3, outer=TRUE, line=-3)
