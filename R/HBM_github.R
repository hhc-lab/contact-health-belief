#############Packages#############
install.packages("dplyr")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("pscl") 
install.packages("boot")
install.packages("readxl")
install.packages("writexl")
install.packages("gee")
install.packages("geepack")
install.packages("texreg")
install.packages("plotly")
install.packages("socialmixr")
install.packages("conmat", repos = "https://idem-lab.r-universe.dev")
install.packages("tidyr")
install.packages("parameters")
install.packages("broom.helpers")
install.packages("Hmisc")
install.packages("car")
library(dplyr)
library(ggplot2)
library(gridExtra)
require(pscl)
require(boot)
library(readxl)
library(writexl)
library(gee)
library(geepack)
library(texreg)
library(plotly)
library(socialmixr)
library(conmat)
library(tidyr)
library(parameters)
library(broom.helpers)
library(Hmisc)
library(car)
library(tidyverse)
#############Load data#############
Longitudinal_r<- read_excel("Longitudinal_clean.xlsx")
Cross<- read_excel("Cross_clean.xlsx")
confirmedcase<- read_excel("Confirmation_2020_2023.xlsx")

#############Functions#############
description_table<- function(x){
  aa<- table(x)
  bb<- aa %>% prop.table()
  tables<- data.frame(name=row.names(aa),number= c(aa),percentage= c(bb))
  return (tables)
}

meanse<- function(x){
  tempmean<- mean(x, na.rm = TRUE)
  tempse<- sd(x, na.rm = TRUE)/sqrt(length(which(is.na(x)==FALSE)))
  output<- paste0(round(tempmean,2)," ± ",round(tempse,2))
  return(output)
}

num_per<- function(x){
  high_number<- length(which(x=="high"))
  high_percentage<- length(which(x=="high"))/length(x[is.na(x)==FALSE])
  low_number<- length(which(x=="alow"))
  low_percentage<- length(which(x=="alow"))/length(x[is.na(x)==FALSE])
  output<- paste0(high_number,"(",scales::percent(high_percentage,0.1),")", low_number,"(",scales::percent(low_percentage,0.1),")")
  return(output)
}

mcnemar<-function(x,y){
  temp<-matrix(c(length(na.omit((Longitudinal_r$ID[x=="high" & y=="high"]))), length(na.omit((Longitudinal_r$ID[x=="high" & y=="alow"]))),
                 length(na.omit((Longitudinal_r$ID[x=="alow" & y=="high"]))), length(na.omit((Longitudinal_r$ID[x=="alow" & y=="alow"])))),
               nrow = 2,dimnames = list("Baseline" = c("High", "Low"),"Outbreak" = c("High", "Low")))
  temp2<-temp/2
  out<-mcnemar.test(temp2, correct = TRUE) 
  return(out$p.value)
}

geetable<- function(x){
  temp<-tidy(x, conf.int = TRUE)
  OddRatio<-round(exp(temp$estimate),2)
  lwr<-round(exp(temp$conf.low),2)
  upr<-round(exp(temp$conf.high),2)
  res<-data.frame(Var = temp$term,
                  B = round(temp$estimate,2),
                  OR_CI = paste0(OddRatio," [",lwr,", ",upr,"]"),
                  p_value = round(temp$p.value,3))
  rownames(res) <- temp$term
  return(res[-1,])
}

chi2<- function(x,y){
  d<- as.table(rbind(table(x),table(y)))
  temp<- chisq.test(d)
  return(temp$p.value)
}

zip_table<- function(x){
  coeffi<-coef(x)
  temp<-x %>% tidy_zeroinfl(exponentiate = TRUE)
  OddRatio<- round(temp$estimate,2)
  lwr<-round(temp$conf.low,2)
  upr<-round(temp$conf.high,2)
  res<-data.frame(Var = temp$original_term,
                  B = round(coeffi,3),
                  OR_CI = paste0(OddRatio," [",lwr,", ",upr,"]"),
                  p_value=round(temp$p.value,4))
  return(res)
}

#############Analyze data#############

#########1-Longitudinal_Description#########
#######1-1 Belief Change#######
##Belief score and group
baselinescore_Sus<- meanse(Longitudinal_r$Perceived_Susceptibility[Longitudinal_r$Time=="baseline"])
baselinescore_Ser<- meanse(Longitudinal_r$Perceived_Severity[Longitudinal_r$Time=="baseline"])
baselinescore_Bar<- meanse(Longitudinal_r$Perceived_Barrier[Longitudinal_r$Time=="baseline"])
baselinescore_Ben<- meanse(Longitudinal_r$Perceived_Benefit[Longitudinal_r$Time=="baseline"])

outbreakscore_Sus<- meanse(Longitudinal_r$Perceived_Susceptibility[Longitudinal_r$Time=="outbreak"])
outbreakscore_Ser<- meanse(Longitudinal_r$Perceived_Severity[Longitudinal_r$Time=="outbreak"])
outbreakscore_Bar<- meanse(Longitudinal_r$Perceived_Barrier[Longitudinal_r$Time=="outbreak"])
outbreakscore_Ben<- meanse(Longitudinal_r$Perceived_Benefit[Longitudinal_r$Time=="outbreak"])

HPSus<- num_per(Longitudinal_r$Susceptibility_group[Longitudinal_r$Time=="baseline"])
HPSer<- num_per(Longitudinal_r$Severity_group[Longitudinal_r$Time=="baseline"])
HPBar<- num_per(Longitudinal_r$Barrier_group[Longitudinal_r$Time=="baseline"])
HPBen<- num_per(Longitudinal_r$Benefit_group[Longitudinal_r$Time=="baseline"])

HPSus2<- num_per(Longitudinal_r$Susceptibility_group[Longitudinal_r$Time=="outbreak"])
HPSer2<- num_per(Longitudinal_r$Severity_group[Longitudinal_r$Time=="outbreak"])
HPBar2<- num_per(Longitudinal_r$Barrier_group[Longitudinal_r$Time=="outbreak"])
HPBen2<- num_per(Longitudinal_r$Benefit_group[Longitudinal_r$Time=="outbreak"])

Sus_group_same<- length(which(Longitudinal_r$Susceptibility_group_same[Longitudinal_r$Time=="baseline"]==TRUE))/length(Longitudinal_r$Susceptibility_group_same[Longitudinal_r$Time=="baseline" & is.na(Longitudinal_r$Susceptibility_group_same[Longitudinal_r$Time=="baseline"])==FALSE])
Sev_group_same<- length(which(Longitudinal_r$Severity_group_same[Longitudinal_r$Time=="baseline"]==TRUE))/length(Longitudinal_r$Severity_group_same[Longitudinal_r$Time=="baseline" & is.na(Longitudinal_r$Severity_group_same[Longitudinal_r$Time=="baseline"])==FALSE])
Bar_group_same<- length(which(Longitudinal_r$Barrier_group_same[Longitudinal_r$Time=="baseline"]==TRUE))/length(Longitudinal_r$Barrier_group_same[Longitudinal_r$Time=="baseline" & is.na(Longitudinal_r$Barrier_group_same[Longitudinal_r$Time=="baseline"])==FALSE])
Ben_group_same<- length(which(Longitudinal_r$Benefit_group_same[Longitudinal_r$Time=="baseline"]==TRUE))/length(Longitudinal_r$Benefit_group_same[Longitudinal_r$Time=="baseline" & is.na(Longitudinal_r$Benefit_group_same[Longitudinal_r$Time=="baseline"])==FALSE])

Sus_HL<- length(which(Longitudinal_r$Susceptibility_group_baseline=="high"& Longitudinal_r$Susceptibility_group_outbreak=="alow"))/length(Longitudinal_r$Susceptibility_group_same[is.na(Longitudinal_r$Susceptibility_group_same)==FALSE])
Sev_HL<- length(which(Longitudinal_r$Severity_group_baseline=="high"& Longitudinal_r$Severity_group_outbreak=="alow"))/length(Longitudinal_r$Severity_group_same[is.na(Longitudinal_r$Severity_group_same)==FALSE])
Bar_HL<- length(which(Longitudinal_r$Barrier_group_baseline=="high"& Longitudinal_r$Barrier_group_outbreak=="alow"))/length(Longitudinal_r$Barrier_group_same[is.na(Longitudinal_r$Barrier_group_same)==FALSE])
Ben_HL<- length(which(Longitudinal_r$Benefit_group_baseline=="high"& Longitudinal_r$Benefit_group_outbreak=="alow"))/length(Longitudinal_r$Benefit_group_same[is.na(Longitudinal_r$Benefit_group_same)==FALSE])

Sus_LH<- length(which(Longitudinal_r$Susceptibility_group_baseline=="alow"&Longitudinal_r$Susceptibility_group_outbreak=="high"))/length(Longitudinal_r$Susceptibility_group_same[is.na(Longitudinal_r$Susceptibility_group_same)==FALSE])
Sev_LH<- length(which(Longitudinal_r$Severity_group_baseline=="alow"&Longitudinal_r$Severity_group_outbreak=="high"))/length(Longitudinal_r$Severity_group_same[is.na(Longitudinal_r$Severity_group_same)==FALSE])
Bar_LH<- length(which(Longitudinal_r$Barrier_group_baseline=="alow"&Longitudinal_r$Barrier_group_outbreak=="high"))/length(Longitudinal_r$Barrier_group_same[is.na(Longitudinal_r$Barrier_group_same)==FALSE])
Ben_LH<- length(which(Longitudinal_r$Benefit_group_baseline=="alow"&Longitudinal_r$Benefit_group_outbreak=="high"))/length(Longitudinal_r$Benefit_group_same[is.na(Longitudinal_r$Benefit_group_same)==FALSE])

#Table 2
HBM_table<- data.frame(belief= c("Susceptibility","Severity","Barrier","Benefit"),
                       baselinescore= c(baselinescore_Sus,baselinescore_Ser,baselinescore_Bar,baselinescore_Ben),
                       baseline_high_proportion= c(HPSus, HPSer, HPBar, HPBen),
                       outbreakscore= c(outbreakscore_Sus,outbreakscore_Ser,outbreakscore_Bar,outbreakscore_Ben),
                       outbreak_high_proportion= c(HPSus2, HPSer2, HPBar2, HPBen2),
                       same= c(Sus_group_same, Sev_group_same, Bar_group_same, Ben_group_same),
                       High_to_Low= c(Sus_HL,Sev_HL,Bar_HL,Ben_HL),
                       Low_to_High= c(Sus_LH,Sev_LH,Bar_LH,Ben_LH))

#McNemar's Chi-squared Test for Contact Data
mcnemar_results<- data.frame(p_value=c(mcnemar(Longitudinal_r$Susceptibility_group[Longitudinal_r$Time=="baseline"], Longitudinal_r$Susceptibility_group[Longitudinal_r$Time=="outbreak"]),
                                       mcnemar(Longitudinal_r$Severity_group_baseline[Longitudinal_r$Time=="baseline"], Longitudinal_r$Severity_group_outbreak[Longitudinal_r$Time=="outbreak"]),
                                       mcnemar(Longitudinal_r$Barrier_group[Longitudinal_r$Time=="baseline"], Longitudinal_r$Barrier_group[Longitudinal_r$Time=="outbreak"]),
                                       mcnemar(Longitudinal_r$Benefit_group[Longitudinal_r$Time=="baseline"], Longitudinal_r$Benefit_group[Longitudinal_r$Time=="outbreak"])))


#######1-2 Belief Correlation#######
#Table S2
new_col<- c("Perceived_Susceptibility","Perceived_Severity","Perceived_Barrier","Perceived_Benefit")

dat_long<- drop_na(Longitudinal_r[,c(1,19:22)])
res_score_long<- rcorr(as.matrix(dat_long[,new_col]),type=c("spearman"))
round(res_score_long$P, 3)
round(res_score_long$r, 3)

dat_cross<- drop_na(Cross[,c(1,11,13,15,17)])
res_score_cross<- rcorr(as.matrix(dat_cross[,new_col]),type=c("spearman"))
round(res_score_cross$P, 3)
round(res_score_cross$r, 3)

dat_all<- rbind(x = dat_long[,new_col], y = dat_cross[,new_col])
res_score_all<- rcorr(as.matrix(dat_all),type=c("spearman"))
round(res_score_all$P, 3)
round(res_score_all$r, 3)

#########2-Longitudinal_GEE#########
#Table 5
del<- which(is.na(Longitudinal_r$Susceptibility_group)==TRUE|is.na(Longitudinal_r$Severity_group)==TRUE|is.na(Longitudinal_r$Barrier_group)==TRUE|is.na(Longitudinal_r$Benefit_group)==TRUE)

gee_all_To<- geeglm(Total_Contact ~ as.numeric(HHsize) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time) + as.factor(Susceptibility_group) + as.factor(Severity_group) + as.factor(Barrier_group) + as.factor(Benefit_group),
                    data = Longitudinal_r[-del,], id = ID, 
                    family = poisson, corstr = "ar1")
vif(gee_all_To)

gee_all_HH<- geeglm(Household_Contact ~ as.numeric(HHsize) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time) + as.factor(Susceptibility_group) + as.factor(Severity_group) + as.factor(Barrier_group) + as.factor(Benefit_group),
                    data = Longitudinal_r[-del,], id = ID, 
                    family = poisson, corstr = "ar1")
vif(gee_all_HH)

gee_all_NHH<- geeglm(Nonhousehold_Contact ~ as.numeric(HHsize) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time) + as.factor(Susceptibility_group) + as.factor(Severity_group) + as.factor(Barrier_group) + as.factor(Benefit_group),
                     data = Longitudinal_r[-del,], id = ID, 
                     family = poisson, corstr = "ar1")
vif(gee_all_NHH)

geetable(gee_all_To)
geetable(gee_all_HH)
geetable(gee_all_NHH)
#Table S4
gee_all_To_X<- geeglm(Total_Contact ~ as.numeric(HHsize) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Time)*as.factor(Severity_group) + as.factor(Time)*as.factor(Barrier_group) + as.factor(Time)*as.factor(Benefit_group),
                      data = Longitudinal_r[-del,], id = ID, 
                      family = poisson, corstr = "ar1")                   
vif(gee_all_To_X)
gee_all_To_X<- geeglm(Total_Contact ~ as.numeric(HHsize) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Severity_group) + as.factor(Time)*as.factor(Barrier_group) + as.factor(Benefit_group),
                      data = Longitudinal_r[-del,], id = ID, 
                      family = poisson, corstr = "ar1")                   
vif(gee_all_To_X)

gee_all_HH_X<- geeglm(Household_Contact ~ as.numeric(HHsize) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Time)*as.factor(Severity_group) + as.factor(Time)*as.factor(Barrier_group) + as.factor(Time)*as.factor(Benefit_group),
                      data = Longitudinal_r[-del,], id = ID, 
                      family = poisson, corstr = "ar1")
vif(gee_all_HH_X)
gee_all_HH_X<- geeglm(Household_Contact ~ as.numeric(HHsize) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Severity_group) + as.factor(Time)*as.factor(Barrier_group) + as.factor(Benefit_group),
                      data = Longitudinal_r[-del,], id = ID, 
                      family = poisson, corstr = "ar1")
vif(gee_all_HH_X)

gee_all_NHH_X<- geeglm(Nonhousehold_Contact ~ as.numeric(HHsize) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Time)*as.factor(Severity_group) + as.factor(Time)*as.factor(Barrier_group) + as.factor(Time)*as.factor(Benefit_group),
                       data = Longitudinal_r[-del,], id = ID, 
                       family = poisson, corstr = "ar1")
vif(gee_all_NHH_X)
gee_all_NHH_X<- geeglm(Nonhousehold_Contact ~ as.numeric(HHsize) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Severity_group) + as.factor(Time)*as.factor(Barrier_group) + as.factor(Benefit_group),
                       data = Longitudinal_r[-del,], id = ID, 
                       family = poisson, corstr = "ar1")
vif(gee_all_NHH_X)

geetable(gee_all_To_X)
geetable(gee_all_HH_X)
geetable(gee_all_NHH_X)

##Predicted value##
#Figure 2A
test_df<- data.frame(
  Time = c(rep("baseline",16),rep("outbreak",16)),
  Susceptibility_group = c(rep("alow",8),rep("high",8),rep("alow",8),rep("high",8)),
  Severity_group = c(rep("alow",4),rep("high",4),rep("alow",4),rep("high",4),
                     rep("alow",4),rep("high",4),rep("alow",4),rep("high",4)),
  Barrier_group = c(rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2)),
  Benefit_group = c("alow", "high", "alow", "high", "alow", "high", "alow", "high",
                    "alow", "high", "alow", "high", "alow", "high", "alow", "high"),
  HHsize = rep(4,16),
  SEX = rep("女",16),
  AGEgroup2 = rep("2middle",16),
  EDU = rep("2大學/專",16),
  Region = rep("1North",16)
)

predict_df<-  data.frame(
  Time = c(rep("baseline",16),rep("outbreak",16)),
  Susceptibility_group = c(rep("alow",8),rep("high",8),rep("alow",8),rep("high",8)),
  Severity_group = c(rep("alow",4),rep("high",4),rep("alow",4),rep("high",4),
                     rep("alow",4),rep("high",4),rep("alow",4),rep("high",4)),
  Barrier_group = c(rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2)),
  Benefit_group = c("alow", "high", "alow", "high", "alow", "high", "alow", "high",
                    "alow", "high", "alow", "high", "alow", "high", "alow", "high"),
  Total= predict(gee_all_To_X, test_df),
  Household = predict(gee_all_HH_X, test_df),
  Non_household = predict(gee_all_NHH_X, test_df)
)

contactorder <- c("Total","Household","Non_household")
predict_barrier<- predict_df[c(14,16,30,32),c(-2,-3,-5)]
predict_barrier$Group<- c("Low-barrier-baseline","High-barrier-baseline","Low-barrier-peak","High-barrier-peak")

predict_barrier%>% 
  pivot_longer(cols = Total:Non_household, names_to = "Contacttype", values_to = "predicted_number") %>% 
  ggplot(mapping = aes(x=factor(Contacttype,contactorder)))+ 
  geom_col(aes(y=predicted_number, fill=Group), width=0.7, position="dodge", color="black") +
  labs(x="Contact type",y="Predicted number of contact")+
  scale_y_continuous(limits = c(0, 2))+
  theme(axis.text.x = element_text(size = 12, face="bold"),
        axis.text.y = element_text(size = 12, face="bold"),
        axis.title.x = element_text(size = 12, face="bold",margin = margin(t = 10)), 
        axis.title.y = element_text(size = 12, face="bold",margin = margin(r = 10)),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) +
  scale_x_discrete(labels = c("Non_household" = "Non-household")) -> 
  f1
f1


#########3-Cross-sectional_Description#########
##Belief score and group
peakscore_Sus<- meanse(Cross$Perceived_Susceptibility[Cross$Time=="zpeak"])
peakscore_Ser<- meanse(Cross$Perceived_Severity[Cross$Time=="zpeak"])
peakscore_Bar<- meanse(Cross$Perceived_Barrier[Cross$Time=="zpeak"])
peakscore_Ben<- meanse(Cross$Perceived_Benefit[Cross$Time=="zpeak"])

troughscore_Sus<- meanse(Cross$Perceived_Susceptibility[Cross$Time=="trough"])
troughscore_Ser<- meanse(Cross$Perceived_Severity[Cross$Time=="trough"])
troughscore_Bar<- meanse(Cross$Perceived_Barrier[Cross$Time=="trough"])
troughscore_Ben<- meanse(Cross$Perceived_Benefit[Cross$Time=="trough"])

HPSus3<- num_per(Cross$Susceptibility_group[Cross$Time=="zpeak"])
HPSer3<- num_per(Cross$Severity_group[Cross$Time=="zpeak"])
HPBar3<- num_per(Cross$Barrier_group[Cross$Time=="zpeak"])
HPBen3<- num_per(Cross$Benefit_group[Cross$Time=="zpeak"])

HPSus4<- num_per(Cross$Susceptibility_group[Cross$Time=="trough"])
HPSer4<- num_per(Cross$Severity_group[Cross$Time=="trough"])
HPBar4<- num_per(Cross$Barrier_group[Cross$Time=="trough"])
HPBen4<- num_per(Cross$Benefit_group[Cross$Time=="trough"])

#Table 3
HBM_table2<- data.frame(belief= c("Susceptibility","Severity","Barrier","Benefit"),
                       peakscore= c(peakscore_Sus,peakscore_Ser,peakscore_Bar,peakscore_Ben),
                       peak_high_proportion= c(HPSus3, HPSer3, HPBar3, HPBen3),
                       troughscore= c(troughscore_Sus,troughscore_Ser,troughscore_Bar,troughscore_Ben),
                       trough_high_proportion= c(HPSus4, HPSer4, HPBar4, HPBen4))

chi_results2<- data.frame(p_value=c(chi2(Cross$Susceptibility_group[Cross$Time=="zpeak"], Cross$Susceptibility_group[Cross$Time=="trough"]),
                                   chi2(Cross$Severity_group[Cross$Time=="zpeak"], Cross$Severity_group[Cross$Time=="trough"]),
                                   chi2(Cross$Barrier_group[Cross$Time=="zpeak"], Cross$Barrier_group[Cross$Time=="trough"]),
                                   chi2(Cross$Benefit_group[Cross$Time=="zpeak"], Cross$Benefit_group[Cross$Time=="trough"])))

#######3-1 Contact Behavior#######
#Table 4
baseline_To<- meanse(Longitudinal_r$Total_Contact[Longitudinal_r$Time=="baseline"])
baseline_HH<- meanse(Longitudinal_r$Household_Contact[Longitudinal_r$Time=="baseline"])
baseline_NHH<- meanse(Longitudinal_r$Nonhousehold_Contact[Longitudinal_r$Time=="baseline"])
outbreak_To<- meanse(Longitudinal_r$Total_Contact[Longitudinal_r$Time=="outbreak"])
outbreak_HH<- meanse(Longitudinal_r$Household_Contact[Longitudinal_r$Time=="outbreak"])
outbreak_NHH<- meanse(Longitudinal_r$Nonhousehold_Contact[Longitudinal_r$Time=="outbreak"])

peak_TO<- meanse(as.numeric(Cross$Total_Contact[Cross$Time=="zpeak"]))
peak_HH<- meanse(Cross$Household_Contact[Cross$Time=="zpeak"])
peak_NHH<- meanse(Cross$Nonhousehold_Contact[Cross$Time=="zpeak"])
peak_NE<- meanse(Cross$Necessary_Contact[Cross$Time=="zpeak"])
peak_UNE<- meanse(Cross$Unnecessary_Contact[Cross$Time=="zpeak"])

trough_TO<- meanse(Cross$Total_Contact[Cross$Time=="trough"])
trough_HH<- meanse(Cross$Household_Contact[Cross$Time=="trough"])
trough_NHH<- meanse(Cross$Nonhousehold_Contact[Cross$Time=="trough"])
trough_NE<- meanse(Cross$Necessary_Contact[Cross$Time=="trough"])
trough_UNE<- meanse(Cross$Unnecessary_Contact[Cross$Time=="trough"])

p_value(t.test(Longitudinal_r$Total_Contact_change[Longitudinal_r$Time=="baseline"], mu= 0, weight=Longitudinal_r$ID_WEIGHT))
p_value(t.test(Longitudinal_r$Household_Contact_change[Longitudinal_r$Time=="baseline"], mu= 0, weight=Longitudinal_r$ID_WEIGHT))
p_value(t.test(Longitudinal_r$Nonhousehold_Contact_change[Longitudinal_r$Time=="baseline"], mu= 0, weight=Longitudinal_r$ID_WEIGHT))

p_value(t.test(Total_Contact~ Time, data = Cross, weight=Cross$ID_WEIGHT, alternative = c("two.sided"), var.equal = FALSE))
p_value(t.test(Household_Contact~ Time, data = Cross, weight=Cross$ID_WEIGHT, alternative = c("two.sided"), var.equal = FALSE))
p_value(t.test(Nonhousehold_Contact~ Time, data = Cross, weight=Cross$ID_WEIGHT, alternative = c("two.sided"), var.equal = FALSE))
p_value(t.test(Necessary_Contact~ Time, data = Cross, weight=Cross$ID_WEIGHT, alternative = c("two.sided"), var.equal = FALSE))
p_value(t.test(Unnecessary_Contact~ Time, data = Cross, weight=Cross$ID_WEIGHT, alternative = c("two.sided"), var.equal = FALSE))


#########4-Cross-sectional_ZIP#########
#Table 6
ZIP_all_To<- zeroinfl(Total_Contact ~ as.numeric(Household_size) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time) + as.factor(Susceptibility_group) + as.factor(Severity_group)+ as.factor(Barrier_group)+ as.factor(Benefit_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_To)
ZIP_all_HH<- zeroinfl(Household_Contact~ as.numeric(Household_size) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time) + as.factor(Susceptibility_group) + as.factor(Severity_group)+ as.factor(Barrier_group)+ as.factor(Benefit_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_HH)
ZIP_all_NHH<- zeroinfl(Nonhousehold_Contact~ as.numeric(Household_size) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time) + as.factor(Susceptibility_group) + as.factor(Severity_group)+ as.factor(Barrier_group)+ as.factor(Benefit_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_NHH)
ZIP_all_NE<- zeroinfl(Necessary_Contact~ as.numeric(Household_size) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time) + as.factor(Susceptibility_group) + as.factor(Severity_group)+ as.factor(Barrier_group)+ as.factor(Benefit_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_NE)
ZIP_all_UNE<- zeroinfl(Unnecessary_Contact~ as.numeric(Household_size) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time) + as.factor(Susceptibility_group) + as.factor(Severity_group)+ as.factor(Barrier_group)+ as.factor(Benefit_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_UNE)

ZIP_all_HH<- zeroinfl(Household_Contact~ as.numeric(Household_size) + as.factor(SEX) + as.factor(Region) + as.factor(Time) + as.factor(Susceptibility_group) + as.factor(Severity_group)+ as.factor(Barrier_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_HH)
ZIP_all_UNE<- zeroinfl(Unnecessary_Contact~ as.numeric(Household_size) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(Region) + as.factor(Time) + as.factor(Susceptibility_group) + as.factor(Severity_group)+ as.factor(Barrier_group)+ as.factor(Benefit_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_UNE)

ZIP_S<- list(zip_table(ZIP_all_To),zip_table(ZIP_all_HH),zip_table(ZIP_all_NHH),zip_table(ZIP_all_NE),zip_table(ZIP_all_UNE))


#Table S5
ZIP_all_To_X<- zeroinfl(Total_Contact ~ as.numeric(Household_size) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Time)*as.factor(Severity_group)+ as.factor(Time)*as.factor(Barrier_group)+ as.factor(Benefit_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_To_X)
ZIP_all_HH_X<- zeroinfl(Household_Contact~ as.numeric(Household_size) + as.factor(SEX) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Time)*as.factor(Severity_group)+ as.factor(Time)*as.factor(Barrier_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_HH_X)
ZIP_all_NHH_X<- zeroinfl(Nonhousehold_Contact~ as.numeric(Household_size) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Time)*as.factor(Severity_group)+ as.factor(Time)*as.factor(Barrier_group)+ as.factor(Benefit_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_NHH_X)
ZIP_all_NE_X<- zeroinfl(Necessary_Contact~ as.numeric(Household_size) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(EDU) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Time)*as.factor(Severity_group)+ as.factor(Time)*as.factor(Barrier_group)+ as.factor(Benefit_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_NE_X)
ZIP_all_UNE_X<- zeroinfl(Unnecessary_Contact~ as.numeric(Household_size) + as.factor(SEX) + as.factor(AGEgroup2) + as.factor(Region) + as.factor(Time)*as.factor(Susceptibility_group) + as.factor(Time)*as.factor(Severity_group)+ as.factor(Time)*as.factor(Barrier_group)+ as.factor(Benefit_group), weights=ID_WEIGHT, data = Cross) 
vif(ZIP_all_UNE_X)

ZIP<- list(zip_table(ZIP_all_To_X),zip_table(ZIP_all_HH_X),zip_table(ZIP_all_NHH_X),zip_table(ZIP_all_NE_X),zip_table(ZIP_all_UNE_X))

##Predicted value##
#Figure 2B
test_df<- data.frame(
  Time = c(rep("zpeak",16),rep("trough",16)),
  Susceptibility_group = c(rep("alow",8),rep("high",8),rep("alow",8),rep("high",8)),
  Severity_group = c(rep("alow",4),rep("high",4),rep("alow",4),rep("high",4),
                     rep("alow",4),rep("high",4),rep("alow",4),rep("high",4)),
  Barrier_group = c(rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2)),
  Benefit_group = c("alow", "high", "alow", "high", "alow", "high", "alow", "high",
                    "alow", "high", "alow", "high", "alow", "high", "alow", "high"),
  Household_size = rep(3,16),
  SEX = rep("female",16),
  AGEgroup2 = rep("2middle",16),
  EDU = rep("2大學/專",16),
  Region = rep("1North",16)
)

predict_df<-  data.frame(
  Time = c(rep("zpeak",16),rep("trough",16)),
  Susceptibility_group = c(rep("alow",8),rep("high",8),rep("alow",8),rep("high",8)),
  Severity_group = c(rep("alow",4),rep("high",4),rep("alow",4),rep("high",4),
                     rep("alow",4),rep("high",4),rep("alow",4),rep("high",4)),
  Barrier_group = c(rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2),
                    rep("alow",2),rep("high",2),rep("alow",2),rep("high",2)),
  Benefit_group = c("alow", "high", "alow", "high", "alow", "high", "alow", "high",
                    "alow", "high", "alow", "high", "alow", "high", "alow", "high"),
  Total = predict(ZIP_all_To_X, test_df),
  Household = predict(ZIP_all_HH_X, test_df),
  Non_household = predict(ZIP_all_NHH_X, test_df),
  Necessary = predict(ZIP_all_NE_X, test_df),
  Unnecessary = predict(ZIP_all_UNE_X, test_df)
)

contactorder <- c("Total", "Household", "Non_household", "Necessary", "Unnecessary")
predict_barrier<- predict_df[c(1,3,17,19),c(-2,-3,-5)]
predict_barrier$Group<- c("Low-barrier-peak","High-barrier-peak","Low-barrier-trough","High-barrier-trough")
predict_barrier%>% 
  pivot_longer(cols = Total:Unnecessary, names_to = "Contact_type", values_to = "predicted_number") %>%
  ggplot(mapping = aes(x=factor(Contact_type,contactorder))) + 
  geom_col(aes(y=predicted_number, fill=Group), width=0.7, position="dodge", color="black") +
  labs(x="Contact type",y="Predicted number of contact")+
  scale_y_continuous(limits = c(0, 4),labels = scales::label_comma())+
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold",margin = margin(t = 10)), 
        axis.title.y = element_text(size = 12, face = "bold",margin = margin(r = 10)),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) +
  scale_x_discrete(labels = c("Non_household" = "Non-household")) -> 
  f2
f2

#########5-Confirmed_case_Figure#########
#Figure 1
colnames(confirmedcase)
names(confirmedcase)[names(confirmedcase)=="New"]<- "Newly_confirmed_cases"
confirmedcase$Newly_confirmed_cases<- as.numeric(confirmedcase$Newly_confirmed_cases)
confirmedcase$Date<- as.Date(confirmedcase$Date)

plot(x=confirmedcase$Date, y=confirmedcase$Newly_confirmed_cases,
     type="l", xlab="Date", ylab="Cases", 
     main="Covid-19 Confirmed Cases", col="aquamarine4")

Rects = data.frame(x1 = c(as.Date("2020-07-18"),as.Date("2021-06-07"),as.Date("2022-09-19"),as.Date("2022-11-26")), 
                   x2 = c(as.Date("2020-09-16"),as.Date("2021-06-14"),as.Date("2022-10-11"),as.Date("2022-11-30")), 
                   y1 = c(-Inf,-Inf), y2= c(Inf,Inf))

ggplot(confirmedcase, aes(x=Date, y=Newly_confirmed_cases))+
  geom_line(color="aquamarine4")+
  ggtitle("Covid-19 Confirmed Cases in Taiwan")+
  labs(x="Date", y="Weekly confirmed cases")+
  scale_x_date(date_breaks = "3 month",date_labels = "%Y/%m",date_minor_breaks = "3 month")+
  geom_rect(data = Rects,
            inherit.aes = FALSE,
            mapping = aes(xmin = x1, xmax = x2,
                          ymin = y1, ymax = y2),
            color = "transparent",
            fill = "red",
            alpha = .2) +
  theme(axis.text = element_text(size=10, face="bold"),
        axis.title.x = element_text(size=12, face="bold",margin = margin(t = 12)), 
        axis.title.y = element_text(size=12, face="bold",margin = margin(r = 12)),
        plot.margin = margin(t = 10, r = 12, b = 10, l = 10, unit = "pt"))  -> p

p


Rects2<-data.frame(x1 = c(as.Date("2021-05-19")), 
                   x2 = c(as.Date("2021-07-27")), 
                   y1 = c(-Inf,-Inf), y2= c(-500,-500))
p2<- p+ 
  geom_rect(data = Rects2,
            inherit.aes = FALSE,
            mapping = aes(xmin = x1, xmax = x2,
                          ymin = y1, ymax = y2),
            color = "transparent",
            fill = "blue4",
            alpha = .2) +
  annotate("text", x=as.Date("2021-06-19"), y=-2000, label="Third-level alert", size=6,color="blue4")+
  
  annotate("segment", x=as.Date("2020-08-18"), y = 6000, xend = as.Date("2020-08-18"), yend = 1000,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  annotate("text", x=as.Date("2020-08-18"), y=9000, label="L-baseline", size=7,color="red")+
  
  annotate("segment", x=as.Date("2021-06-10"), y = 7000, xend = as.Date("2021-06-10"), yend = 2000,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  annotate("text", x=as.Date("2021-06-10"), y=10000, label="L-peak", size=7,color="red")+
  
  annotate("segment", x=as.Date("2022-09-30"), y = 60000, xend = as.Date("2022-09-30"), yend = 55000,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  annotate("text", x=as.Date("2022-09-30"), y=65000, label="C-peak", size=7,color="red")+
  
  annotate("segment", x=as.Date("2022-11-28"), y = 25000, xend = as.Date("2022-11-28"), yend = 20000,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  annotate("text", x=as.Date("2023-01-15"), y=30000, label="C-trough", size=7,color="red")
p2


#########6-Social_Contact_Matrix#########
#Figure S2

#######6-1 peak########
rm(contactagedata,contact_info,contact_info_temp)
par<- Cross[Cross$Time=="zpeak",c(1:9)]

set.seed(123)
for (i in c(1:nrow(par))){
  if(par$AGE[i]=="20-24歲"){
    par$AGE[i]<- sample(c(20:24),1)
  }else if(par$AGE[i]=="25-29歲"){
    par$AGE[i]<- sample(c(25:29),1)
  }else if(par$AGE[i]=="30-34歲"){
    par$AGE[i]<- sample(c(30:34),1)
  }else if(par$AGE[i]=="35-39歲"){
    par$AGE[i]<- sample(c(35:39),1)
  }else if(par$AGE[i]=="40-44歲"){
    par$AGE[i]<- sample(c(40:44),1)
  }else if(par$AGE[i]=="45-49歲"){
    par$AGE[i]<- sample(c(45:49),1)
  }else if(par$AGE[i]=="50-54歲"){
    par$AGE[i]<- sample(c(50:54),1)
  }else if(par$AGE[i]=="55-59歲"){
    par$AGE[i]<- sample(c(55:59),1)
  }else if(par$AGE[i]=="60-64歲"){
    par$AGE[i]<- sample(c(60:64),1)
  }else if(par$AGE[i]=="65-69歲"){
    par$AGE[i]<- sample(c(65:69),1)
  }else if(par$AGE[i]=="70-74歲"){
    par$AGE[i]<- sample(c(70:74),1)
  }else if(par$AGE[i]=="75歲以上"){
    par$AGE[i]<- sample(c(75:80),1)
  }
}

summary(as.factor(par$AGE))
names(par)[names(par)=="ID"]<- "part_id"
names(par)[names(par)=="AGE"]<- "part_age"
names(par)[names(par)=="SEX"]<- "part_gender"
par$country<-rep("Taiwan")

contactagedata<- Cross[Cross$Time=="zpeak",c(1,39:50)]
contact_info<- pivot_longer(
  contactagedata, 
  cols = starts_with("Q14A_age"), 
  names_to = "contact_age", 
  values_to = "cnt_age_exact",
  values_drop_na = TRUE
)

if(!require('sjmisc')) {
  install.packages('sjmisc')
  library('sjmisc')
}

contact_info$cnt_age_est_min<-rep(NA)
contact_info$cnt_age_est_max<-rep(NA)

contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="0-5歲",0)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="0-5歲",5)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="6-10歲",6)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="6-10歲",10)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="11-15歲",11)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="11-15歲",15)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="16-20歲",16)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="16-20歲",20)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="21-25歲",21)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="21-25歲",25)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="26-30歲",26)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="26-30歲",30)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="31-35歲",31)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="31-35歲",31)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="36-40歲",36)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="36-40歲",40)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="41-45歲",41)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="41-45歲",45)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="46-50歲",46)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="46-50歲",50)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="51-55歲",51)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="51-55歲",55)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="56-60歲",56)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="56-60歲",60)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="61-65歲",61)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="61-65歲",61)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="66-70歲",66)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="66-70歲",70)
contact_info$cnt_age_est_min <- replace(contact_info$cnt_age_est_min,contact_info$cnt_age_exact=="70歲以上",73)
contact_info$cnt_age_est_max <- replace(contact_info$cnt_age_est_max,contact_info$cnt_age_exact=="70歲以上",78)
names(contact_info)[names(contact_info)=="ID"]<- "part_id"
contact_info$cnt_age_exact<-rep(NA)
##Total
peak_survey<- survey(par,contact_info)
cm_peak <- contact_matrix(peak_survey, age.limits = c(0,20,30,40,50,60,70),
                          symmetric = FALSE,
                          missing.contact.age = "remove")
mx_peak <- cm_peak$matrix
matrix_plot(mx_peak, main = "Total contacts in Peak period",
            max.legend=1.5,cex.lab=1,cex.axis= 0.8,legend.width=1)

#######6-2 trough#########
rm(contactagedata2,contact_info2,contact_info_temp2)
par2<- Cross[Cross$Time=="trough",c(1:9)]

set.seed(123)
for (i in c(1:nrow(par2))){
  if(par2$AGE[i]=="20-24歲"){
    par2$AGE[i]<- sample(c(20:24),1)
  }else if(par2$AGE[i]=="25-29歲"){
    par2$AGE[i]<- sample(c(25:29),1)
  }else if(par2$AGE[i]=="30-34歲"){
    par2$AGE[i]<- sample(c(30:34),1)
  }else if(par2$AGE[i]=="35-39歲"){
    par2$AGE[i]<- sample(c(35:39),1)
  }else if(par2$AGE[i]=="40-44歲"){
    par2$AGE[i]<- sample(c(40:44),1)
  }else if(par2$AGE[i]=="45-49歲"){
    par2$AGE[i]<- sample(c(45:49),1)
  }else if(par2$AGE[i]=="50-54歲"){
    par2$AGE[i]<- sample(c(50:54),1)
  }else if(par2$AGE[i]=="55-59歲"){
    par2$AGE[i]<- sample(c(55:59),1)
  }else if(par2$AGE[i]=="60-64歲"){
    par2$AGE[i]<- sample(c(60:64),1)
  }else if(par2$AGE[i]=="65-69歲"){
    par2$AGE[i]<- sample(c(65:69),1)
  }else if(par2$AGE[i]=="70-74歲"){
    par2$AGE[i]<- sample(c(70:74),1)
  }else if(par2$AGE[i]=="75歲以上"){
    par2$AGE[i]<- sample(c(75:80),1)
  }
}

summary(as.factor(par2$AGE))
names(par2)[names(par2)=="ID"]<- "part_id"
names(par2)[names(par2)=="AGE"]<- "part_age"
names(par2)[names(par2)=="SEX"]<- "part_gender"
par2$country<-rep("Taiwan")

contactagedata2<- Cross[Cross$Time=="trough",c(1,39:50)]
contact_info2<- pivot_longer(
  contactagedata2, 
  cols = starts_with("Q14A_age"), 
  names_to = "contact_age", 
  values_to = "cnt_age_exact",
  values_drop_na = TRUE
)

contact_info2$cnt_age_est_min<-rep(NA)
contact_info2$cnt_age_est_max<-rep(NA)

contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==1,0)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==1,5)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==2,6)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==2,10)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==3,11)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==3,15)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==4,16)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==4,20)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==5,21)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==5,25)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==6,26)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==6,30)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==7,31)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==7,31)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==8,36)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==8,40)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==9,41)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==9,45)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==10,46)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==10,50)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==11,51)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==11,55)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==12,56)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==12,60)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==13,61)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==13,61)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==14,66)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==14,70)
contact_info2$cnt_age_est_min <- replace(contact_info2$cnt_age_est_min,contact_info2$cnt_age_exact==15,73)
contact_info2$cnt_age_est_max <- replace(contact_info2$cnt_age_est_max,contact_info2$cnt_age_exact==15,78)
names(contact_info2)[names(contact_info2)=="ID"]<- "part_id"
contact_info2$cnt_age_exact<-rep(NA)

##Total
trough_survey<- survey(par2,contact_info2)
cm_trough <- contact_matrix(trough_survey, age.limits = c(0,20,30,40,50,60,70),
                            symmetric = FALSE,
                            missing.contact.age = "remove")
mx_trough <- cm_trough$matrix
matrix_plot(mx_trough, main = "Total contacts in Trough period",
            max.legend=1.5,cex.lab=1,cex.axis= 0.8,legend.width=1)

