setwd("~/MASTER Sciences de la Mer/Stage M2/projets/david_shanafelt3.0/analyse_R/TL2")





library(reshape2)
library(corrplot)
library(ggplot2)
library(tidyr)
library(cowplot)
library(latex2exp) # Permet d'écrire en latex => Tex($\\..)
library(scales)
library(dplyr) 




######PLOT OPTION #####

windowsFonts(Times=windowsFont("TT Times New Roman")) # Need pour changer la police
# theme
theme<-theme_bw()+
  theme(text = element_text(size=20, family='Times'),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"))


theme_dec<-theme_bw()+
  theme(text = element_text(size=20, family='Times'),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"),
        panel.spacing.x = unit(6,"mm"))

theme_big<-theme_bw()+
  theme(text = element_text(size=25, family='Times'),
        axis.text = element_text(size=25),
        axis.line = element_line(),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"),
        panel.spacing.x = unit(4,"mm"))


theme_matrix<-theme+theme(panel.grid.major = element_blank(),
                          axis.ticks=element_blank(),
                          axis.title=element_blank())



size_guide <- guides(shape = guide_legend(override.aes = list(size = 1.5)))



lam0 <- TeX('$\\lambda = 0')
lam05 <- TeX('$\\lambda = 0.5')
lam1 <- TeX('$\\lambda = 1')
lam045 <- TeX('$\\lambda = 0.45')
lam09 <- TeX('$\\lambda = 0.9')
D0 = TeX(("$D = 0"))
D00001 = TeX(("$D = 0.0001"))
D0001 = TeX(("$D = 0.001"))
D001 = TeX(("$D = 0.01"))
D01= TeX(("$D = 0.1"))
D1 = TeX(("$D = 1"))
D05 = TeX(("$D = 0.5"))


leg_recy_en =TeX("Recycling efficiency $\\lambda$")
leg_inva_en =TeX("Invariability value  $\\S^*_i$")
leg_bio_en = TeX("Biomass at equilibrium $\\B^*_i$")
leg_var_en = TeX("Standard deviation $\\sigma_i^*$")
leg_d_en = TeX("Self-regulation coefficient $\\D$")
leg_in_en = TeX("$Stability$")

N1_s_en = TeX(("Mineral nutrient $\\B_0$"))
N2_s_en = TeX(("Primary producer $\\B_1$"))
N3_s_en = TeX(("Herbivore $\\B_2$"))





y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))

x_axis_log10<-scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))

y_axis_log10<-scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))





####### PARTIE RECYCLAGE ###########

datfix = read.table(file = "data_recy_TL2_pertuB2.txt" , fill = TRUE, sep = ";", header = TRUE)
datfix = datfix[,-82]


# Graphique Inva/Bio/Var

all_recy_fix = datfix[,c(11,12,23:25,35,39,43)]


all_recy_fix = all_recy_fix[which((all_recy_fix$D == 0)|(all_recy_fix$D==0.0001)|(all_recy_fix$D==0.1)),]

colnames(all_recy_fix) = c("lambda","D","1","2","3","11","22","33")

CV1 = all_recy_fix$`1` / sqrt(all_recy_fix$`11`)
CV2 = all_recy_fix$`2` / sqrt(all_recy_fix$`22`)
CV3 = all_recy_fix$`3` / sqrt(all_recy_fix$`33`)

all_recy_fix = cbind(all_recy_fix,CV1,CV2,CV3)
# Inverse de l'écart-type
all_recy_fix$`11` = sqrt(all_recy_fix$`11`)
all_recy_fix$`22` = sqrt(all_recy_fix$`22`)
all_recy_fix$`33` = sqrt(all_recy_fix$`33`)



comparti = rep(c(rep(1,length(all_recy_fix$lambda)),rep(2,length(all_recy_fix$lambda)),rep(3,length(all_recy_fix$lambda))),3)
forme = c(rep(1,length(all_recy_fix$lambda)*3),rep(2,length(all_recy_fix$lambda)*3),rep(3,length(all_recy_fix$lambda)*3))


all_recy_fix = melt(all_recy_fix,id.vars = c("lambda", "D"))
all_recy_fix = cbind(all_recy_fix,comparti,forme)


all_recy_fix$D = as.factor(all_recy_fix$D)
levels(all_recy_fix$D)<-c(D0,D00001,D01)

all_recy_fix$forme = as.factor(all_recy_fix$forme)
all_recy_fix$forme = factor(all_recy_fix$forme,levels(all_recy_fix$forme)[c(3,2,1)])

all_recy_fix$comparti = as.factor(all_recy_fix$comparti)
levels(all_recy_fix$comparti) = c(N1_s_en,N2_s_en,N3_s_en)


plot_all_recy_fixe = ggplot(data = all_recy_fix,aes( x = lambda , y = value, color = forme, linetype = forme))+
  geom_line(size = 1.5)+
  
  scale_linetype_discrete(name = "",
                          labels = c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))))+
  
  facet_grid(D~comparti, labeller = label_parsed,scales = "free")+
  scale_colour_manual(name="",values=c("blue", "peru","slategrey"),
                      labels = c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))) ) + # Set legend title
  xlab(leg_recy_en) + ylab(leg_in_en) + # Set axis labels
  theme_dec + y_axis_log10 


  # ggsave("all_recy_TL2_pertuB3_ifixe_changesize.pdf", width = 14.5, height = 8.69)


#### Explication INVA B1/B2 via la jacobienne

ex_recy_fix = datfix[,c(11,12,65,66,67,69)]
# ex_recy_fix = ex_recy_fix[which((ex_recy_fix$D == 0)|(ex_recy_fix$D==0.0001)|(ex_recy_fix$D==0.1)|(ex_recy_fix$D==1)),]
ex_recy_fix = ex_recy_fix[which(ex_recy_fix$D==0.1),]


expli_J23J12_B1 =  ex_recy_fix$J23 * ex_recy_fix$J12
expli_J13_B1 =    ex_recy_fix$J13

expli_J13J21_B2 = ex_recy_fix$J13 * ex_recy_fix$J21
expli_J23_B2 = ex_recy_fix$J23


ex_recy_fix = ex_recy_fix[,-(3:6)]
ex_recy_fix = cbind(ex_recy_fix,expli_J23J12_B1,expli_J13_B1, expli_J13J21_B2, expli_J23_B2)


comparti = c(rep(rep(1,length(ex_recy_fix$lambda)),2),rep(rep(2,length(ex_recy_fix$lambda)),2))


ex_recy_fix = melt(ex_recy_fix, id.vars = c("lambda","D"))
ex_recy_fix = cbind(ex_recy_fix,comparti)


ex_recy_fix$D = as.factor(ex_recy_fix$D)
# levels(ex_recy_fix$D)<-c(D0,D00001,D01,D1)
levels(ex_recy_fix$D)<-c(D01)

ex_recy_fix$comparti = as.factor(ex_recy_fix$comparti)
levels(ex_recy_fix$comparti) = c(TeX("Effect of $\ B_2$ on $\ B^_0$"),TeX("Effect of $\ B_2$ on $\ B^_1$"))

plot_explijaco_recy_fix = ggplot(data = ex_recy_fix,aes( x = lambda , y = value, color = variable))+
  geom_line(size = 1.5)+
  facet_grid(D~comparti, labeller = label_parsed,scales = "free")+
  
  scale_linetype_discrete(name = "",guide = FALSE) +
  scale_colour_manual(name="",values=c("red","green4","blue","goldenrod3"), 
                      labels = c(expression(paste(J[12]*J[0][1])), expression(paste(J[0][2])),
                                 expression(paste(J[0][2]*J[10])), expression(paste(J[12])))) + # Set legend title
  xlab(leg_recy_en) + ylab("Net direct or indirect effect") + # Set axis labels
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed", size = 1)+
  theme

  # ggsave("expli_recy_B1B2_jacobienne_TL2_pertuB3_ifixe_changesize.pdf", width = 12, height = 5.54)




##########################################################################################################

######### PARTIE SELF REGULATION #######


dafix = read.table(file = "data_D_TL2_pertuB2.txt" , fill = TRUE, sep = ";", header = TRUE)
dafix = dafix[,-82]


# Graphique Inva/Bio/Var

all_D_fix = dafix[,c(8,11,12,23:25,35,39,43)]

colnames(all_D_fix) = c("I","lambda","D","1","2","3","11","22","33")


CV1 = all_D_fix$`1` / sqrt(all_D_fix$`11`)
CV2 = all_D_fix$`2` / sqrt(all_D_fix$`22`)
CV3 = all_D_fix$`3` / sqrt(all_D_fix$`33`)

all_D_fix = cbind(all_D_fix,CV1,CV2,CV3)
# Inverse de l'écart-type
all_D_fix$`11` = sqrt(all_D_fix$`11`)
all_D_fix$`22` = sqrt(all_D_fix$`22`)
all_D_fix$`33` = sqrt(all_D_fix$`33`)



comparti = rep(c(rep(1,length(all_D_fix$lambda)),rep(2,length(all_D_fix$lambda)),rep(3,length(all_D_fix$lambda))),3)
forme = c(rep(1,length(all_D_fix$lambda)*3),rep(2,length(all_D_fix$lambda)*3),rep(3,length(all_D_fix$lambda)*3))


all_D_fix = melt(all_D_fix,id.vars = c("I","lambda", "D"))
all_D_fix = cbind(all_D_fix,comparti,forme)

all_D_fix$lambda = as.factor(all_D_fix$lambda)
levels(all_D_fix$lambda)<-c(lam0,lam05,lam1)
all_D_fix$lambda = factor(all_D_fix$lambda,levels(all_D_fix$lambda)[c(3,2,1)])

all_D_fix$forme = as.factor(all_D_fix$forme)
all_D_fix$forme = factor(all_D_fix$forme,levels(all_D_fix$forme)[c(3,2,1)])

all_D_fix$comparti = as.factor(all_D_fix$comparti)
levels(all_D_fix$comparti) = c(N1_s_en,N2_s_en,N3_s_en)



plot_all_D_fix = ggplot(data = all_D_fix,aes( x = D , y = value, color = forme, linetype = forme))+
  geom_line(size = 1.5)+
  
  scale_linetype_discrete(name = "",
                          labels = c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))))+
  
  facet_grid(lambda~comparti, labeller = label_parsed, scales = "fixed")+
  scale_colour_manual(name="",values=c("blue","peru","slategrey"), 
                      labels = c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i])))) + # Set legend title
  xlab(leg_d_en) + ylab(leg_in_en) + # Set axis labels
  theme + x_axis_log10_short + y_axis_log10

 # ggsave("all_D_TL2_pertuB3_ifixe_changesize.pdf", width = 12, height = 9.46)




### Explication INVA B1/B2 via la jacobienne

ex_d_fix = dafix[,c(11,12,65,66,67,69)]




expli_J23J12_B1 =  ex_d_fix$J23 * ex_d_fix$J12
expli_J13_B1 =    ex_d_fix$J13

expli_J13J21_B2 = ex_d_fix$J13 * ex_d_fix$J21
expli_J23_B2 = ex_d_fix$J23


ex_d_fix = ex_d_fix[,-(3:6)]
ex_d_fix = cbind(ex_d_fix,expli_J23J12_B1,expli_J13_B1, expli_J13J21_B2, expli_J23_B2)


comparti = c(rep(rep(1,length(ex_d_fix$lambda)),2),rep(rep(2,length(ex_d_fix$lambda)),2))


ex_d_fix = melt(ex_d_fix, id.vars = c("lambda","D"))
ex_d_fix = cbind(ex_d_fix,comparti)


ex_d_fix$lambda = as.factor(ex_d_fix$lambda)
levels(ex_d_fix$lambda)<-c(lam0,lam05,lam1)
ex_d_fix$lambda = factor(ex_d_fix$lambda,levels(ex_d_fix$lambda)[c(3,2,1)])
ex_d_fix$comparti = as.factor(ex_d_fix$comparti)
levels(ex_d_fix$comparti) = c(TeX("Effect of $\ $ $\\B_2$ on $\\B^_0$"),TeX("Effect of $\ $ $\\B_2$ on $\\B^_1$"))

plot_explijaco_d_fix = ggplot(data = ex_d_fix,aes( x = D , y = value, color = variable))+
  geom_line(size = 1.5)+
  facet_grid(lambda~comparti, labeller = label_parsed,scales = "fixed")+
  
  scale_linetype_discrete(name = "",guide = FALSE) +
  scale_colour_manual(name="",values=c("red","green4","blue","goldenrod3"), 
                      labels = c(expression(paste(J[12]*J[0][1])), expression(paste(J[0][2])),
                                 expression(paste(J[0][2]*J[10])), expression(paste(J[12])))) + # Set legend title
  xlab(leg_d_en) + ylab("Net direct or indirect effect") + # Set axis labels
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed", size = 1)+
  theme + x_axis_log10

# ggsave("expli_D_B1B2_jacobienne_TL2_pertuB3_ifixe_changesize.pdf", width = 10, height = 7.69)

















