setwd("~/MASTER Sciences de la Mer/Stage M2/projets/david_shanafelt3.0/analyse_R/TL1")


library(reshape2)
library(corrplot)
library(ggplot2)
library(tidyr)
library(cowplot)
library(latex2exp) # Permet d'écrire en latex => Tex($\\..)
library(scales)
library(dplyr) 
library(wesanderson)
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


theme_big<-theme_bw()+
  theme(text = element_text(size=25, family='Times'),
        axis.text = element_text(size=25),
        axis.line = element_line(),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"))

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

leg_net_effect = TeX(("Net effect of $\ $ $\\B_1$ on $\\B_0$"))






y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))

x_axis_log10<-scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))

y_axis_log10<-scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))




######RECYCLAGE iFIXE/COMPENSé #################


# Graphique avec Inva/Bio/Var

datfixe = read.table(file = "data_recy_TL1_pertuB1.txt" , fill = TRUE, sep = ";", header = TRUE)
datfixe = datfixe[,-54]


all_recy_fixe = datfixe[,c(11,12,23:24,31,34)]
# all_recy_fixe = all_recy_fixe[which((all_recy_fixe$D == 0)|(all_recy_fixe$D==0.0001)|(all_recy_fixe$D==0.1)|(all_recy_fixe$D==1)),]
all_recy_fixe = all_recy_fixe[which((all_recy_fixe$D == 0)|(all_recy_fixe$D==1)),]

colnames(all_recy_fixe) = c("lambda","D","1","2","11","22")


CV1 = all_recy_fixe$`1` / sqrt(all_recy_fixe$`11`)
CV2 = all_recy_fixe$`2` / sqrt(all_recy_fixe$`22`)


all_recy_fixe = cbind(all_recy_fixe,CV1,CV2)
all_recy_fixe$`11` = sqrt(all_recy_fixe$`11`)
all_recy_fixe$`22` = sqrt(all_recy_fixe$`22`)

comparti = rep(c(rep(1,length(all_recy_fixe$lambda)),rep(2,length(all_recy_fixe$lambda))),3)
forme = c(rep(1,length(all_recy_fixe$lambda)*2),rep(2,length(all_recy_fixe$lambda)*2),rep(3,length(all_recy_fixe$lambda)*2))



all_recy_fixe = melt(all_recy_fixe,id.vars = c("lambda", "D"))
all_recy_fixe = cbind(all_recy_fixe,comparti,forme)


all_recy_fixe$D = as.factor(all_recy_fixe$D)
levels(all_recy_fixe$D)<-c(D0,D1)

all_recy_fixe$forme = as.factor(all_recy_fixe$forme)
all_recy_fixe$forme = factor(all_recy_fixe$forme,levels(all_recy_fixe$forme)[c(3,2,1)])
all_recy_fixe$comparti = as.factor(all_recy_fixe$comparti)
levels(all_recy_fixe$comparti) = c(N1_s_en,N2_s_en)


plot_all_recyfixe = ggplot(data = all_recy_fixe,aes( x = lambda , y = value,color = forme, linetype = forme))+
  geom_line(size = 1.5)+
  
  scale_linetype_discrete(name = "", 
                          labels = c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))) ) +
  
  facet_grid(D~comparti, labeller = label_parsed,scales = "free")+
  scale_colour_manual(name="",values=c("blue","peru","slategrey"),
                      labels = c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))) ) + # Set legend title
  xlab(leg_recy_en) + ylab(leg_in_en) + # Set axis labels
  theme + y_axis_log10

  # ggsave("all_recy_TL1_pertuB2_comp_ifixe_short_changesize.pdf",width = 10.5, height = 7.92)




### J12

j12fixe = datfixe[,c(11,12,47)]

# j12fixe = j12fixe[which((j12fixe$D == 0)|(j12fixe$D==0.0001)|(j12fixe$D==0.1)|(j12fixe$D==1)),]
j12fixe = j12fixe[which((j12fixe$D == 0)|(j12fixe$D==1)),]
j12fixe$D = as.factor(j12fixe$D)
# levels(j12fixe$D)<-c(D0,D00001,D01,D1)
levels(j12fixe$D)<-c(D0,D1)

plot_j12_recyfixe = ggplot(data = j12fixe,aes( x = lambda , y = J12))+
  geom_line(size = 1.5)+
  
  facet_grid(D~., labeller = label_parsed ,scales = "fixed")+
  scale_colour_manual(name="Niveaux trophiques",values=c("black", "blue","red","green4","purple"), labels = c(expression(B[1]),expression(B[2]))) + # Set legend title
  xlab(leg_recy_en) + ylab(leg_net_effect) + # Set axis labels
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed", size = 1)+
  theme_big

 # ggsave("j12_recy_TL1_pertuB2_comp_ifixe_short_changesize.pdf", width = 9,height = 6.92)





####### SELF-REGULATION COMPENSé##########


dafixe = read.table(file = "data_D_TL1_pertuB1.txt" , fill = TRUE, sep = ";", header = TRUE)
dafixe = dafixe[,-54]


# Graphique Inva/Bio/Var


all_Dfixe = dafixe[,c(11,12,23:24,31,34)]

all_Dfixe = all_Dfixe[which((all_Dfixe$lambda == 0)|(all_Dfixe$lambda == 0.9)),]

colnames(all_Dfixe) = c("lambda","D","1","2","11","22")


CV1 = all_Dfixe$`1` / sqrt(all_Dfixe$`11`)
CV2 = all_Dfixe$`2` / sqrt(all_Dfixe$`22`)

all_Dfixe = cbind(all_Dfixe,CV1,CV2)
all_Dfixe$`11` = sqrt(all_Dfixe$`11`)
all_Dfixe$`22` = sqrt(all_Dfixe$`22`)

comparti = rep(c(rep(1,length(all_Dfixe$lambda)),rep(2,length(all_Dfixe$lambda))),3)
forme = c(rep(1,length(all_Dfixe$lambda)*2),rep(2,length(all_Dfixe$lambda)*2),rep(3,length(all_Dfixe$lambda)*2))



all_Dfixe = melt(all_Dfixe,id.vars = c("lambda", "D"))

all_Dfixe = cbind(all_Dfixe, comparti,forme)

all_Dfixe$lambda = as.factor(all_Dfixe$lambda)
levels(all_Dfixe$lambda)<-c(lam0,lam09)
all_Dfixe$lambda = factor(all_Dfixe$lambda,levels(all_Dfixe$lambda)[c(3,2,1)])

all_Dfixe$forme = as.factor(all_Dfixe$forme)
all_Dfixe$forme = factor(all_Dfixe$forme,levels(all_Dfixe$forme)[c(3,2,1)])
all_Dfixe$comparti = as.factor(all_Dfixe$comparti)
levels(all_Dfixe$comparti) = c(N1_s_en,N2_s_en)


plot_inva_Dfixe = ggplot(data = all_Dfixe,aes( x = D , y = value, color = forme, linetype =forme))+
  geom_line(size = 1.5)+
  
  scale_linetype_discrete(name = "", 
                          labels = c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))) ) +
  
  facet_grid(lambda~comparti, labeller = label_parsed, scales = "fixed")+
  scale_colour_manual(name="",values=c("blue","peru","slategrey"), 
                      labels =  c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))) ) + # Set legend title
  xlab(leg_d_en) + ylab(leg_in_en) + # Set axis labels
  theme + x_axis_log10 + y_axis_log10

 # ggsave("all_D_TL1_pertuB2_comp_ifixe_changesize.pdf", width = 13, height = 8.69)








