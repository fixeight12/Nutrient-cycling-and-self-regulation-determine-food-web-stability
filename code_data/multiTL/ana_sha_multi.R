setwd("~/MASTER Sciences de la Mer/Stage M2/projets/david_shanafelt3.0/analyse_R/multiTL")


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

theme_matrix<-theme+theme(panel.grid.major = element_blank(),
                          axis.ticks=element_blank(),
                          axis.title=element_blank())



size_guide <- guides(shape = guide_legend(override.aes = list(size = 1.5)))



lam0 <- TeX('$\\lambda = 0')
lam05 <- TeX('$\\lambda = 0.5')
lam1 <- TeX('$\\lambda = 1')
lam045 <- TeX('$\\lambda = 0.45')
lam09 <- TeX('$\\lambda = 0.9')
lam07 <- TeX('$\\lambda = 0.7')
D0 = TeX(("$D = 0"))
D00001 = TeX(("$D = 0.0001"))
D0001 = TeX(("$D = 0.001"))
D001 = TeX(("$D = 0.01"))
D003 = TeX(("$D = 0.03"))
D005 = TeX(("$D = 0.05"))
D01= TeX(("$D = 0.1"))
D1 = TeX(("$D = 1"))
D05 = TeX(("$D = 0.5"))

leg_recy_en =TeX("Recycling efficiency $\\lambda$")
leg_inva_en =TeX("Invariability value $\\S_{i}$")
leg_inva_en2 = TeX(string = "Invariability value $\ S_{i}$")
leg_inva_en3 = expression(paste("Invariability value " , sep = " "  , S[i]))


leg_bio_en = TeX("biomass at equilibrium $\\B^*_i$")
leg_var_en = TeX("standard deviation $\\sigma_i^*$")
leg_d_en = TeX("self-regulation coefficient $\\D$")
leg_in_en = TeX("$\\B_i^*$    $\\sigma_i^*$    $\\S_i^*$")

N1_s_en = TeX(("Mineral nutrient $\\B_0$"))
N2_s_en = TeX(("Primary producer $\\B_1$"))
N3_s_en = TeX(("Herbivores $\\B_2$"))
N4_s_en = TeX(("Carnivores $\\B_3$"))

T0 = TeX(("TL0"))
T1 = TeX(("TL1"))
T2 = TeX(("TL2"))
T3 = TeX(("TL3"))





y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))

x_axis_log10<-scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))

y_axis_log10<-scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))




############### PARTIE iFIXE ###################


dat0 =  data.frame (n  = c(0,0,0,0,0,0),
                  lambda = c(1,0.7,0.0,1,0.7,0.0),
                  D = c(0.0,0.0,0.0,0.1,0.1,0.1),
                  Eq_1 = c(1,1,1,1,1,1),
                  V11= c(0.0125,0.0125,0.0125,0.0125,0.0125,0.0125)
)



dat1 = read.table(file = "data_multiTL1.txt" , fill = TRUE, sep = ";", header = TRUE)
dat1 = dat1[,-54]

dat2 = read.table(file = "data_multiTL2.txt" , fill = TRUE, sep = ";", header = TRUE)
dat2 = dat2[,-82]

dat3 = read.table(file = "data_multiTL3.txt" , fill = TRUE, sep = ";", header = TRUE)
dat3 = dat3[,-118]


################ Partie ifixe mais pour LAMBDA = 0.75########.


################# Graphique idée de Pierre#################

# ATTENTION MODIFICATION DES STA ICI !!!!!!!!!

## mise en place de la mesure d'invariabilité TL0

sta0 = dat0
colnames(sta0) = c("n","lambda","D","1","11")

CV1 = sta0$`1` / sqrt(sta0$`11`)

sta0 = cbind(sta0,CV1)

## mise en place de la mesure d'invariabilité TL1

sta1 = dat1
sta1 = sta1[,c(3,11,12,23,24,31,34)]
colnames(sta1) = c("n","lambda","D","1","2","11","22")

CV1 = sta1$`1` / sqrt(sta1$`11`)
CV2 = sta1$`2` / sqrt(sta1$`22`)
sta1 = cbind(sta1,CV1,CV2)


## mise en place de la mesure d'invariabilité TL2

sta2 = dat2
sta2 = sta2[,c(3,11,12,23,24,25,35,39,43)]
colnames(sta2) = c("n","lambda","D","1","2","3","11","22","33")

CV1 = sta2$`1` / sqrt(sta2$`11`)
CV2 = sta2$`2` / sqrt(sta2$`22`)
CV3 = sta2$`3` / sqrt(sta2$`33`)
sta2 = cbind(sta2,CV1,CV2,CV3)


## mise en place de la mesure d'invariabilité TL3

sta3 = dat3
sta3 = sta3[,c(3,11,12,23,24,25,26,39,44,49,54)]
colnames(sta3) = c("n","lambda","D","1","2","3","4","11","22","33","44")

CV1 = sta3$`1` / sqrt(sta3$`11`)
CV2 = sta3$`2` / sqrt(sta3$`22`)
CV3 = sta3$`3` / sqrt(sta3$`33`)
CV4 = sta3$`4` / sqrt(sta3$`44`)
sta3 = cbind(sta3,CV1,CV2,CV3,CV4)



# Pour les NUTRIMENTS

other_0 = rbind(sta0[,c(1,2,3,6)],sta1[,c(1,2,3,8)],sta2[,c(1,2,3,10)],sta3[,c(1,2,3,12)])

other_0$D = as.factor(other_0$D)
levels(other_0$D) = c(D0,D01)
# other_0$D = factor(other_0$D,levels(other_0$D)[c(2,1)])

other_0$n = as.factor(other_0$n)
levels(other_0$n) = c(T0,T1,T2,T3)

### Texte a ajouter pour TL3, D=0

ann_text_0 <- data.frame(lambda = 0.5 ,CV1 = 10**1.35,lab = "No equilibrium point",
                       D = factor(D0,levels = c(D0,D01)),
                       n = factor(T3, levels = c(T0,T1,T2,T3)))




plot_pierre_0 = ggplot(other_0, aes(x = lambda, y = CV1, color =n))+
  geom_line(size = 1.5)+
  facet_grid(D~., labeller = label_parsed, scales = "free")+
  scale_colour_manual(name="",values=c("black","blue","red","green4","purple"),
                      labels = c(expression(paste("TL0")),expression(paste("TL1")),expression(paste("TL2")),expression(paste("TL3"))) ) +
  xlab(leg_recy_en)+ ylab(leg_inva_en3)+
  geom_text(data = ann_text_0,label = "No equilibrium point", size = 5, show.legend = FALSE)+
  theme_dec + y_axis_log10


 # ggsave("stabi_across_chain_B0_precis_ensemble_changesize.pdf",width = 10, height = 7.69)

# Pour les PLANTES

other_1 = rbind(sta1[,c(1,2,3,9)],sta2[,c(1,2,3,11)],sta3[,c(1,2,3,13)])

other_1$D = as.factor(other_1$D)
levels(other_1$D) = c(D0,D01)
# other_1$D = factor(other_1$D,levels(other_1$D)[c(2,1)])

other_1$n = as.factor(other_1$n)
levels(other_1$n) = c(T1,T2,T3)

### Texte a ajouter pour TL3, D=0

ann_text_1 <- data.frame(lambda = 0.5 ,CV2 = 10**1.5,lab = "No equilibrium point",
                       D = factor(D0,levels = c(D0,D01)),
                       n = factor(T3, levels = c(T1,T2,T3) )) #,
                              #    show.legend = NA))




plot_pierre_1 = ggplot(other_1, aes(x = lambda, y = CV2, color =n))+
  geom_line(size = 1.5)+
  facet_grid(D~., labeller = label_parsed, scales = "fixed")+
  scale_colour_manual(name="",values=c("blue","red","green4","purple"),
                      labels = c(expression(paste("TL1")),expression(paste("TL2")),expression(paste("TL3"))) ) +
  xlab(leg_recy_en)+ ylab(leg_inva_en3)+
  geom_text(data = ann_text_1,label = "No equilibrium point", size = 6, show.legend = FALSE)+
  theme_dec + y_axis_log10

 # ggsave("stabi_across_chain_B1_precis_ensemble_changesize.pdf",width =9, height = 6.92)










