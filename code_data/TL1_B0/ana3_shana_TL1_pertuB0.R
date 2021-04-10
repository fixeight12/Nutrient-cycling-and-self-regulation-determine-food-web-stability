setwd("~/MASTER Sciences de la Mer/Stage M2/projets/david_shanafelt3.0/analyse_R/TL1_B0")


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
leg_in_en = TeX("Stability")

N1_s_en = TeX(("Mineral nutrient $\\B_0$"))
N2_s_en = TeX(("Primary producer $\\B_1$"))
N3_s_en = TeX(("Herbivore $\\B_2$"))

N1_a = TeX("$\\B_0$")
N2_a = TeX("$\\B_1$")




y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))

x_axis_log10<-scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))

y_axis_log10<-scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))




######RECYCLAGE COMPENSé #################


# Graphique avec Inva/Bio/Var

datcom = read.table(file = "data_recy_TL1_pertuB0.txt" , fill = TRUE, sep = ";", header = TRUE)
datcom = datcom[,-54]


# Grapgique avec Inva/Bio/Var

all_recycom = datcom[,c(11,12,23:24,31,34)]
all_recycom = all_recycom[which((all_recycom$D == 0)|(all_recycom$D==1)),]

colnames(all_recycom) = c("lambda","D","1","2","11","22")


CV1 = all_recycom$`1` / sqrt(all_recycom$`11`)
CV2 = all_recycom$`2` / sqrt(all_recycom$`22`)


all_recycom = cbind(all_recycom,CV1,CV2)
all_recycom$`11` = sqrt(all_recycom$`11`)
all_recycom$`22` = sqrt(all_recycom$`22`)

comparti = rep(c(rep(1,length(all_recycom$lambda)),rep(2,length(all_recycom$lambda))),3)
forme = c(rep(1,length(all_recycom$lambda)*2),rep(2,length(all_recycom$lambda)*2),rep(3,length(all_recycom$lambda)*2))



all_recycom = melt(all_recycom,id.vars = c("lambda", "D"))
all_recycom = cbind(all_recycom,comparti,forme)


####

all_recycom$D = as.factor(all_recycom$D)
levels(all_recycom$D)<-c(D0,D1)

all_recycom$forme = as.factor(all_recycom$forme)
all_recycom$forme = factor(all_recycom$forme,levels(all_recycom$forme)[c(3,2,1)])
all_recycom$comparti = as.factor(all_recycom$comparti)
levels(all_recycom$comparti) = c(N1_a,N2_a)




plot_all_recycom = ggplot(data = all_recycom,aes( x = lambda , y = value,color = forme, linetype = forme))+
  geom_line(size = 1.5)+
  
 scale_linetype_discrete(name = "", 
                          labels = c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))) ) +
  
  facet_grid(D~comparti, labeller = label_parsed,scales = "free")+
  scale_colour_manual(name="",values=c("blue","peru","slategrey"),
                      labels = c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))) ) + # Set legend title
  xlab(leg_recy_en) + ylab(leg_in_en) + # Set axis labels
  theme_dec + y_axis_log10_short + scale_x_continuous(breaks = c(0,0.5,0.9))

# ggsave("all_recy_TL1_pertuB1_comp.pdf",width = 13, height = 10)


####### SELF-REGULATION COMPENSé##########


dacom = read.table(file = "data_D_TL1_pertuB0.txt" , fill = TRUE, sep = ";", header = TRUE)
dacom = dacom[,-54]




# Graphique Inva/Bio/Var


all_Dcom = dacom[,c(11,12,23:24,31,34)]
all_Dcom = all_Dcom[which((all_Dcom$lambda == 0)|(all_Dcom$lambda == 0.9)),]
colnames(all_Dcom) = c("lambda","D","1","2","11","22")


CV1 = all_Dcom$`1` / sqrt(all_Dcom$`11`)
CV2 = all_Dcom$`2` / sqrt(all_Dcom$`22`)

all_Dcom = cbind(all_Dcom,CV1,CV2)
all_Dcom$`11` = sqrt(all_Dcom$`11`)
all_Dcom$`22` = sqrt(all_Dcom$`22`)

comparti = rep(c(rep(1,length(all_Dcom$lambda)),rep(2,length(all_Dcom$lambda))),3)
forme = c(rep(1,length(all_Dcom$lambda)*2),rep(2,length(all_Dcom$lambda)*2),rep(3,length(all_Dcom$lambda)*2))



all_Dcom = melt(all_Dcom,id.vars = c("lambda", "D"))

all_Dcom = cbind(all_Dcom, comparti,forme)

all_Dcom$lambda = as.factor(all_Dcom$lambda)
levels(all_Dcom$lambda)<-c(lam0,lam09)
all_Dcom$lambda = factor(all_Dcom$lambda,levels(all_Dcom$lambda)[c(3,2,1)])

all_Dcom$forme = as.factor(all_Dcom$forme)
all_Dcom$forme = factor(all_Dcom$forme,levels(all_Dcom$forme)[c(3,2,1)])
all_Dcom$comparti = as.factor(all_Dcom$comparti)
levels(all_Dcom$comparti) = c(N1_a,N2_a)


plot_inva_Dcom = ggplot(data = all_Dcom,aes( x = D , y = value, color = forme, linetype =forme))+
  geom_line(size = 1.5)+
  
  scale_linetype_discrete(name = "", 
                          labels = c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))) ) +
  
  facet_grid(lambda~comparti, labeller = label_parsed, scales = "fixed")+
  scale_colour_manual(name="",values=c("blue","peru","slategrey"), 
                      labels =  c(expression(paste("mean biomass", sep = "; "  ,B[i]^"*")),expression(paste("standard deviation  of biomass" , sep = "; "  , sigma[i])),expression(paste("invariability of biomass", sep = "; " ,S[i]))) ) + # Set legend title
  xlab(leg_d_en) + ylab(leg_in_en) + # Set axis labels
  theme + x_axis_log10_short + y_axis_log10_short

# ggsave("all_D_TL1_pertuB1_comp.pdf", width = 13, height = 10)


# Put both on one

all_pertuB0 = plot_grid(plot_all_recycom + theme(legend.position="none")
                       , plot_inva_Dcom + theme(legend.position="none")
                       , labels = c("A", "B"), label_size = 24)



legend <- get_legend(
  # create some space to the left of the legend
  plot_all_recycom + theme(legend.box.margin = margin(0, 0, 0, 7))
)

all_pertuB0 = plot_grid(all_pertuB0, legend, rel_widths = c(1.5, 0.5))

  # ggsave("all_pertuB0.pdf", width = 18,height = 8)



