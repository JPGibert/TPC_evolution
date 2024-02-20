library("nls.multstart")
library("dplyr")
library("ggplot2")
library("tidyr")
library("corrplot")
library("reshape2")
library("gridExtra")
library("expm")
library("forcats")
library("tictoc")
library("rstatix")
library("yhat")
library("statgenGxE")
library("inti")
library("ellipse")


# Vignette; https://cran.r-project.org/web/packages/statgenGxE/vignettes/statgenGxE.html
tpcs <- read.csv("~/Desktop/Gibert_lab_project/data/TPC_data.csv")
head(tpcs)
Clones <- unique(tpcs$Clone)

tpcs %>%
	mutate(unique_rep=paste(tpcs$Clone,".",tpcs$Rep,sep=""))

## CALCULATE G, E and GXE
## Crate object structure for statgenGxE
dropsTD <- statgenSTA::createTD(data = tpcs, genotype = "Clone", trial = "Temp")

plot(dropsTD, plotType = "box", traits = "r", colorTrialBy = "genotype",
     orderBy = "ascending")

dropsVarComp <- gxeVarComp(TD = dropsTD, trait = "r")
summary(dropsVarComp) ## With genotype and genotype:temperature as random effects
vc(dropsVarComp)

plot(dropsVarComp)

herit(dropsVarComp)
#The calculation does:
0.57/(0.57+1.2/7+1.11/(6*7))
#which is G/(G+GxE/7+residual/(6*7))


## HERITABILITY From package inti
hr <- H2cal(data = tpcs
            , trait = "r"
            , gen.name = "Clone"
            , rep.n = 6
            , fixed.model = "0 + (1|Temp) + Clone"
            , random.model = "1 + (1|Temp) + (1|Clone)"
            , emmeans = TRUE
            , plot_diag = TRUE
            , outliers.rm = TRUE
            )
hr
## Trait is very heritable! See inti package vignette for details on calculation of heritability
## Standard H^2=0.746, H^2 (Cullis)=0.911, H^2 (Piepho) = 0.953


head(tpcs)
# Preliminary stats (linear model)
mod <- lm(r~Temp*Clone, data=tpcs)
summary(mod)
anova(mod)

# ANOVA for G, E and GXE
# E
tpcs %>%
  group_by(Clone) %>%
  anova_test(r ~ Temp)

# G
tpcs %>%
  group_by(Temp) %>%
  anova_test(r ~ Clone)

# GxE
tpcs %>%
  anova_test(
    r ~ Clone*Temp
  )

###########################################################################################
##### TPC FITS
tpcs <- tpcs %>%
	filter(Final>=1) %>%
	mutate(r = log(Final/Initial)) %>%
	group_by(Clone) %>%
	mutate(r_scale=10,
         log_r=log(r+r_scale)) %>%
  	ungroup


TPC_fits <- tpcs %>%
 	group_by(Clone) %>%
  	do(TPC_fit = nls_multstart(log_r ~ a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temp+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temp+273.15)))),
	data = .,
	iter = 500,
	start_lower = c(a=-10, E_a=0.1, E_d=0.5, Th=285),
	start_upper = c(a=10, E_a=4, E_d=10, Th=330),
	supp_errors = 'Y',
	na.action = na.omit,
	lower = c(a=-10, E_a=0, E_d=0, Th=0))) %>%
  rowwise() %>%
  mutate(a=coef(TPC_fit)[[1]], E_a=coef(TPC_fit)[[2]], E_d=coef(TPC_fit)[[3]], Th=coef(TPC_fit)[[4]], T_opt=E_d*Th/(E_d+8.6e-5*Th*log(E_d/E_a-1)))
  
pars <- dplyr::select(TPC_fits, -TPC_fit)
  
TPC_eqn<-function(a, E_a, E_d, Th, Temperature, r_scale){exp(a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temperature+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temperature+273.15)))))-r_scale}
       

TPC_predicted<-expand.grid(Clone=Clones, Temperature=seq(0, 50, length.out=500)) %>%
  left_join(dplyr::select(TPC_fits, -TPC_fit)) %>%
  left_join(distinct(dplyr::select(tpcs, Clone, r_scale))) %>%
  mutate(r=TPC_eqn(a, E_a, E_d, Th, Temperature, r_scale)) #%>%
  #filter(r>=-1)

TPC_summary_spread<-TPC_predicted %>%
  group_by(Clone) %>%
  mutate(CT_min=ifelse(lag(r)<0 & r>0, Temperature, NA),
         CT_max=ifelse(lag(r)>0 & r<0, Temperature, NA),
         r_peak=max(r)) %>%
  filter(!is.na(CT_min) | !is.na(CT_max)) %>%
  dplyr::select(Clone, CT_min, CT_max, r_peak) %>%
  ungroup() %>%
  gather(param, param_val, -Clone) %>%
  drop_na %>%
  distinct %>%
  spread(param, param_val) %>%
  left_join(dplyr::select(TPC_fits, Clone, E_a, E_d, T_opt)) %>%
  mutate(T_opt=T_opt-273.15, T_range=CT_max-CT_min, TPC_asymmetry=abs((T_opt-CT_min)-(CT_max-T_opt))) %>%
  arrange(r_peak)        

        
# plot all species' TPCs together
TPC_predicted$Clone<-factor(TPC_predicted$Clone)

mean_fitness <- TPC_predicted %>%
	group_by(Temperature) %>%               
    summarize(meanFit=mean(r), SD=sd(r), dissim=mean(dist(r)))                            

g_1<-
ggplot()+
  geom_hline(yintercept=0, color="gray", size=0.7, linetype=2)+
  geom_line(data=TPC_predicted, aes(Temperature, r, color=Clone), size=1)+
  geom_line(data=mean_fitness, aes(Temperature, meanFit), size=1,linetype = "dashed") +
  #scale_color_manual(values=rev(spp_color_palette)) +
  #geom_point(data=tpcs, aes(Temp, r, color=Clone))+
  scale_x_continuous(limits=c(10, 38)) +
    scale_y_continuous(limits=c(-4,12)) +
  labs(x="Temperature (C)", y=Intrinsic~growth~rate~(r)~(cells~cell^-1~d^-1)) +
  theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(0,0,0,1), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="right",
        legend.text=element_text(size=12)
        )

g_2<-ggplot() +
  geom_hline(yintercept=0, color="gray", size=0.7, linetype=2) +
  geom_line(data=TPC_predicted, aes(Temperature, r, color=Clone), size=0.65) +
  geom_line(data=mean_fitness, aes(Temperature, meanFit),linetype = "dashed") +
  geom_point(data=tpcs, aes(Temp, r, color=Clone), size=1.5, shape=1) +
  #geom_jitter(data=tpcs, aes(Temp, r, color=Clone), size=2, shape=1, alpha=1/2) +
  #facet_grid(. ~ Clone) +
  facet_wrap(.~ Clone, ncol = 7) +
  scale_x_continuous(limits=c(10, 38)) +
  scale_y_continuous(limits=c(-4,12)) +
  labs(x="Temperature (ºC)", y=Intrinsic~growth~rate~(r)~(d^-1)) +
  theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=16), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        strip.text = element_text(size=11),
        aspect.ratio=0.7,
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12, face="italic")
        )

#ggsave("~/Desktop/JP/Papers_in_progress/Megan_TPCs/Figures/Fig1.pdf",g_2, width = 30,  height = 30,  units = "cm")
#ggsave("~/Desktop/JP/Papers_in_progress/Megan_TPCs/Figures/Fig2.pdf",g_1, width = 20,  height = 10,  units = "cm")



#####################################################
### GENETIC CORRELATIONS BETWEEN SHAPE PARAMETERS
## Here we only keep parameters E_a, E_d, T_opt, CT_min and E_a
p.mat <- cor.mtest(TPC_summary_spread[,c(-1,-8,-9)][-2,]) #We eliminate B2192.III for this
dev.new()
corrplot(cor(as.data.frame(TPC_summary_spread)[,c(-1,-8,-9)][-2,]),method="color", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat$p, sig.level = 0.05, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE) 

## Only for Ctmin, Ea, Topt, r_peak
p.mat <- cor.mtest(TPC_summary_spread[,c(-1,-2,-6,-8,-9)][-2,]) #We eliminate B2192.III for this
dev.new()
corrplot(cor(as.data.frame(TPC_summary_spread)[,c(-1,-2,-6,-8,-9)][-2,]),method="color", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat$p, sig.level = 0.05, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE) 

#####################################################
### EVOLVABILITY AND G-MATRIX

#"...useful measure of evolvability, because it corresponds to the predicted evolutionary response, in the percentage of the trait mean, to an episode of unit-strength selection (Opedal et al 2022 PNAS)"

## UNIVARIATE CASE
#Create Z-score function
evolvability <- function(array){(sd(array, na.rm=TRUE)^2)/(mean(array, na.rm=TRUE)^2)}

#minus_mean <- function(array){array-mean(array, na.rm=TRUE)^2}

to_bar <- TPC_summary_spread %>% 
	summarise(across(c(2:9),evolvability)) %>%
	gather("Param","evol",1:8) %>%
	mutate(Param=factor(Param,levels=
		c("E_d","E_a","r_peak","TPC_asymmetry","T_range","CT_min","CT_max","T_opt"))) %>%
	arrange(desc(evol))


#pdf("~/Desktop/Gibert_Lab/Grants/NSF/2022/CAREER/Figures/plot.pdf",width = 6,  height = 6)
boxplot(evol~Param, data=to_bar, ylim=c(0,1), las=TRUE)
box(type="l")
dev.off()


## IMPORTANT: Bootstrap clones to reduce effect of bad fits on evolvability calculation
Mat <- TPC_summary_spread[-2,] ## IF clone B2912 III is not eliminated for this analysis, E_a estimates are biased and evolvability analyses are extremely skwewed (mutimodal or disjoint distribution) 
Btstrap_data <- list()
for(i in 1:1000){ ## Commented out are the versions that include all TPC parameters
	re_mat <- Mat[sample(nrow(Mat), replace=TRUE),]
	Btstrap_data[[i]] <- re_mat %>% 
	#summarise(across(c(2:9),evolvability)) %>%
	summarise(across(c(3,4,5,7),evolvability)) %>%
	gather("Param","evol",1:4) %>%
	mutate(Param=factor(Param,levels=
		#c("E_d","r_peak","E_a","TPC_asymmetry","T_range","CT_min","CT_max","T_opt"))) 
		c("r_peak","E_a","CT_min","T_opt")))
}

## Calculate means in list
tomean_dat <- do.call(rbind, Btstrap_data)
mean_dat <- tomean_dat %>%
		group_by(Param) %>%
		summarise(mean_evol=mean(evol,na.rm=TRUE)) 

for(i in 1:1000){
	if(i==1){ 
		boxplot(mean_evol~c(1,2,3,4), data=mean_dat, ylim=c(0,0.2), las=TRUE, col="red", axes=FALSE, ylab="", xlab="")
		box(lwd=2)
		axis(1,c(1,2,3,4), tck=0.015)
		axis(2,c(0,0.05,0.1,0.15,0.2), tck=0.015, las=TRUE)
		# the order (3,1,2,4) here to avoid converting to factor and passing the order 
		points(evol~jitter(c(3,1,2,4)), data=Btstrap_data[[i]], pch=1, col=rgb(0.2,0.2,0.2,0.25))
	}else{
		points(evol~jitter(c(3,1,2,4)), data=Btstrap_data[[i]], pch=1, col=rgb(0.2,0.2,0.2,0.25))
	}	
}


## MULTIVARIATE CASE (see Opedal et al 2022 PNAS)
# Z-Transform the data
	# Here the "scale" function does (x-mean)/sd
Z_TPC_summary <- data.frame(lapply(TPC_summary_spread[,2:7],scale))

#If we want more TPC parameters skip this line
Z_TPC_summary <- Z_TPC_summary %>%
			select(CT_min, E_a, T_opt, r_peak) %>%
			filter(E_a<2 & T_opt>-2)

# GENETIC CORRELATIONS
p.mat <- cor.mtest(Z_TPC_summary) #We eliminate B2192.III for this
dev.new()
corrplot(cor(Z_TPC_summary), method="ellipse", order="FPC", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat$p, sig.level = 0.05, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE) 

corrplot.mixed(cor(Z_TPC_summary),upper="ellipse", order = 'AOE',p.mat = p.mat$p, sig.level = 0.05, insig = "blank")


corM <- cor(Z_TPC_summary)
corM[which(p.mat$p>0.05)]<- NA
plotcorr(cor(Z_TPC_summary))
    
eig <- eigen(cor(Z_TPC_summary))

par(mfrow=c(1,3))
plot(r_peak~E_a, data=Z_TPC_summary, ylim=c(-2,2), xlim=c(-1,0.6))
#abline(lm(CT_min~E_a, data=Z_TPC_summary))

arrows(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$r_peak)
		,-0.575+mean(Z_TPC_summary$E_a),-0.561+mean(Z_TPC_summary$r_peak)
		)  
arrows(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$r_peak)
		,+0.575+mean(Z_TPC_summary$E_a),+0.561+mean(Z_TPC_summary$r_peak)
		)
points(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$r_peak),pch=16)
 		
arrows(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$r_peak)
		,-0.705+mean(Z_TPC_summary$E_a),0.59+mean(Z_TPC_summary$r_peak)
		)  
arrows(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$r_peak)
		,0.705+mean(Z_TPC_summary$E_a),-0.59+mean(Z_TPC_summary$r_peak)
		)  
## Needs to figure out the magnitude of the vectors

plot(CT_min~E_a, data=Z_TPC_summary,ylim=c(-3,3), xlim=c(-1,0.6))
#abline(lm(CT_min~E_a, data=Z_TPC_summary))

arrows(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$CT_min)
		,-0.575+mean(Z_TPC_summary$E_a),-0.505+mean(Z_TPC_summary$CT_min)
		)  
arrows(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$CT_min)
		,+0.575+mean(Z_TPC_summary$E_a),+0.505+mean(Z_TPC_summary$CT_min)
		)
points(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$CT_min),pch=16)
 		
arrows(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$CT_min)
		,-0.705+mean(Z_TPC_summary$E_a),0.29+mean(Z_TPC_summary$CT_min)
		)  
arrows(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$CT_min)
		,0.705+mean(Z_TPC_summary$E_a),-0.29+mean(Z_TPC_summary$CT_min)
		) 
### PUT CT_min as the x axis instead!
plot(CT_min~T_opt, data=Z_TPC_summary,ylim=c(-3,3), xlim=c(-1.5,1.5))
#abline(lm(CT_min~E_a, data=Z_TPC_summary))

arrows(mean(Z_TPC_summary$T_opt),mean(Z_TPC_summary$CT_min)
		,-0.505+mean(Z_TPC_summary$T_opt),-0.31+mean(Z_TPC_summary$CT_min)
		)  
arrows(mean(Z_TPC_summary$T_opt),mean(Z_TPC_summary$CT_min)
		,0.505+mean(Z_TPC_summary$T_opt),0.31+mean(Z_TPC_summary$CT_min)
		)
points(mean(Z_TPC_summary$T_opt),mean(Z_TPC_summary$CT_min),pch=16)
 		
arrows(mean(Z_TPC_summary$T_opt),mean(Z_TPC_summary$CT_min)
		,0.29+mean(Z_TPC_summary$T_opt),-0.24+mean(Z_TPC_summary$CT_min)
		)  
arrows(mean(Z_TPC_summary$T_opt),mean(Z_TPC_summary$CT_min)
		,-0.29+mean(Z_TPC_summary$T_opt),0.24+mean(Z_TPC_summary$CT_min)
		)

arrows(mean(Z_TPC_summary$E_a),mean(Z_TPC_summary$CT_min),0.399,-0.40)
arrows(0,0,-0.575,-0.505)  
arrows(0,0,0.399,-0.40)

         
ggplot(TPC_summary_spread, aes(x=CT_min, y=T_opt)) +
  geom_point() +
  geom_smooth(method=lm)         
lm(CT_min~T_opt,data=TPC_summary_spread)%>%summary()
           


#########################################################
## QUANTIFYING SELECTION

# filter out only necessary temperatures 
TPC_predicted_a<-expand.grid(Clone=Clones, Temperature= c(13.5000000,16.5000000,19.5000000, 22.0000000,25.0000000,28.0000000,31.0000000,34.0000000,37.0000000)) %>%
  left_join(dplyr::select(TPC_fits, -TPC_fit)) %>%
  left_join(distinct(dplyr::select(tpcs, Clone, r_scale))) %>%
  mutate(r=TPC_eqn(a, E_a, E_d, Th, Temperature, r_scale)) %>%
#filter(r>=-1)
  mutate(Temp_range = ifelse(Temperature < 20, "Low", ifelse(Temperature > 30, "High", "Medium")))

# TPC_summary_spread 
# TPC_predicted 1 are extracted value

TPC_predicted_a_LOW <- TPC_predicted_a %>% 
        filter(Temp_range == "Low")
TPC_predicted_a_MED <- TPC_predicted_a %>% 
  filter(Temp_range == "Medium")
TPC_predicted_a_HIGH <- TPC_predicted_a %>% 
  filter(Temp_range == "High")

Test <- merge(TPC_predicted_a,TPC_summary_spread,by="Clone")

# T_opt graph
TPC_predicted_mod <- TPC_predicted_a %>% 
						filter(T_opt>306) 
						
T_opt_mod <- lm(formula = r ~ (T_opt+I(T_opt^2))*Temp_range , data = TPC_predicted_mod)
summary(T_opt_mod)

t_opt_df_low <- data.frame(T_opt=seq(306.5,310.5,length.out=200), Temp_range=rep("Low",200))
t_opt_df_med <- data.frame(T_opt=seq(306.5,310.5,length.out=200), Temp_range=rep("Medium",200))
t_opt_df_high <- data.frame(T_opt=seq(306.5,310.5,length.out=200), Temp_range=rep("High",200))
pred_low <- predict(T_opt_mod,t_opt_df_low,se=TRUE)
pred_med <- predict(T_opt_mod,t_opt_df_med,se=TRUE)
pred_high <- predict(T_opt_mod,t_opt_df_high,se=TRUE)
pred_2_low<- cbind(t_opt_df_low,pred_low$fit)
pred_2_med<- cbind(t_opt_df_med,pred_med$fit)
pred_2_high<- cbind(t_opt_df_high,pred_high$fit)
colnames(pred_2_low) <- c("T_opt", "Temp_range","fit")
colnames(pred_2_med) <- c("T_opt", "Temp_range","fit")
colnames(pred_2_high) <- c("T_opt", "Temp_range","fit")

T_opt_graph <-plot(r~T_opt,data=TPC_predicted_mod,col=c("red","blue","green")[as.factor(TPC_predicted_a$Temp_range)], pch=16, xlim=c(306,311))
lines(fit~T_opt,data=pred_2_low, col="blue", lwd=2)
lines(fit~T_opt,data=pred_2_med, col="green", lwd=2)
lines(fit~T_opt,data=pred_2_high, col="red", lwd=2)

# r_peak graph

TPC_predicted_mod <- merge(TPC_predicted_a,TPC_summary_spread,by="Clone")

summary(TPC_predicted_mod$r_peak)
						
T_opt_mod <- lm(formula = r ~ (r_peak+I(r_peak^2))*Temp_range , data = TPC_predicted_mod)
T_opt_mod <- lm(formula = r ~ (r_peak)*Temp_range , data = TPC_predicted_mod)
summary(T_opt_mod)

r_peak_df_low <- data.frame(r_peak=seq(3.205,11.646,length.out=200), Temp_range=rep("Low",200))
r_peak_df_med <- data.frame(r_peak=seq(3.205,11.646,length.out=200), Temp_range=rep("Medium",200))
r_peak_df_high <- data.frame(r_peak=seq(3.205,11.646,length.out=200), Temp_range=rep("High",200))
pred_low <- predict(T_opt_mod,r_peak_df_low,se=TRUE)
pred_med <- predict(T_opt_mod,r_peak_df_med,se=TRUE)
pred_high <- predict(T_opt_mod,r_peak_df_high,se=TRUE)
pred_2_low<- cbind(r_peak_df_low,pred_low$fit)
pred_2_med<- cbind(r_peak_df_med,pred_med$fit)
pred_2_high<- cbind(r_peak_df_high,pred_high$fit)
colnames(pred_2_low) <- c("r_peak", "Temp_range","fit")
colnames(pred_2_med) <- c("r_peak", "Temp_range","fit")
colnames(pred_2_high) <- c("r_peak", "Temp_range","fit")

plot(r~r_peak,data=TPC_predicted_mod,col=c("red","blue","green")[as.factor(TPC_predicted_mod$Temp_range)], pch=16, xlim=c(3,12))
lines(fit~r_peak,data=pred_2_low, col="blue")
lines(fit~r_peak,data=pred_2_med, col="green")
lines(fit~r_peak,data=pred_2_high, col="red")


# E_a graph

TPC_predicted_mod_ea <- TPC_predicted_a %>% 
						filter(E_a<0.4)
						
## Quadratic term non-significant so not needed
E_a_mod <- lm(formula = r ~(E_a+I(E_a^2))*Temp_range, data = TPC_predicted_mod_ea)
summary(E_a_mod)
## We keep only linear terms
E_a_mod <- lm(formula = r ~(E_a)*Temp_range, data = TPC_predicted_mod_ea)
summary(E_a_mod)

E_a_df_low <- data.frame(E_a=seq(0,1.25,length.out=200), Temp_range=rep("Low",200))
E_a_df_med <- data.frame(E_a=seq(0,1.25,length.out=200), Temp_range=rep("Medium",200))
E_a_df_high <- data.frame(E_a=seq(0,1.25,length.out=200), Temp_range=rep("High",200))
pred_low_E_a <- predict(E_a_mod,E_a_df_low,se=TRUE)
pred_med_E_a <- predict(E_a_mod,E_a_df_med,se=TRUE)
pred_high_E_a <- predict(E_a_mod,E_a_df_high,se=TRUE)
pred_2_low_E_a <- cbind(E_a_df_low,pred_low_E_a$fit)
pred_2_med_E_a <- cbind(E_a_df_med,pred_med_E_a$fit)
pred_2_high_E_a <- cbind(E_a_df_high,pred_high_E_a$fit)
colnames(pred_2_low_E_a) <- c("E_a", "Temp_range","fit")
colnames(pred_2_med_E_a) <- c("E_a", "Temp_range","fit")
colnames(pred_2_high_E_a) <- c("E_a", "Temp_range","fit")

# Could be that there is indeed something more complicated going on, most likely than not, the data is crazy above 0.5
plot(r~E_a,data=TPC_predicted_mod_ea,col=c("red","blue","green")[as.factor(TPC_predicted_mod_ea$Temp_range)], pch=16, xlim=c(0.05,0.5))
lines(fit~E_a,data=pred_2_med_E_a, col="green")
lines(fit~E_a,data=pred_2_high_E_a, col="red")
lines(fit~E_a,data=pred_2_low_E_a, col="blue")


## CTmin graph

TPC_predicted_mod <- merge(TPC_predicted_a,TPC_summary_spread,by="Clone")

summary(TPC_predicted_mod$CT_min)
						
CT_min_mod <- lm(formula = r ~ (CT_min+I(CT_min^2))*Temp_range , data = TPC_predicted_mod)
CT_min_mod <- lm(formula = r ~ (CT_min)*Temp_range , data = TPC_predicted_mod)
summary(CT_min_mod)

CT_min_df_low <- data.frame(CT_min=seq(12.32,21.98,length.out=200), Temp_range=rep("Low",200))
CT_min_df_med <- data.frame(CT_min=seq(12.32,21.98,length.out=200), Temp_range=rep("Medium",200))
CT_min_df_high <- data.frame(CT_min=seq(12.32,21.98,length.out=200), Temp_range=rep("High",200))
pred_low <- predict(CT_min_mod,CT_min_df_low,se=TRUE)
pred_med <- predict(CT_min_mod,CT_min_df_med,se=TRUE)
pred_high <- predict(CT_min_mod,CT_min_df_high,se=TRUE)
pred_2_low<- cbind(CT_min_df_low,pred_low$fit)
pred_2_med<- cbind(CT_min_df_med,pred_med$fit)
pred_2_high<- cbind(CT_min_df_high,pred_high$fit)
colnames(pred_2_low) <- c("CT_min", "Temp_range","fit")
colnames(pred_2_med) <- c("CT_min", "Temp_range","fit")
colnames(pred_2_high) <- c("CT_min", "Temp_range","fit")

CT_min_graph <-plot(r~CT_min,data=TPC_predicted_mod,col=c("red","blue","green")[as.factor(TPC_predicted_mod$Temp_range)], pch=16, xlim=c(12,22))
lines(fit~CT_min,data=pred_2_low, col="blue")
lines(fit~CT_min,data=pred_2_med, col="green")
lines(fit~CT_min,data=pred_2_high, col="red")



# E_d graph (Not great parameter fits so likely not needed)

TPC_predicted_mod_ed <- TPC_predicted_a %>% 
						filter(E_d<10)
## Quadratic term not needed
E_d_mod <- lm(formula = r ~(E_d+I(E_d^2))*Temp_range, data = TPC_predicted_mod_ed)
summary(E_d_mod)
## We keep the linear term
E_d_mod <- lm(formula = r ~(E_d)*Temp_range, data = TPC_predicted_mod_ed)
summary(E_d_mod)


E_d_df_low <- data.frame(E_d=seq(0,30,length.out=200), Temp_range=rep("Low",200))
E_d_df_med <- data.frame(E_d=seq(0,30,length.out=200), Temp_range=rep("Medium",200))
E_d_df_high <- data.frame(E_d=seq(0,30,length.out=200), Temp_range=rep("High",200))
pred_low_E_d <- predict(E_d_mod,E_d_df_low,se=TRUE)
pred_med_E_d <- predict(E_d_mod,E_d_df_med,se=TRUE)
pred_high_E_d <- predict(E_d_mod,E_d_df_high,se=TRUE)
pred_2_low_E_d <- cbind(E_d_df_low,pred_low_E_d$fit)
pred_2_med_E_d <- cbind(E_d_df_med,pred_med_E_d$fit)
pred_2_high_E_d<- cbind(E_d_df_high,pred_high_E_d$fit)
colnames(pred_2_low_E_d) <- c("E_d", "Temp_range","fit")
colnames(pred_2_med_E_d) <- c("E_d", "Temp_range","fit")
colnames(pred_2_high_E_d) <- c("E_d", "Temp_range","fit")

plot(r~E_d,data=TPC_predicted_mod_ed,col=c("red","blue","green")[as.factor(TPC_predicted_mod_ed$Temp_range)], pch=16, xlim=c(0,10))
lines(fit~E_d,data=pred_2_low_E_d, col="blue")
lines(fit~E_d,data=pred_2_med_E_d, col="green")
lines(fit~E_d,data=pred_2_high_E_d, col="red")


### FIGURE 3 (QUANTIFYING SELECTION)

TPC_predicted_mod <- merge(TPC_predicted_a,TPC_summary_spread,by="Clone")
Sel_fig <- TPC_predicted_mod %>% 
						filter(E_a.x<0.4 & T_opt.x>306) %>% # Lose 18 data points
						select(Temp_range, r, E_a.x, T_opt.x, CT_min, r_peak) %>%
						mutate(E_a=E_a.x, T_opt=T_opt.x, .keep = "unused") %>%
						relocate(E_a, r_peak, CT_min,T_opt, .after=r) 
Sel_fig <- cbind(Sel_fig[,-c(3:6)], data.frame(lapply(Sel_fig[,3:6],scale)))						

Sel_fig_gathered <- Sel_fig %>%
  	pivot_longer(cols= E_a:T_opt, names_to='variables') %>%
  mutate(variables=factor(variables, levels=c("r_peak", "E_a", "CT_min", "T_opt")))

## Models
T_opt_mod <- lm(formula = r ~ (T_opt+I(T_opt^2))*Temp_range , data = Sel_fig)
E_a_mod <- lm(formula = r ~ (E_a+I(E_a^2))*Temp_range , data = Sel_fig)
CT_min_mod <- lm(formula = r ~ (CT_min+I(CT_min^2))*Temp_range , data = Sel_fig)
r_peak_mod <- lm(formula = r ~ (r_peak+I(r_peak^2))*Temp_range , data = Sel_fig)
summary(T_opt_mod)
summary(E_a_mod)
summary(CT_min_mod)
summary(r_peak_mod)


pred_func<-function(Var, Temp){
  x_vals<-seq(min(Sel_fig_gathered$value), max(Sel_fig_gathered$value), length.out=200)
  temp_mod<-eval(parse(text=paste0(Var, "_mod")))
  temp_df<-data.frame(Variable=x_vals, Temp_range=Temp)
  colnames(temp_df)[1]<-Var
  temp_predict<-predict(temp_mod, temp_df, se=T)
  return(data.frame(Variable=x_vals, Value=temp_predict$fit, variables=Var, Temp_range=Temp))
}
pred_func("T_opt", "Low")

pred_df<-bind_rows(pred_func("T_opt", "Low"),
                   pred_func("T_opt", "Medium"),
                   pred_func("T_opt", "High"),
                   pred_func("CT_min", "Low"),
                   pred_func("CT_min", "Medium"),
                   pred_func("CT_min", "High"),
                   pred_func("E_a", "Low"),
                   pred_func("E_a", "Medium"),
                   pred_func("E_a", "High"),
                   pred_func("r_peak", "Low"),
                   pred_func("r_peak", "Medium"),
                   pred_func("r_peak", "High")) %>% 
  mutate(variables=factor(variables, levels=c("r_peak", "E_a", "CT_min", "T_opt")))
  
  
## Figure (NEEDS TO OVERLAY MODEL PREDICTIONS STILL)
Fig2<-ggplot() +
  geom_hline(yintercept=0, color="gray", size=0.7, linetype=2) +
  #geom_line(data=TPC_predicted, aes(Temperature, r, color=Clone), size=0.65) +
  #geom_line(data=mean_fitness, aes(Temperature, meanFit),linetype = "dashed") +
  geom_point(data=Sel_fig_gathered, aes(x=value, y=r, color=Temp_range), size=1.5, shape=16) +
  
  geom_line(data=pred_df, aes(Variable, Value, color=Temp_range)) +
  
  scale_color_manual(values=c("orangered", "skyblue", "orange")) +
  #geom_jitter(data=tpcs, aes(Temp, r, color=Clone), size=2, shape=1, alpha=1/2) +
  facet_wrap(. ~ variables, scales = "free", ncol=2) +
  #scale_x_continuous(limits=c(10, 38)) +
  scale_y_continuous(limits=c(-6,13)) +
  labs(x="Variable", y=Intrinsic~growth~rate~(r)~(d^-1)) +
  theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=16), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        strip.text = element_text(size=11),
        aspect.ratio=0.7,
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12, face="italic")
        ) 
  







     
	# Is variability in r increasing with temperature?         
SD_plot <- ggplot(mean_fitness, aes(x=Temperature, y=SD)) +
  	geom_line() +
  	geom_smooth(method=lm) + 	
  	scale_x_continuous(limits=c(10, 38), sec.axis = dup_axis()) +
  	theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.05, "cm"),
        plot.margin=unit(c(0,0,0,0), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12, face="italic")
        )
## Increased variability in r across genotypes with temperature that is roughly linear with temp-

SD_plot <- ggplot(mean_fitness, aes(x=Temperature, y=SD)) +
  	geom_line() +
  	geom_smooth(method=lm) + 	
  	scale_x_continuous(limits=c(10, 38)) +
  	scale_y_continuous(limits=c(0,2.5)) +
  	theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(0,0,0,0), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12, face="italic")
        )

Diss_plot <- ggplot(mean_fitness, aes(x=Temperature, y=dissim)) +
  	geom_line() +
  	geom_smooth(method=lm) + 	
  	scale_x_continuous(limits=c(10, 38)) +
  	scale_y_continuous(limits=c(0,2.5)) +
  	theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(0,0,0,0), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12, face="italic")
        )

grid.arrange(SD_plot,Diss_plot)

## Quantify dominance of each genotype (likelihood they will persist across environmental conditions)
## Calculate selection relative ro Clone A*V
            
TPC_forselection_1 <- TPC_predicted %>%
	filter(Clone=="A*V") %>%
	select(Temperature,r) %>%
	mutate(r_star=r, r=NULL)

TPC_forselection <- merge(TPC_predicted,TPC_forselection_1,by="Temperature") %>%
	group_by(Clone) %>%
	mutate(Diff=(r+4)/(r_star+4), sel=1-Diff)	
	
g_3<-
ggplot()+
  geom_hline(yintercept=1, color="gray", size=0.7, linetype=2)+
  geom_line(data=TPC_forselection, aes(Temperature, Diff, color=Clone), size=1)+
  #scale_color_manual(values=rev(spp_color_palette)) +
  #geom_point(data=tpcs, aes(Temp, r, color=Clone))+
  scale_x_continuous(limits=c(10, 37)) +
  scale_y_continuous(limits=c(0.25,1.5))+
  labs(x="Temperature (C)", y="Relative Fitness") +
  theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(0,0,0,1), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12)
        )
g_4<-
ggplot()+
  geom_hline(yintercept=0, color="gray", size=0.7, linetype=2)+
  geom_line(data=TPC_forselection, aes(Temperature, sel, color=Clone), size=1)+
  #scale_color_manual(values=rev(spp_color_palette)) +
  #geom_point(data=tpcs, aes(Temp, r, color=Clone))+
  scale_x_continuous(limits=c(15, 37)) +
  scale_y_continuous(limits=c(-1,1))+
  labs(x="Temperature (C)", y="Selection coefficient") +
  theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(0,0,0,1), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12)
        )


lay=rbind(c(2,2,1,1),c(2,2,3,3),c(2,2,4,4))

grid.arrange(g_1,g_2,SD_plot, g_4, layout_matrix=lay)


# Plot descriptotrs of the TPCS

g_6 <- TPC_summary_spread %>%
  ggplot(aes(x=reorder(Clone,T_opt),y=T_opt)) +
  geom_hline(yintercept=30, color="gray", size=0.7, linetype=2)+
  geom_point(size = 5, aes(color=Clone)) + 
  geom_segment(aes(xend = Clone, yend =30,color=Clone), size = 1.2)+
  #scale_x_continuous(limits=c(15, 37)) +
  scale_y_continuous(limits=c(30,38)) +
  #facet_wrap(Parameter~., ncol = 2) +
  #labs(x="Temperature (C)", y="Selection coefficient") +
  theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(0,0,0,1), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="Right",
        legend.text=element_text(size=12)
        )                

g_7 <- TPC_summary_spread %>%
  ggplot(aes(x=reorder(Clone,T_opt),y=E_a)) +
  geom_hline(yintercept=0.1, color="gray", size=0.7, linetype=2)+
  geom_point(size = 5, aes(color=Clone)) + 
  geom_segment(aes(xend = Clone, yend = 0.1,color=Clone), size = 1.2)+
  #scale_x_continuous(limits=c(15, 37)) +
  scale_y_continuous(limits=c(0.1,0.3)) +
  #facet_wrap(Parameter~., ncol = 2) +
  #labs(x="Temperature (C)", y="Selection coefficient") +
  theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(0,0,0,1), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="Right",
        legend.text=element_text(size=12)
        ) 

g_8 <- TPC_summary_spread %>%
  ggplot(aes(x=reorder(Clone,T_opt),y=r_peak)) +
  geom_hline(yintercept=3, color="gray", size=0.7, linetype=2)+
  geom_point(size = 5, aes(color=Clone)) + 
  geom_segment(aes(xend = Clone, yend = 3,color=Clone), size = 1.2)+
  #scale_x_continuous(limits=c(15, 37)) +
  scale_y_continuous(limits=c(3,8)) +
  #facet_wrap(Parameter~., ncol = 2) +
  #labs(x="Temperature (C)", y="Selection coefficient") +
  theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(0,0,0,1), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="Right",
        legend.text=element_text(size=12)
        ) 

g_9 <- TPC_summary_spread %>%
  ggplot(aes(x=reorder(Clone,T_opt),y=CT_min)) +
  geom_hline(yintercept=10, color="gray", size=0.7, linetype=2)+
  geom_point(size = 5, aes(color=Clone)) + 
  geom_segment(aes(xend = Clone, yend = 10,color=Clone), size = 1.2)+
  #scale_x_continuous(limits=c(15, 37)) +
  scale_y_continuous(limits=c(10,20)) +
  #facet_wrap(Parameter~., ncol = 2) +
  #labs(x="Temperature (C)", y="Selection coefficient") +
  theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=14), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(0,0,0,1), "cm"),
        aspect.ratio=0.55,
        legend.key=element_blank(),
        legend.position="Right",
        legend.text=element_text(size=12)
        ) 

dev.new()
plo2 <- grid.arrange(g_6,g_7,g_8,g_9)
ggsave("~/Desktop/Gibert_Lab/Grants/NSF/2022/CAREER/Figures/plot2.pdf",plo2, width = 23,  height = 14,  units = "cm")

####
## Model predictions of genetic frequencies with temperature
###. 

head(TPC_predicted)

TPC_predicted <- merge(TPC_predicted,mean_fitness,by="Temperature")
TPC_predicted$rel_fitness <- (TPC_predicted$r+10)/(TPC_predicted$meanFit+10)
Temp <- unique(TPC_predicted$Temperature)

fit <- TPC_predicted %>%
	filter(Temperature==Temp[250]) %>%
	select(r)

## Add clonal_lines
Clone_order <- TPC_predicted %>%
	filter(Temperature==Temp[250]) %>%
	select(Clone)

## Define function for model
## Temperature has to be passed as Temp[nmbr]
evo_mod <- function(temp,noise,time_steps){
	## Select appropriate temperature data
	fit <- TPC_predicted %>%
		filter(Temperature==temp) %>%
		select(r)

	## Run model
	N_clones <- 16
	time <- seq(1,time_steps,1)
	mat_freq <- matrix(rep(0,N_clones*time_steps),nrow=N_clones,ncol=time_steps)
	mat_freq[,1] <- rep(1/N_clones,N_clones)	
	for(i in 1:(time_steps-1)){
		if(i==1){ #We add +4 to r to avoid dealing with negative vaues of r and how it impacts the calculation of rel.fit.
			new_freq <- pmax(diag((fit$r+4)/(mean(fit$r)+4),N_clones,N_clones)%*%rep(1/N_clones,N_clones) + runif(N_clones,min=-noise,max=noise),0) ## pmax turn 					negative frequencies into 0s
			mat_freq[,2] <- new_freq/sum(new_freq) # Dividing my the sum ensures that sum(frequencies)=1. It rescales all 					frequencies so that losing alleles doesn't lower the sum(frequencies). Needed due to randomness 
		}else{
			mat_freq[,i+1] <-pmax(diag((fit$r+4)/sum(mat_freq[,i-1]*(fit$r+4)),N_clones,N_clones)%*%mat_freq[,i]+ runif(12,min=-noise,max=noise),0)
			mat_freq[,i+1] <- mat_freq[,i+1]/sum(mat_freq[,i+1])
		}
	}
	return(rowMeans(mat_freq[,1:500]))
	#return(mat_freq[,1:500]) 
	# If change for return(mat_freq) it ispossible to plot the dynamics of the model against time
}

## Try
mat <- evo_mod(Temp[200],0.05,500)

## This is good to see the dynamics of the model against time
matplot(time,t(mat[,1:200]), type = 'l', lty=1) 

## Actual runs
runs <- sapply(Temp,evo_mod,0.02,500)
matplot(Temp,t(runs), type = 'l', lty=1)
## Prep data for plotting

runs_clone <- cbind(Clone_order,runs)
colnames(runs_clone) <- c("Clone",Temp)

new_runs <- runs_clone %>%
	gather("Temp","Freq",2:501) %>%
	mutate(Temp=as.numeric(Temp))

## Plot
g_10 <- ggplot(new_runs,aes(x=Temp, y=Freq, fill=Clone)) + 
    		geom_area(size=0.1, colour="white") +
    		scale_x_continuous(limits=c(12, 38)) +
  			theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        		panel.border=element_rect(color = "black", size=1, fill=NA),
        		axis.text=element_text(size=14), axis.title=element_text(size=16),
		        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
		        axis.ticks.length.x.top=unit(0, "cm"),
		        axis.ticks.length=unit(-0.15, "cm"),
		        plot.margin=unit(c(1,1,1,1), "cm"),
		        aspect.ratio=0.55,
		        legend.key=element_blank(),
		        legend.position="Right",
		        legend.text=element_text(size=12)
		        ) 
plot3 <- grid.arrange(g_1,g_10)

ggsave("~/Desktop/Gibert_Lab/Grants/NSF/2022/CAREER/Figures/plot3.pdf",plot3, width = 14,  height = 23,  units = "cm")

####
## Model predictions of change in genetic frequencies with fluctuating temperatures
###.

# Function to find closest value
closest<-function(array,value){
	#find position
	which(abs(array-value)==min(abs(array-value)))}

evo_mod_Temp_noise <- function(temp,noise,temp_noise,time_steps){
	
	# Generate array of temperatures with noise
	noise_Temp <- Temp[closest(Temp,temp)]+runif(500,-temp_noise,temp_noise)
	# Generate array of fitness at the temp of interest to be used
	temp_ID <- sapply(noise_Temp,closest,Temp)
	
	## Run model
	time <- seq(1,time_steps,1)
	mat_freq <- matrix(rep(0,10*time_steps),nrow=10,ncol=time_steps)
	mat_freq[,1] <- rep(1/10,10)	
	for(i in 1:(time_steps-1)){
			## Select appropriate temperature data
	fit <- TPC_predicted %>%
		filter(Temperature==Temp[temp_ID[i]]) %>%
		select(r)
		
		if(i==1){ 
			#We add +4 to r to avoid dealing with negative vaues of r and how it impacts the calculation of rel.fit.
			new_freq <- pmax(diag((fit$r+4)/(mean(fit$r)+4),10,10)%*%rep(1/10,10) + runif(10,min=-noise,max=noise),0) ## pmax turn 					negative frequencies into 0s
			mat_freq[,2] <- new_freq/sum(new_freq) # Dividing my the sum ensures that sum(frequencies)=1. It rescales all 					frequencies so that losing alleles doesn't lower the sum(frequencies). Needed due to randomness 
		}else{
			mat_freq[,i+1] <-pmax(diag((fit$r+4)/sum(mat_freq[,i-1]*(fit$r+4)),10,10)%*%mat_freq[,i]+ runif(10,min=-noise,max=noise),0)
			mat_freq[,i+1] <- mat_freq[,i+1]/sum(mat_freq[,i+1])
		}
	}
	return(rowMeans(mat_freq[,1:time_steps]))
	#return(mat_freq[,1:time_steps]) 
	# If change for return(mat_freq) it ispossible to plot the dynamics of the model against time
}

runs <- evo_mod_Temp_noise(32,0.01,5,300)
matplot(time[1:300],t(runs[,1:300]), type = 'l', lty=1)

## Actual runs
tic()
runs <- sapply(Temp,evo_mod_Temp_noise,0.01,5,300)
toc()

runs_clone <- cbind(Clone_order,runs)
colnames(runs_clone) <- c("Clone",Temp)

new_runs_2 <- runs_clone %>%
	gather("Temp","Freq",2:501) %>%
	mutate(Temp=as.numeric(Temp))

## Plot
g_11 <- ggplot(new_runs_2,aes(x=Temp, y=Freq, fill=Clone)) + 
    		geom_area(size=0.1, colour="white") +
    		scale_x_continuous(limits=c(10, 38)) +
  			theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        		panel.border=element_rect(color = "black", size=1, fill=NA),
        		axis.text=element_text(size=14), axis.title=element_text(size=16),
		        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
		        axis.ticks.length.x.top=unit(0, "cm"),
		        axis.ticks.length=unit(-0.15, "cm"),
		        plot.margin=unit(c(1,1,1,1), "cm"),
		        aspect.ratio=0.55,
		        legend.key=element_blank(),
		        legend.position="Right",
		        legend.text=element_text(size=12)
		        )
g_2		        
dev.new()
grid.arrange(g_10,g_11)




         # THE END