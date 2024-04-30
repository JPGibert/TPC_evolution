### LOAD PACKAGES----------
install.packages("nls.multstart")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("corrplot")
install.packages("reshape2")
install.packages("gridExtra")
install.packages("expm")
install.packages("forcats")
install.packages("tictoc")
install.packages("rstatix")
install.packages("yhat")
install.packages("statgenGxE")
install.packages("inti")
install.packages("ellipse")
install.packages("ggpubr")
install.packages("gridExtra")
install.packages("grid")
install.packages("ggforce")
install.packages("evolqg")
install.packages("foreach")

library("nls.multstart")
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
library("tidyverse")
library("ggpubr")
library("grid")
library("lme4")
library("ggforce")
library("evolqg")
library("foreach")
library("dplyr")
 
### LOAD AND PROCESS DATA -----------

## Load TPC data from initial 22 genotypes
tpcs <- read.csv("~/Desktop/Gibert_lab_project/data/TPC_data.csv")
Clones <- unique(tpcs$Clone)
tpcs %>%
  mutate(unique_rep=paste(tpcs$Clone,".",tpcs$Rep,sep=""))

## Load TPC data for genotypes AXS and CU4106, used in competition experiment
tpc_exp_pred <- read.csv("~/Desktop/Gibert_lab_project/data/two_clone_TPC.csv") %>% mutate( r_scale=10,log_r=log(r+r_scale)) 
tpc_exp_pred$Clone[tpc_exp_pred$Clone =="CU416"] <- "CU4106"

## Load data from competition experiment between genotypes AXS and CU4106. 
# Non-compiled data available in file "Flow_Cyt_Data.zip"
# YFP indicates individuals that were "gated" or captured as strain AXS by flow cytometry
# AutoF indicates individuals that were "gated" or captured as strain CU4106 by flow cytometry
mar10data <- read.csv("~/Desktop/Gibert_lab_project/data/FlowCytData.csv")
mar10data <- mar10data[-c(1)]
colnames(mar10data) <- c("id", "FSC_H", "SSC_H","FITC_H", "PE_H", "PerCP_H", "APC_H","FSC_A","SSC_A","FITC_A", "PE_A", "PerCP_A",
                         "APC_A","Width","Time", "temp", "ab", "rep","gate", "controlrep")
mar10_NAdata <- read.csv("~/Desktop/Gibert_lab_project/data/Flow_Cyt_Data/NA_Values.csv")
colnames(mar10_NAdata) <-  c("id", "FSC_H", "SSC_H","FITC_H", "PE_H", "PerCP_H", "APC_H","FSC_A","SSC_A","FITC_A", "PE_A", "PerCP_A",
                             "APC_A","Width","Time", "temp", "ab", "rep","gate", "controlrep" )

### FIGURE 2 --------------------------

## Process TPC data from 22 genotypes
# Calculate heritability and G/E/GxE using statgenGxE, vignette available here: https://cran.r-project.org/web/packages/statgenGxE/vignettes/statgenGxE.html
# Create object structure for "statgenGxE" package
dropsTD <- statgenSTA::createTD(data = tpcs, genotype = "Clone", trial = "Temp")

# Visualize r by genotype in ascending order
plot(dropsTD, plotType = "box", traits = "r", colorTrialBy = "genotype",
     orderBy = "ascending")

dropsVarComp <- gxeVarComp(TD = dropsTD, trait = "r")
summary(dropsVarComp) 
# With genotype and genotype:temperature as random effects, display amount of observed variation attributable to each variable
vc(dropsVarComp)
plot(dropsVarComp)
herit(dropsVarComp)

# Calculate heritability using "inti" package
hr <- H2cal(data = tpcs
            , trait = "r"
            , gen.name = "Clone"
            , rep.n = 6
            , fixed.model = "0 + (1|Temp) + Clone"
            , random.model = "1 + (1|Temp) + (1|Clone)"
            , emmeans = TRUE
            , plot_diag = TRUE
            , outliers.rm = TRUE)

# Reveal table with heritabilities
hr$tabsmr %>% table()
# Calculations show trait is very heritable; Standard H^2=0.846, H^2 (Cullis)=0.930, H^2 (Piepho) = 0.971. See inti package vignette for details.

# Create preliminary stats for a linear model
mod <- lm(r~Temp*Clone, data=tpcs)
summary(mod)
anova(mod)

# Run ANOVA for E, G, and GxE
# E
tpcs %>%
  group_by(Clone) %>%
  anova_test(r ~ Temp)
# G
tpcs %>%
  group_by(Temp2=as.factor(Temp)) %>%
  anova_test(r ~ Clone)
# GxE
tpcs %>%
  mutate(Temp2=as.factor(Temp))%>%
  anova_test(
    r ~ Clone*Temp2)

# Generate full data set and calculate log_r from available data
tpcs <- tpcs %>%
  filter(Final>=1) %>%
  mutate(r = log(Final/Initial)) %>%
  dplyr::group_by(Clone) %>%
  mutate(r_scale=10,
         log_r=log(r+r_scale)) %>%
  ungroup

# Use nls.multstart package to fit the curve to TPC data
TPC_fits <- tpcs %>%
  dplyr::group_by(Clone) %>%
  do(TPC_fit = nls_multstart(log_r ~ a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temp+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temp+273.15)))),
                             data = .,
                             iter = 500,
                             start_lower = c(a=-10, E_a=0.1, E_d=0.5, Th=285),
                             start_upper = c(a=10, E_a=4, E_d=10, Th=330),
                             supp_errors = 'Y',
                             na.action = na.omit,
                             lower = c(a=-10, E_a=0, E_d=0, Th=0))) %>%
  rowwise() %>%
  dplyr::mutate(a=coef(TPC_fit)[[1]], 
         E_a=coef(TPC_fit)[[2]], 
         E_d=coef(TPC_fit)[[3]], 
         Th=coef(TPC_fit)[[4]], 
         T_opt=E_d*Th/(E_d+8.6e-5*Th*log(E_d/E_a-1)))

pars <- dplyr::select(TPC_fits, -TPC_fit)

# set TPC equation
TPC_eqn<-function(a, E_a, E_d, Th, Temperature, r_scale){exp(a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temperature+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temperature+273.15)))))-r_scale}

TPC_predicted<-expand.grid(Clone=Clones, Temperature=seq(0, 50, length.out=500)) %>%
  left_join(dplyr::select(TPC_fits, -TPC_fit)) %>%
  left_join(distinct(dplyr::select(tpcs, Clone, r_scale))) %>%
  mutate(r=TPC_eqn(a, E_a, E_d, Th, Temperature, r_scale))

TPC_summary_spread<-TPC_predicted %>%
  dplyr::group_by(Clone) %>%
  dplyr::mutate(CT_min=ifelse(lag(r)<0 & r>0, Temperature, NA),
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
  dplyr::mutate(T_opt=T_opt-273.15, T_range=CT_max-CT_min, TPC_asymmetry=abs((T_opt-CT_min)-(CT_max-T_opt))) %>%
  arrange(r_peak)        

# Plot all species' TPCs together
TPC_predicted$Clone<-factor(TPC_predicted$Clone)

# Calculate mean fitness curve for all genotypes 
mean_fitness <- TPC_predicted %>%
  dplyr::group_by(Temperature) %>%               
  dplyr::summarize(meanFit=mean(r), SD=sd(r), dissim=mean(dist(r)))     

# Figure 2; all TPCs overlaid
g_1<-
  ggplot()+
  geom_hline(yintercept=0, color="gray", linewidth=0.7, linetype=2)+
  geom_line(data=TPC_predicted, aes(Temperature, r, color=Clone), linewidth=1)+
  geom_line(data=mean_fitness, aes(Temperature, meanFit), linewidth=1,linetype = "dashed")+ 
  scale_x_continuous(limits=c(10, 38)) +
  scale_y_continuous(limits=c(-4,12)) +
  labs(x="Temperature (C)", y=Intrinsic~growth~rate~(r)~(cells~cell^-1~d^-1)) +
  theme(plot.background=element_blank(), panel.background=element_blank(),
        panel.border=element_rect(color = "black", linewidth=1, fill=NA),
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

# Alternative to Figure 2, displaying all 22 individual TPCs 
g_2<-ggplot() +
  geom_hline(yintercept=0, color="gray", linewidth=0.7, linetype=2) +
  geom_line(data=TPC_predicted, aes(Temperature, r, color=Clone), linewidth=0.65) +
  geom_line(data=mean_fitness, aes(Temperature, meanFit),linetype = "dashed") +
  geom_point(data=tpcs, aes(Temp, r, color=Clone), size=1.5, shape=1)+
  facet_wrap(.~ Clone, ncol = 7) +
  scale_x_continuous(limits=c(10, 38)) +
  scale_y_continuous(limits=c(-4,12)) +
  labs(x="Temperature (ºC)", y=Intrinsic~growth~rate~(r)~(d^-1)) +
  theme(plot.background=element_blank(), panel.background=element_blank(),
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

### FIGURE 3A, 3B, 3C, 3D ------------------

# GENETIC CORRELATION PLOT BETWEEN SHAPE PARAMETERS
  
# Here we only keep relevant parameters E_a, E_d, T_opt, CT_min
p.mat <- cor.mtest(TPC_summary_spread[,c(-1,-8,-9)][-2,]) #We eliminate genotype B2192.III for this
dev.new()
# generate correlation plot between all parameters
corrplot(cor(as.data.frame(TPC_summary_spread)[,c(-1,-8,-9)][-2,]),method="color", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat$p, sig.level = 0.05, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE) 

# correlation plot for only CT_min, E_a, T_opt, r_peak
p.mat <- cor.mtest(TPC_summary_spread[,c(-1,-2,-6,-8,-9)][-2,]) #We eliminate B2192.III for this
dev.new()
corrplot(cor(as.data.frame(TPC_summary_spread)[,c(-1,-2,-6,-8,-9)][-2,]),method="color", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat$p, sig.level = 0.05, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE) 

# QUANTIFYING SELECTION

# Filter out only necessary temperatures
TPC_predicted_a<-expand.grid(Clone=Clones, Temperature= c(13.5000000,16.5000000,19.5000000, 22.0000000,25.0000000,28.0000000,31.0000000,34.0000000,37.0000000)) %>%
  left_join(dplyr::select(TPC_fits, -TPC_fit)) %>%
  left_join(distinct(dplyr::select(tpcs, Clone, r_scale))) %>%
  dplyr::mutate(r=TPC_eqn(a, E_a, E_d, Th, Temperature, r_scale)) %>%
  #designate temperature bins Low, Medium, High
  dplyr::mutate(Temp_range = ifelse(Temperature < 20, "Low", ifelse(Temperature > 30, "High", "Medium")))

# Split data into three temperature treatments, designated as Low, Medium, and High
TPC_predicted_a_LOW <- TPC_predicted_a %>%
  filter(Temp_range == "Low")
TPC_predicted_a_MED <- TPC_predicted_a %>%
  filter(Temp_range == "Medium")
TPC_predicted_a_HIGH <- TPC_predicted_a %>%
  filter(Temp_range == "High")

TPC_predicted_mod <- merge(TPC_predicted_a,TPC_summary_spread,by="Clone")
Sel_fig <- TPC_predicted_mod %>% 
  filter(E_a.x<0.4 & T_opt.x>306) %>% 
  select(Temp_range, r, E_a.x, T_opt.x, CT_min, r_peak) %>%
  mutate(E_a=E_a.x, T_opt=T_opt.x, .keep = "unused") %>%
  relocate(E_a, r_peak, CT_min,T_opt, .after=r) 
Sel_fig <- cbind(Sel_fig[,-c(3:6)], data.frame(lapply(Sel_fig[,3:6],scale)))						

Sel_fig_gathered <- Sel_fig %>%
  pivot_longer(cols= E_a:T_opt, names_to='variables') %>%
  mutate(variables=factor(variables, levels=c("r_peak", "E_a", "CT_min", "T_opt")))

## Generate models for each parameter of interest
T_opt_mod <- lm(formula = r ~ (T_opt+I(T_opt^2))*Temp_range , data = Sel_fig)
E_a_mod <- lm(formula = r ~ (E_a+I(E_a^2))*Temp_range , data = Sel_fig)
CT_min_mod <- lm(formula = r ~ (CT_min+I(CT_min^2))*Temp_range , data = Sel_fig)
r_peak_mod <- lm(formula = r ~ (r_peak+I(r_peak^2))*Temp_range , data = Sel_fig)

# Predict r for each parameter of interest across temperature interests
pred_func<-function(Var, Temp){
  x_vals<-seq(min(Sel_fig_gathered$value), max(Sel_fig_gathered$value), length.out=200)
  temp_mod<-eval(parse(text=paste0(Var, "_mod")))
  temp_df<-data.frame(Variable=x_vals, Temp_range=Temp)
  colnames(temp_df)[1]<-Var
  temp_predict<-predict(temp_mod, temp_df, se=T)
  return(data.frame(Variable=x_vals, Value=temp_predict$fit, variables=Var, Temp_range=Temp))
}

# Bind all combinations so data can be displayed in one graph
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


# Figure 3, showing GxE framework across thermal performance parameters
Fig3<-ggplot() +
  geom_hline(yintercept=0, color="gray", linewidth=0.7, linetype=2) +
  geom_point(data=Sel_fig_gathered, aes(x=value, y=r, color=Temp_range), size=1.5, shape=16) +
  geom_line(data=pred_df, aes(Variable, Value, color=Temp_range)) +
  scale_color_manual(values=c("orangered", "skyblue", "orange")) +
  facet_wrap(. ~ variables, scales = "free", ncol=2) +
  scale_y_continuous(limits=c(-6,13)) +
  labs(x="Variable", y=Intrinsic~growth~rate~(r)~(d^-1)) +
  theme(plot.background=element_blank(), panel.background=element_blank(),
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

### FIGURE 3E, 3F, 3G ----------



#G-matrix estimate
#creating reduced dataset for analyses
TPC_summary_spread_r<- TPC_summary_spread %>% select(Clone,CT_min, E_a, T_opt, r_peak)
TPC_summary_spread_r<-TPC_summary_spread_r[-c(2,10),] #filtering outliers

# Obtaining "test" temperatures
# itemps<-round(seq(1,381, length.out=6)[-1]/10)*10
# temps<-unique(TPC_predicted$Temperature)[itemps]
temps<-c(13, 19, 22, 25.1, 30, 32, 38)

# Building Gw (sensu Stinchcombe et al 2014; https://doi.org/10.1111/evo.12321)
# it requires obtaining the temperature specific rs (fitness) for each lineage
# and calculating the covariance among all traits, fitness included
Gwdist<-
  foreach(i=temps) %do% {
    joined<-right_join(TPC_summary_spread_r,
                       subset(TPC_predicted, round(Temperature,1)==i & 
                                Clone %in% TPC_summary_spread$Clone) %>% 
                         select(Clone,r))
    joined<-na.omit(joined)
    joined$W<-exp(joined$r)
    joined$W<-joined$W/mean(joined$W)
    # joined[,-1]<-sweep(joined[,-1],MARGIN = 2,STATS = colMeans(joined[,-1]),FUN = "/")
    out<-lm(cbind(CT_min,E_a,T_opt,r_peak,W)~1, data=joined) %>%
      evolqg::BayesianCalculateMatrix(.,samples=1000)
    out$Ps*0.5
  }
names(Gwdist)<-c(13, 19, 22, 25, 30, 32, 38)

traits<-c("r_peak","E_a", "CT_min", "T_opt")
colTraits<-RColorBrewer::brewer.pal(4, "Accent")

#Figure 3G, plotting the evolutionary responses according to Price equation
deltaz.plot<-
  llply(Gwdist, function(x)  {
    foreach(i=1:1000,.combine = "rbind") %do%{
      dz<-x[i,traits,"W"]
      G<-x[i,traits,traits]
      (dz)*sqrt(diag(G))
    }
  }) %>%
  melt %>% 
  subset(., value>-5 &value<5) %>%
  mutate(., Var2=factor(Var2, traits)) %>%
  ggplot(., aes(value, L1))+
  geom_vline(aes(xintercept=0), linetype=2)+
  ggridges::geom_density_ridges(aes(fill=Var2), alpha=0.4, show.legend = F,
                      quantile_lines = TRUE,quantiles=c(0.025,0.975))+
  facet_grid(.~Var2, scales = "free", labeller = label_parsed)+
  scale_fill_manual(name="Variable",values=colTraits)+
  ylim(names(Gwdist))+
  ylab(expression(paste(Temp,"(",C^o,")",)))+
  xlab(expression(paste("Predicted evolutionary change (", Delta,"z)")))+
  theme(axis.ticks.length=unit(-0.15, "cm"),
        strip.text.x = element_blank(),
        strip.background = element_blank())

#Figure 3F, plotting the estimated gradient of selection according to the breeder's equation
beta.plot<-
  llply(Gwdist, function(x)  {
    foreach(i=1:1000,.combine = "rbind") %do%{
      dz<-x[i,traits,"W"]
      G<-x[i,traits,traits]
      (dz%*%solve(G))*sqrt(diag(G))
    }
  }) %>%
  melt %>%
  mutate(Var2=factor(Var2,levels=traits)) %>%
  mutate(Var2=recode_factor(Var2,
                            CT_min="CT[min]",
                            T_opt ="T[opt]",
                            E_a   ="E[a]",
                            r_peak="r[peak]")) %>%
  mutate(., Var2=factor(Var2, c("r[peak]","E[a]","CT[min]","T[opt]"))) %>%
  ggplot(., aes(value, L1))+
  geom_vline(aes(xintercept=0), linetype=2)+
  ggridges::geom_density_ridges(aes(fill=Var2), alpha=0.4, show.legend = F,
                      quantile_lines = TRUE,quantiles=c(0.025,0.975))+
  facet_grid(.~Var2, scales = "free", labeller = label_parsed)+
  # xlim(c(-10,10))+
  scale_fill_manual(name="Variable",values=colTraits)+
  ylim(names(Gwdist))+
  ylab(expression(paste(Temp,"(",C^o,")",)))+
  xlab(expression(paste("Multivariate selection gradient (",beta,")")))+
  theme(strip.background = element_blank(),
        axis.ticks.length=unit(-0.15, "cm"))

#Figure 3E, trait plots of genetic correlation
ea_ctmin.plot<-
  ggplot(TPC_summary_spread_r, aes(E_a,CT_min))+ 
  stat_ellipse(geom = "polygon",fill = "black", alpha = 0.25)+
  geom_point(aes(color=Clone), show.legend = FALSE)+
  xlab(expression(E[a]))+
  ylab(expression(CT[min]))

topt_ctmin.plot<-
  ggplot(TPC_summary_spread_r, aes(T_opt,CT_min))+
  stat_ellipse(geom = "polygon",fill = "black", alpha = 0.25)+
  geom_point(aes(color=Clone), show.legend = FALSE)+
  xlab(expression(T[opt]))+
  ylab(expression(CT[min]))

ea_rpeak.plot<-
  ggplot(TPC_summary_spread_r, aes(E_a,r_peak))+   
  stat_ellipse(geom = "polygon",fill = "black", alpha = 0.25)+
  geom_point(aes(color=Clone), show.legend = FALSE)+
  xlab(expression(E[a]))+
  ylab(expression(r[peak]))

### FIGURE 4D, 4E, 4F -------------------------------------

# Process data from mar10data and mar10_NAdata based on replicant number/control condition
mar10dt <- rbind(mar10data, mar10_NAdata) 
mar10dt$temp = as.numeric(mar10dt$temp)
mar10dt$gate[which(mar10dt$gate=="AutoF")] <- "Auto"
mar10dt$gate[which(mar10dt$gate=="YFP.")] <- "YFP"
mar10dt$controlrep[which(mar10dt$controlrep=="2.cs")] <- "2"
mar10dt$controlrep[which(mar10dt$controlrep=="3.cs")] <- "3"
mar10dt$controlrep[which(mar10dt$controlrep=="4.cs")] <- "4"
mar10dt$controlrep[which(mar10dt$controlrep=="5.cs")] <- "5"
mar10dt$controlrep[which(mar10dt$controlrep=="6.cs")] <- "6"
mar10dt$controlrep[which(mar10dt$controlrep=="7.cs")] <- "7"
mar10dt$controlrep[which(mar10dt$controlrep=="8.cs")] <- "8"
mar10dt$controlrep[is.na(mar10dt$controlrep)] <- "1"

### Data analysis

# Group data by temperature, replicant number, genotype "gate", antibiotic condition, and control replicant number
# "gate" indicates what strain each individual is captured/designated as by flow cytometry imaging.
# Individuals were designated as either falling in gate "AutoF" or gate "YFP"
dt_sum0 <- mar10dt %>% 
  dplyr::group_by(temp, ab, rep, gate, controlrep) %>%
  # Length of the variable FSC_H indicates the number of individuals captured per microcosm in each gate
  # For example, a count of 16 means that 16 individuals were detected as auto fluorescing ("AutoF") CU4106 individuals when grown from a CU1406 microcosm at 19C in antibiotic+ conditions
  dplyr::summarize(count = length(FSC_H), 
            FSC_H_mean = mean(FSC_H), 
            SSC_H_mean = mean(SSC_H), 
            FITC_H_mean = mean(FITC_H), 
            PE_H_mean = mean(PE_H),
            PerCP_H_mean = mean(PerCP_H),
            APC_H_mean = mean(APC_H),
            FSC_A_mean = mean(FSC_A),
            SSC_A_mean = mean(SSC_A),
            FITC_A_mean = mean(FITC_A),
            PE_A_mean = mean(PE_A),
            PerCP_A_mean = mean(PerCP_A),
            APC_A_mean = mean(APC_A))
  # NA values mean that zero individuals were detected for the given gate, mark as such
dt_sum0$count[is.na(dt_sum0$FSC_H_mean)] <- 0
dt_sum0$rep <- as.factor(dt_sum0$rep)
dt_sum0$ab <- as.factor(dt_sum0$ab)
# Remove irrelevant fluorescent microscopy variables collected
dt_sum <- dt_sum0[-c(7:18)] 

### Controls: proportions 
# Calculate proportions of each genotype (AXS, denoted as "ControlYFP" and CU4106, denoted as "ControlCU")
# We use our control strain proportions to accurately calculate the true proportion of CU4106 and AXS in our competition microcosms
# We will use these proportions to adjust our experimental counts to ensure they accurately reflect control conditions
ctrlcounts <- dt_sum %>% 
  filter(rep == "ControlCU" | rep == "ControlYFP")
ctrlcounts <- ctrlcounts %>% 
  dplyr::group_by(temp, rep, ab, controlrep) %>% 
  # Sum is the total number of individuals detected for each level of treatment in both gates combined
  dplyr::mutate(sum = sum(count)) %>% 
  # X is the proportion of individuals that are one strain; for example, a proportion of 0.95 indicates that 95% of individuals grown from a single-strain CU4106 control microcosm were detected as auto-fluorescing individuals
  mutate(X = count/sum)
# Ensure all variables are numeric
ctrlcounts$sum <- as.numeric(ctrlcounts$sum)
ctrlcounts$X <- as.numeric(ctrlcounts$X)

# Remove all rows of replicated proportions (for example, row 5 is irrelevant because it provides the same information as row 21)
# Leave only the rows that include AXS individuals detected under YFP gating
ctrlcounts$X[ctrlcounts$gate=="Auto" & ctrlcounts$rep=="ControlYFP"] <- NA
ctrlcounts$X[ctrlcounts$gate=="Auto" & ctrlcounts$rep=="ControlCU"] <- NA
ctrlcounts$X[ctrlcounts$gate=="YFP" & ctrlcounts$rep=="ControlCU"] <- NA
prp <- ctrlcounts %>% na.omit() %>%
  # Find the mean proportion x of each treatment condition, which indicates on average how many YFP (i.e., AXS) individuals were detected
  dplyr::group_by(temp, rep, ab) %>% 
  dplyr::mutate(x = mean(X))

# Now that we have the average proportion, we don't need all control replicant data
prp1 <- prp %>% 
  filter(controlrep==1)

# Filter for only YFP gate and combine with average proportion dataset "prp1" so each microcosm has its corresponding control proportion assigned
dt_sum2 <- dt_sum %>% 
  filter(rep != "ControlCU", rep != "ControlYFP", gate != "Auto") %>% 
  arrange(temp, ab, rep, gate, controlrep) %>% ungroup %>% mutate(proportion = rep(prp1$x, each = 7))

# Generate final adjusted experimental counts of YFP/AXS based on control proportions
dt_adjCountsYFP <- dt_sum2 %>% 
  dplyr::mutate(adjCount=dt_sum2$count/dt_sum2$proportion)
# Create dataset of final experimental counts of CU4106
dt_adjCountsCU <- dt_sum %>% 
  filter(gate == "Auto", rep != "ControlCU", rep != "ControlYFP")

# Create merged data set where we can see the proportions and counts of each genotype
total <- merge(dt_adjCountsCU,dt_adjCountsYFP,by=c("temp","rep", "ab")) %>% 
  mutate(total=count.x+count.y, CU= total-adjCount) %>% 
  dplyr::rename("adjCountYFP" = "adjCount", "adjCountCU" = "CU") %>% 
  mutate(proportionCU=adjCountCU/total) %>% 
  mutate(proportionYFP=adjCountYFP/total) %>% 
  # We filter for counts greater than 0
  filter(adjCountCU>0)

# Clean up data set variable names
total$gate.y[which(total$gate.y=="YFP.")] <- "YFP"
totalYFP <- total[ , c("temp",  "rep", "ab", "gate.y", "count.y","adjCountYFP", "total",  "proportionYFP")]
totalYFP <- dplyr::rename(totalYFP, 
                   count = count.y, 
                   gate = gate.y, 
                   adjCount = adjCountYFP, 
                   proportion = proportionYFP)
totalCU <- total[ , c("temp",  "rep", "ab", "gate.x", "count.x", "total", "adjCountCU", "proportionCU")]
totalCU <- dplyr::rename(totalCU, 
                  count = count.x, 
                  gate = gate.x, 
                  adjCount = adjCountCU, 
                  proportion = proportionCU)
total_final <- rbind(totalCU, totalYFP)
total_final$temp<-as.factor(total_final$temp)

total_final2 <- total_final %>% 
  dplyr::group_by(ab, gate, temp) %>% 
  dplyr::summarise(avg = mean(proportion)) %>% 
  mutate(antibiotics = case_when(ab == "AB" ~ "Antibiotics", ab == "NoAB" ~ "No Antibiotics"))

## Appendix S7 Table
AppendixS7 <- total %>% 
  select(temp, rep, ab, proportion, count.x, adjCountCU, count.y, adjCountYFP, total)

# Split data set by antibiotic condition
total_final_ab <- total_final2 %>% 
  filter(antibiotics == "Antibiotics")
total_final_noab <- total_final2 %>% 
  filter(antibiotics == "No Antibiotics")

# Figure 4e displaying genotype proportion in antibiotic + conditions
Fig4e <- ggplot(total_final_ab, aes(fill=gate, y=avg, x=temp)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#de95c0","#C3D48A")) +
  scale_y_continuous(limits=c(0, 1), breaks =c(0, 0.5,1)) +
  labs(x="Temperature (C)", y="Frequency") +
  theme(plot.title = element_text(size=18, face = "bold"),
        plot.background=element_blank(), panel.background=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.text=element_text(size =16), axis.title=element_text(size=18),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        legend.position = c(30, 30),
        strip.background =element_blank(),
        strip.text.x = element_text(size = 16))

# Format Figure 4e for use in Appendix S8
Fig4e1 <- ggplot(total_final_ab, aes(fill=gate, y=avg, x=temp)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#de95c0","#C3D48A")) +
  scale_y_continuous(limits=c(0, 1), breaks =c(0, 0.5,1)) +
  labs(title = "Antibiotics", x="Temperature (C)", y="Frequency") +
  theme(plot.title = element_text(size=18, face = "bold"),
        plot.background=element_blank(), panel.background=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.text=element_text(size=16), axis.title=element_text(size=18),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        legend.position = c(30, 30),
        strip.background =element_blank(),
        strip.text.x = element_text(size = 16))

# Figure 4h displaying genotype proportion in antibiotic - conditions, for use in Appendix S8 Panel D
Fig4h <- ggplot(total_final_noab, aes(fill=gate, y=avg, x=temp)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#de95c0","#C3D48A")) +
  scale_y_continuous(limits=c(0, 1), breaks =c(0, 0.5,1)) +
  labs(title = "No Antibiotics", x="Temperature (C)", y="Frequency") +
  theme(plot.title = element_text(size=18, face = "bold"),
        plot.background=element_blank(), panel.background=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.text=element_text(size=16), axis.title=element_text(size=18),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        legend.position = c(30, 30),
        strip.background =element_blank(),
        strip.text.x = element_text(size = 16))

# Display antibiotic treatment conditions side by side
AppendixFigure <- ggarrange(Fig4e1, Fig4h)

##Figure 4D, displayed model predictions of frequencies across temperatures

# Fit model to the TPCs for clones AXS and CU4106
TPC_fits_exp <- tpc_exp_pred %>%
  dplyr::group_by(Clone) %>%
  do(TPC_fit = nls_multstart(log_r ~ a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temp+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temp+273.15)))),
                             data = .,
                             iter = 500,
                             start_lower = c(a=-10, E_a=0.1, E_d=0.5, Th=285),
                             start_upper = c(a=10, E_a=4, E_d=10, Th=330),
                             supp_errors = 'Y',
                             na.action = na.omit,
                             lower = c(a=-10, E_a=0, E_d=0, Th=0))) %>%
  rowwise() %>%
  dplyr::mutate(a=coef(TPC_fit)[[1]], E_a=coef(TPC_fit)[[2]], E_d=coef(TPC_fit)[[3]], Th=coef(TPC_fit)[[4]], T_opt=E_d*Th/(E_d+8.6e-5*Th*log(E_d/E_a-1)))

pars <- dplyr::select(TPC_fits_exp, -TPC_fit)  
TPC_eqn<-function(a, E_a, E_d, Th, Temperature, r_scale){exp(a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temperature+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temperature+273.15)))))-r_scale}

TPC_predicted_exp<-expand.grid(Clone=c("AXS", "CU4106"), Temperature=seq(0, 50, length.out=500)) %>%
  left_join(dplyr::select(TPC_fits_exp, -TPC_fit)) %>%
  left_join(distinct(dplyr::select(tpc_exp_pred, Clone, r_scale))) %>%
  mutate(r=TPC_eqn(a, E_a, E_d, Th, Temperature, r_scale))

# Plot all species' TPCs together
TPC_predicted_exp$Clone<-factor(TPC_predicted_exp$Clone)

# Use the TPCs to predict genotypic frequencies across temperatures

# Calculate mean_fitness and prepare data 
mean_fitness <- TPC_predicted_exp %>%
  dplyr::group_by(Temperature) %>%               
  dplyr::summarize(meanFit=mean(r), SD=sd(r), dissim=mean(dist(r))) 

TPC_predicted_exp <- merge(TPC_predicted_exp,mean_fitness,by="Temperature")
TPC_predicted_exp$rel_fitness <- (TPC_predicted_exp$r+10)/(TPC_predicted_exp$meanFit+10)
Temp <- unique(TPC_predicted_exp$Temperature)

## Define function for model
## Temperature has to be passed as Temp[nmbr] because it has to match existing temperatures in the TPC_predict dataset that contains the TPCs that will be used to make the average fitness predictions.
evo_mod <- function(temp,noise,time_steps){
  ## Select appropriate temperature data
  fit <- TPC_predicted_exp %>%
    filter(Temperature==temp) %>%
    select(r)
  ## Run model
  N_clones <- 2
  time <- seq(1,time_steps,1)
  mat_freq <- matrix(rep(0,N_clones*time_steps),nrow=N_clones,ncol=time_steps)
  mat_freq[,1] <- rep(1/N_clones,N_clones)	
  for(i in 1:(time_steps-1)){
    if(i==1){ #We add +4 to r to avoid dealing with negative vaues of r and how it impacts the calculation of rel.fit.
      new_freq <- pmax(diag((fit$r+4)/(mean(fit$r)+4),N_clones,N_clones)%*%rep(1/N_clones,N_clones) + runif(N_clones,min=-noise,max=noise),0) # pmax turns	negative frequencies into 0s
      mat_freq[,2] <- new_freq/sum(new_freq) # Dividing by the sum ensures that sum(frequencies)=1. Rescales all frequencies so that losing alleles doesn't lower the sum(frequencies). Needed due to randomness 
    }else{
      mat_freq[,i+1] <-pmax(diag((fit$r+4)/sum(mat_freq[,i-1]*(fit$r+4)),N_clones,N_clones)%*%mat_freq[,i]+ runif(2,min=-noise,max=noise),0) # Here 2 is the number of clones
      mat_freq[,i+1] <- mat_freq[,i+1]/sum(mat_freq[,i+1])
    }
  }
  return(rowMeans(mat_freq[,1:time_steps]))
}

# Check that function is working
mat <- evo_mod(Temp,0.05,100)
runs <- sapply(Temp,evo_mod,0.00,20)

# Prepare data for plotting
runs_clone <- as.data.frame(runs) %>% mutate(Clone=c("AXS", "CU4106")) %>% relocate("Clone")
colnames(runs_clone) <- c("Clone", Temp)

new_runs <- runs_clone %>%
  gather("Temp","Freq",2:501) %>%
  mutate(Temp=as.numeric(Temp)) 

new_runs$Clone <- as.factor(new_runs$Clone)
new_runs$Clone <- relevel(new_runs$Clone, "CU4106")

# Plot Figure 4f displaying model predictions in changes in frequency across temperature
Fig4f <- ggplot(new_runs,aes(x=Temp, y=Freq, fill=Clone)) + 
  geom_area(linewidth=0.1, colour="white") +
  scale_x_continuous(limits=c(10, 38)) +
  scale_fill_manual(name = "Clone", values = c("#de95c0", "#C3D48A"))+
  labs(x="Temperature (C)", y="Frequency") +
  theme(plot.background=element_blank(), panel.background=element_blank(),#panel.grid=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=16), axis.title=element_text(size=18),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        legend.position="none")

TPC_summary_spread<-TPC_predicted_exp %>%
  dplyr::group_by(Clone) %>%
  dplyr::mutate(CT_min=ifelse(lag(r)<0 & r>0, Temperature, NA),
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

#Remove all irrelevant strains from data set

TPC_forselection_1 <- TPC_predicted_exp %>%
  filter(Clone=="CU4106") %>%
  select(Temperature,r) %>%
  mutate(r_star=r, r=NULL)

TPC_forselection <- merge(TPC_predicted_exp,TPC_forselection_1,by="Temperature") %>%
  group_by(Clone) %>%
  mutate(Diff=(r+4)/(r_star+4), sel=1-Diff)	

# Remove all irrelevant TPC parameters from data set
Fig4ddata <- TPC_forselection[, -c(3:14, 16)]

# Plot Figure 4d inset displaying relative fitness of AXS to CU4106
Fig4dinset <- ggplot(data=Fig4ddata, aes(x=Temperature, y=Diff, color=Clone))+
  geom_line(linewidth =1) +
  scale_x_continuous(limits=c(10, 37)) +
  scale_y_continuous(limits=c(0.75,1.1)) +
  scale_color_manual(name = "Clone", values = c("AXS" ="#C3D48A", "CU4106" ="#de95c0")) +
  labs(x="Temperature (C)", y="Relative Fitness") +
  theme(plot.background=element_blank(), panel.background=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=16), axis.title=element_text(size=18),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        legend.position=c(0.80, 0.28),
        legend.background = element_rect(color = NULL),
        legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size=14),
        legend.title= element_text(size=15))

## Figure 4d

Clones <- unique(tpc_exp_pred$Clone)
tpc_exp_pred <- tpc_exp_pred %>%
  mutate(unique_rep=paste(tpc_exp_pred$Clone,".",tpc_exp_pred$Rep,sep="")) %>% 
  mutate(Clone = replace(Clone, Clone =="YFP_young", "AXS")) %>% 
  mutate(Clone = replace(Clone, Clone == "CU4106_young", "CU4106"))
# Calculate r across temperatures and generate data set
tpc_exp_pred <- tpc_exp_pred %>%
  filter(Final>=1) %>%
  mutate(r = log(Final/Initial)) %>%
  group_by(Clone) %>%
  mutate(r_scale=10,
         log_r=log(r+r_scale)) %>%
  ungroup
tpc_exp_pred = subset(tpc_exp_pred, select = -c(unique_rep) )
tpc_exp_pred$Clone <- as.factor(tpc_exp_pred$Clone)

# Use nls.multstart package to fit the curve to TPC data
TPC_fits_exp <- tpc_exp_pred %>%
  dplyr::group_by(Clone) %>%
  do(TPC_fit = nls_multstart(log_r ~ a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temp+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temp+273.15)))),
                             data = .,
                             iter = 500,
                             start_lower = c(a=-10, E_a=0.1, E_d=0.5, Th=285),
                             start_upper = c(a=10, E_a=4, E_d=10, Th=330),
                             supp_errors = 'Y',
                             na.action = na.omit,
                             lower = c(a=-10, E_a=0, E_d=0, Th=0))) %>%
  rowwise() %>%
  dplyr::mutate(a=coef(TPC_fit)[[1]], 
         E_a=coef(TPC_fit)[[2]], 
         E_d=coef(TPC_fit)[[3]], 
         Th=coef(TPC_fit)[[4]], 
         T_opt=E_d*Th/(E_d+8.6e-5*Th*log(E_d/E_a-1)))

pars <- dplyr::select(TPC_fits_exp, -TPC_fit)

TPC_eqn<-function(a, E_a, E_d, Th, Temperature, r_scale){exp(a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temperature+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temperature+273.15)))))-r_scale}

TPC_predicted_exp<-expand.grid(Clone=TPC_fits_exp$Clone, Temperature=seq(0, 50, length.out=500)) %>%
  left_join(dplyr::select(TPC_fits_exp, -TPC_fit)) %>%
  left_join(distinct(dplyr::select(tpc_exp_pred, Clone, r_scale))) %>%
  mutate(r=TPC_eqn(a, E_a, E_d, Th, Temperature, r_scale))

TPC_summary_spread_exp<-TPC_predicted_exp %>%
  dplyr::group_by(Clone) %>%
  dplyr::mutate(CT_min=ifelse(lag(r)<0 & r>0, Temperature, NA),
         CT_max=ifelse(lag(r)>0 & r<0, Temperature, NA),
         r_peak=max(r)) %>%
  filter(!is.na(CT_min) | !is.na(CT_max)) %>%
  dplyr::select(Clone, CT_min, CT_max, r_peak) %>%
  ungroup() %>%
  gather(param, param_val, -Clone) %>%
  drop_na %>%
  distinct %>%
  spread(param, param_val) %>%
  left_join(dplyr::select(TPC_fits_exp, Clone, E_a, E_d, T_opt)) %>%
  mutate(T_opt=T_opt-273.15, T_range=CT_max-CT_min, TPC_asymmetry=abs((T_opt-CT_min)-(CT_max-T_opt))) %>%
  arrange(r_peak)        

# Plot all species' TPCs together
TPC_predicted_exp$Clone<-factor(TPC_predicted_exp$Clone)

# Calculate mean fitness curve 
mean_fitness <- TPC_predicted_exp %>%
  group_by(Temperature) %>%               
  summarize(meanFit=mean(r), SD=sd(r), dissim=mean(dist(r)))     

# Figure 4d, demonstrating the TPC for genotypes AXS and CU4106. Combine with above Figure 4d inset for final presented version.
Fig4d<-
  ggplot()+
  geom_hline(yintercept=0, color="gray", linewidth=0.7, linetype=2)+
  geom_line(data=TPC_predicted_exp, aes(Temperature, r, color=Clone), linewidth=2.5)+
  scale_color_manual(name = "Clone", values = c("AXS" ="#C3D48A", "CU4106" ="#de95c0")) +
  scale_x_continuous(limits=c(10, 38)) +
  scale_y_continuous(limits=c(-2,20)) +
  labs(x="Temperature (C)", y=Intrinsic~growth~rate~(r)~(cells~cell^-1~d^-1)) +
  theme(plot.background=element_blank(), panel.background=element_blank(),
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

# Generate Figure 4 with blank inset for microscopy images in 4c, 4d inset, 4e, and 4f
blank <- grid.rect(gp=gpar(col="white"))
pdf("/Users/meganliu/Desktop/Gibert_lab_project/Figures/Figure2", width = 10,  height = 4)
Fig4 <- ggarrange(blank, Fig4d, Fig4e, Fig4f,
                  nrow = 2,
                  ncol = 2,
                  labels = c("c)", "d)", "e)", "f)"),
                  heights = c(1, 1, 1, 1, 1),
                  widths = c(1, 1, 1, 1, 1))
print(Fig4)
dev.off()

# Finalize Appendix Figure S8
Fig4appendix <- ggarrange(Fig4d, Fig4f, Fig4e1, Fig4h, 
                          ncol = 2, 
                          nrow = 2,
                          labels = c("a)", "b)", "c)", "d)"),
                          heights = c(1, 1, 1, 1),
                          widths = c(1, 1, 1, 1))

pdf("/Users/meganliu/Desktop/Gibert_lab_project/Figures/Figure4appendix", width = 10,  height = 7)
print(Fig4appendix)
dev.off() 
