#################################################################
########## Code for Liu et al Communications Biology ############
################################################################# 

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
install.packages("ggnewscale")
install.packages("MCMCglmm")
install.packages("rsample")
install.packages("purrr")
install.packages("furrr")
install.packages("rsample")


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
library("purrr")
library("furrr")
library("MCMCglmm")
library("patchwork")
library("rsample")


 
### LOAD AND PROCESS DATA -----------

## Load TPC data from initial 22 genotypes
tpcs <- read.csv("~/Desktop/JP/Papers_in_review_submitted/2_Resubmitted/Megan_TPCs/MS/Communications_Biology/Resubmission/New_Code/Analyses/Data/TPC_data.csv")
Clones <- unique(tpcs$Clone)
tpcs %>%
  mutate(unique_rep=paste(tpcs$Clone,".",tpcs$Rep,sep=""))

## Load TPC data for genotypes AXS and CU4106, used in rapid evolution experiment 
tpc_exp_pred <- read.csv("~/Desktop/JP/Papers_in_review_submitted/2_Resubmitted/Megan_TPCs/MS/Communications_Biology/Resubmission/New_Code/Analyses/Data/two_clone_TPC.csv") %>% mutate( r_scale=10,log_r=log(r+r_scale)) 
tpc_exp_pred$Clone[tpc_exp_pred$Clone =="CU416"] <- "CU4106"

## Load data from competition experiment between genotypes AXS and CU4106. 
# Non-compiled data available in file "Flow_Cyt_Data.zip"
# YFP indicates individuals that were "gated" or captured as strain AXS by flow cytometry
# AutoF indicates individuals that were "gated" or captured as strain CU4106 by flow cytometry
mar10data <- read.csv("~/Desktop/JP/Papers_in_review_submitted/2_Resubmitted/Megan_TPCs/MS/Communications_Biology/Resubmission/New_Code/Analyses/Data/FlowCytData.csv")
mar10data <- mar10data[-c(1)]
colnames(mar10data) <- c("id", "FSC_H", "SSC_H","FITC_H", "PE_H", "PerCP_H", "APC_H","FSC_A","SSC_A","FITC_A", "PE_A", "PerCP_A",
                         "APC_A","Width","Time", "temp", "ab", "rep","gate", "controlrep")
mar10_NAdata <- read.csv("~/Desktop/JP/Papers_in_review_submitted/2_Resubmitted/Megan_TPCs/MS/Communications_Biology/Resubmission/New_Code/Analyses/Data/NA_Values.csv")
colnames(mar10_NAdata) <-  c("id", "FSC_H", "SSC_H","FITC_H", "PE_H", "PerCP_H", "APC_H","FSC_A","SSC_A","FITC_A", "PE_A", "PerCP_A",
                             "APC_A","Width","Time", "temp", "ab", "rep","gate", "controlrep" )

###############################################################################################
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

# Calculations show trait is very heritable; Standard H^2=0.746, H^2 (Cullis)=0.911, H^2 (Piepho) = 0.953. See inti package vignette for details.

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


## CODE TO TEST EFFECT OF SCALE AND CODE USED FOR RESULTS IN MAIN TEXT
tpcs <- tpcs %>%
  filter(Final>=1) %>%
  mutate(r = log(Final/Initial)) %>%
  dplyr::group_by(Clone) %>%
  mutate(r_scale=2,
         log_r=log(r+r_scale)) %>%
  ungroup

## Code to fit TPCs in log space to which 
TPC_fits <- tpcs %>%
  dplyr::group_by(Clone) %>%
  do(
  TPC_fit = nls_multstart(
  log_r ~ a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temp+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temp+273.15)))),
  #TPC_fit = nls_multstart(log_r ~ a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temp+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temp+273.15)))),
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
         Th=coef(TPC_fit)[[4]]#, 
         #T_opt=E_d*Th/(E_d+8.6e-5*Th*log(E_d/E_a-1))
         )

####-----------------------------------------------------------------------------
pars <- TPC_fits
pars <- dplyr::select(TPC_fits, -TPC_fit) # Use this for fits with scale

# set TPC equation
TPC_eqn<-function(a, E_a, E_d, Th, Temperature, r_scale){exp(a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temperature+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temperature+273.15)))))-r_scale}

TPC_predicted<-expand.grid(Clone=Clones, Temperature=seq(0, 50, length.out=500)) %>%
  left_join(pars) %>%
  left_join(distinct(dplyr::select(tpcs, Clone, r_scale))) %>%
  mutate(r=TPC_eqn(a, E_a, E_d, Th, Temperature, r_scale))

TPC_summary_spread<-TPC_predicted %>%
  dplyr::group_by(Clone) %>%
  dplyr::mutate(CT_min=ifelse(lag(r)<0 & r>0, Temperature, NA),
         CT_max=ifelse(lag(r)>0 & r<0, Temperature, NA),
         r_peak=max(r),
    		 T_opt = ifelse(Temperature[which.max(r)]<45, Temperature[which.max(r)], NA)
         ) %>%
  filter(!is.na(CT_min) | !is.na(CT_max)) %>%
  dplyr::select(Clone, CT_min, CT_max, r_peak, T_opt) %>%
  ungroup() %>%
  gather(param, param_val, -Clone) %>%
  drop_na %>%
  distinct %>%
  spread(param, param_val) %>%
  left_join(dplyr::select(TPC_fits, Clone, E_a, E_d, Th, a)) %>%
  left_join(dplyr::select(TPC_fits, Clone, E_a, E_d, Th)) %>%
  #dplyr::mutate(T_opt=T_opt, T_range=CT_max-CT_min, TPC_asymmetry=abs((T_opt-CT_min)-(CT_max-T_opt))) %>%
  dplyr::mutate(T_opt=T_opt) %>%
  arrange(r_peak)        

TPC_summary_spread_2 <- TPC_summary_spread
print(TPC_summary_spread, n=22)

kB <- 8.6e-5            # Boltzmann constant, eV K⁻¹
ea_numeric <- tpcs %>%
  dplyr::inner_join(dplyr::select(TPC_summary_spread, Clone, CT_min, T_opt),
             by = "Clone") %>%
  dplyr::filter(Temp >= CT_min,
         Temp <= T_opt,
         r           >  0) %>%           # log() is safe now
  dplyr::group_by(Clone) %>%
  dplyr::summarise(
    Ea_emp = if (dplyr::n() < 3) NA_real_ else {
      slope <- coef(
        lm(log(r) ~ I(1 / (Temp + 273.15)))
      )[2]
      -kB * slope
    },
    .groups = "drop"
  )

# add the empirical Ea (Ea_emp) to summary table
TPC_summary_spread <- TPC_summary_spread %>%
  left_join(ea_numeric, by = "Clone")

pl <- ggplot(TPC_summary_spread,
       aes(x = E_a, y = Ea_emp, colour = Clone, label = Clone)) +

  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60") +            # 1 : 1 reference

  geom_point(size = 3) +                      # dots for each clone
  # ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +  # tidy labels
  geom_text(vjust = -0.6, size = 3, show.legend = FALSE) +    # simple labels

  scale_colour_viridis_d(option = "turbo") +  # palette for Clone
  coord_equal() +                             # square plotting area

  labs(x = expression(Model~E[a]~"(eV)"),
       y = expression(Empirical~E[a]~"(eV)"),
       colour = "Clone") +

  theme_bw() +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "none")      


stats_df <- data.frame(
  parameter = names(TPC_summary_spread)[-1],          # drop Clone
  mean      = sapply(TPC_summary_spread[-1], mean, na.rm = TRUE),
  sd        = sapply(TPC_summary_spread[-1],  sd,   na.rm = TRUE)
)

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
  labs(x="Temperature (ºC)", y=Intrinsic~growth~rate~(r)~(cells~cell^-1~d^-1)) +
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
#############################################################################################


##############################################################################################
######## APPENDIX ONE-TO-ONE PLOTS OF SCALE COMPARISON TO TEST DIFFERENCES IN SCALE (APPENDIX FIGURES)\
		## This is not needed to make main text figures.

# Variables to compare 
vars <- c("CT_min", "r_peak", "E_a", "T_opt")

# Assemble a single data frame with side-by-side columns
combo <- TPC_summary_spread %>%                                   
  select(Clone, all_of(vars)) %>% 
  rename_with(~ paste0(.x, "_d1"), -Clone) %>%                    
  inner_join(                                                     
    TPC_summary_spread_2 %>%                                      
      select(Clone, all_of(vars)) %>% 
      rename_with(~ paste0(.x, "_d2"), -Clone),
    by = "Clone"
  )

## Plot of paramter estimates to compare two different model fits
square_limits <- function(x, y, pad = 0.05) {
  rng  <- range(c(x, y), na.rm = TRUE)          # min / max over both axes
  span <- diff(rng)
  rng + c(-1, 1) * span * pad                  # small padding on both ends
}

plots <- imap(vars, function(v, i) {
  x  <- combo[[paste0(v, "_d1")]]
  y  <- combo[[paste0(v, "_d2")]]
  lims <- square_limits(x, y)
  p <- ggplot(
        combo,
        aes_string(x = paste0(v, "_d1"),
                   y = paste0(v, "_d2"),
                   colour = "Clone",
                   label  = "Clone")
      ) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = "grey50") +
      geom_point(size = 3, na.rm = TRUE) +
      geom_text(vjust = -0.5, size = 3, show.legend = FALSE, na.rm = TRUE) +
      coord_equal(xlim = lims, ylim = lims, expand = FALSE) +  
      scale_colour_viridis_d(option = "turbo") +
      labs(title = v,
           x = paste(v, "(scale = 10)"),
           y = paste(v, "(scale = 2)")) +
      theme_bw() +
      theme(plot.title      = element_text(face = "bold"),
            legend.position = "none")

  # add linear model only for the 3rd variable
  if (i == 3) {
    p <- p + geom_smooth(method = "lm", se = FALSE,
                         colour = "black", linewidth = 0.8)
  }
  p
})

# 2x2 grid 
#(CAUTION, FROM NOW ON, THIS ONLY WORKS WHEN MODELS WITH DIFFERENT SCALES 
#OR PARAMETERS HAVE BEEN RAN)
final_plot <-
  wrap_plots(plots, nrow = 2) +                        
  plot_layout(guides  = "collect",                     
              widths  = c(1, 1),                       
              heights = c(1, 1))                       

### Test of rank-order preservation:
# Combine and rank the values per dataset 
ranked_combo <- TPC_summary_spread %>%
  select(Clone, all_of(vars)) %>%
  rename_with(~ paste0(.x, "_d1"), -Clone) %>%
  inner_join(
    TPC_summary_spread_2 %>%
      select(Clone, all_of(vars)) %>%
      rename_with(~ paste0(.x, "_d2"), -Clone),
    by = "Clone"
  ) %>%
  rowwise() %>%
  mutate(across(ends_with("_d1") | ends_with("_d2"), as.numeric)) %>%
  ungroup()

# Compute ranks within each variable for both datasets 
for (v in vars) {
  ranked_combo[[paste0(v, "_rank_d1")]] <- rank(ranked_combo[[paste0(v, "_d1")]], na.last = "keep")
  ranked_combo[[paste0(v, "_rank_d2")]] <- rank(ranked_combo[[paste0(v, "_d2")]], na.last = "keep")
}

# Plot rank–rank comparisons for each variable
rank_plots <- map(vars, function(v) {
  ggplot(ranked_combo,
         aes_string(x = paste0(v, "_rank_d1"),
                    y = paste0(v, "_rank_d2"),
                    label = "Clone",
                    colour = "Clone")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(size = 3, na.rm = TRUE) +
    geom_text(vjust = -0.5, size = 3, show.legend = FALSE, na.rm = TRUE) +
    coord_equal() +
    scale_colour_viridis_d(option = "turbo") +
    labs(title = paste("Rank comparison:", v),
         x = paste("Rank (scale = 10)"),
         y = paste("Rank (scale= 2)")) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
})

# Arrange plots in a grid 
(rank_plots[[1]] + rank_plots[[2]]) /
  (rank_plots[[3]] + rank_plots[[4]])


############################################################################################
##### FIGURE 3
#### Shape parameter G and correlation with r-----

# Generate full dataset and calculate log_r from available data
tpcs <- tpcs %>%
  filter(Final >= 1) %>%
  mutate(r = log(Final / Initial)) %>%
  dplyr::group_by(Clone) %>%
  mutate(r_scale = 2,
         log_r = log(r + r_scale)) %>%
  ungroup()

# Define the TPC equation
TPC_eqn <- function(a, E_a, E_d, Th, Temperature, r_scale) {
  exp(a + (E_a / (8.6 * 10^-5)) * (1 / 298.15 - 1 / (Temperature + 273.15)) -
        log(1 + exp((E_d / (8.6 * 10^-5)) * (1 / Th - 1 / (Temperature + 273.15))))) - r_scale
}

# Function to fit TPC and extract parameters
fit_tpc <- function(data) {
  tryCatch({
    fit <- nls_multstart(
      log_r ~ a + (E_a / (8.6 * 10^-5)) * (1 / 298.15 - 1 / (Temp + 273.15)) -
        log(1 + exp((E_d / (8.6 * 10^-5)) * (1 / Th - 1 / (Temp + 273.15)))),
      data = data,
      iter = 500,
      start_lower = c(a = -10, E_a = 0.1, E_d = 0.5, Th = 285),
      start_upper = c(a = 10, E_a = 4, E_d = 10, Th = 330),
      supp_errors = 'Y',
      na.action = na.omit,
      lower = c(a = -10, E_a = 0, E_d = 0, Th = 0)
    )
    tibble(
      a = coef(fit)[[1]],
      E_a = coef(fit)[[2]],
      E_d = coef(fit)[[3]],
      Th = coef(fit)[[4]]#,
      #T_opt = coef(fit)[[3]] * coef(fit)[[4]] / (coef(fit)[[3]] + 8.6e-5 * coef(fit)[[4]] * log(coef(fit)[[3]] / coef(fit)[[2]] - 1))
    )
  }, error = function(e) {
    tibble(a = NA, E_a = NA, E_d = NA, Th = NA)
  })
}


# Bootstrap and fit TPCs
plan(multisession, workers = parallel::detectCores() - 1)

# Perform bootstrapping and TPC fitting
split_data <- tpcs %>%
  group_by(Clone) %>%
  group_split()  # Creates a list of data frames, one per Clone

# Perform bootstrapping and TPC fitting
bootstrapped_tpcs <- future_map_dfr(split_data, ~ {
  clone_data <- .x
  clone_name <- unique(clone_data$Clone) 
  #map_dfr(1:100, ~ {
  	map_dfr(1:300, ~ {
    boot_data <- clone_data %>%
      slice_sample(n = nrow(clone_data), replace = TRUE)
    
    tpc_fit <- fit_tpc(boot_data)
    
    tpc_fit %>%
      mutate(Replicate = .x, Clone = clone_name)
  })
}, .options = furrr_options(seed = TRUE))  # for reproducibility


TPC_predictions <- bootstrapped_tpcs %>%
  group_by(Replicate, Clone) %>%
  mutate(Temperature = list(seq(13, 38, by = 0.1))) %>%  # Temperature range for prediction
  unnest(Temperature) %>%
  mutate(
    r = TPC_eqn(a, E_a, E_d, Th, Temperature, 2)  # Predicted r values
  )


### THIS PIECE OF CODE RECREATES APPENDIX FIGURE
## Check fits
# Plot the TPCs
ggplot(TPC_predictions, aes(x = Temperature, y = r, group = interaction(Clone, Replicate), color = Clone)) +
  geom_line(alpha = 0.3) +  # Use alpha for transparency to visualize overlapping lines
  facet_wrap(~ Clone) +     # Separate plots for each Clone
  theme_minimal() +
  labs(
    title = "Bootstrapped TPC Curves",
    x = "Temperature (°C)",
    y = "Intrinsic Growth Rate (r)",
    color = "Clone"
  )

## Now we create dataset for Fig 3 of main text:
TPC_summary_spread2 <- TPC_predictions %>%
  dplyr::group_by(Clone, Replicate) %>%
  arrange(Temperature) %>%  # Ensure data is sorted
  dplyr::summarise(
    CT_min = first(na.omit(ifelse(lag(r)<0 & r>0, Temperature, NA))),
    r_peak = max(r, na.rm = TRUE),
    E_a = median(E_a),   # Assuming E_a is constant for each Clone and Replicate
    T_opt = Temperature[which.max(r)],  # Assuming T_opt is constant for each Clone and Replicate
    r = median(r)
  ) %>% 
  ungroup()

TPC_summary_spread2 <- TPC_summary_spread2 %>%
  distinct(Clone, Replicate, CT_min, r_peak, E_a, T_opt, r) %>%
  dplyr::mutate(
    CT_min_scaled = scale(CT_min),  # Center and scale CT_min
    T_opt_scaled = scale(T_opt),    # Center and scale T_opt
    r_peak_scaled = scale(r_peak),  # Center and scale r_peak
    E_a_scaled = scale(E_a),         # Center and scale E_a
    r_scaled = scale(r)
  )%>%
  arrange(Clone, Replicate)

## Check for outliers
TPC_summary_spread2 %>%
  dplyr::group_by(Clone) %>%
  dplyr::summarise(
    CT_min_mean = mean(CT_min, na.rm = TRUE),
    CT_min_sd   = sd(CT_min, na.rm = TRUE),
    r_peak_mean = mean(r_peak, na.rm = TRUE),
    r_peak_sd   = sd(r_peak, na.rm = TRUE),
    E_a_mean    = mean(E_a, na.rm = TRUE),
    E_a_sd      = sd(E_a, na.rm = TRUE),
    T_opt_mean  = mean(T_opt, na.rm = TRUE),
    T_opt_sd    = sd(T_opt, na.rm = TRUE),
    .groups = "drop"
  ) %>% print(.,n = Inf, width = Inf)
## big outlier: CU427.4 (i=14) in Ea (4 times and ten times larger)

TPC_summary_spread2 <- TPC_summary_spread2 %>%
  						dplyr::group_by(Clone) %>% 
  						filter(Clone!="CU427.4")
# If needed, but removing these does not alter results

TPC_summary_long <- TPC_summary_spread2 %>%
  pivot_longer(cols = c(CT_min_scaled, r_peak_scaled, T_opt_scaled, E_a_scaled, r_scaled), names_to = "Variable", values_to = "Value")

# Create boxplots for each variable
prueba <- ggplot(TPC_summary_long, aes(x = Variable, y = Value)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Variation in TPC Parameters Across Clones",
    x = "TPC Parameters",
    y = "Values"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Calculate prior means
prior_means <- colMeans(TPC_summary_spread2[, c("CT_min_scaled", "T_opt_scaled", "r_peak_scaled", "E_a_scaled", "r_scaled")], na.rm = TRUE)

# Calculate prior covariance matrix
prior_cov <- cov(TPC_summary_spread2[, c("CT_min_scaled", "T_opt_scaled", "r_peak_scaled", "E_a_scaled", "r_scaled")], use = "complete.obs")
scaled_prior_cov <- prior_cov / mean(diag(prior_cov))

# Informative priors
prior <- list(
  R = list(V = scaled_prior_cov, nu = 5),  # Residual covariance
  G = list(G1 = list(V = scaled_prior_cov, nu = 5))  # Random effects covariance
)

#NON-Informative Priors
#prior <- list(
#  R = list(V = diag(5), nu = 5+1),  # Residual covariance: uninformative prior
#  G = list(
#    G1 = list(V = diag(5), nu = 5+1)  # Random effects covariance: uninformative prior
#  )
#)


#TPC_summary_spread2 <- TPC_summary_spread %>% 
#  select(-E_a) %>%          # remove the Sharpe–Schoolfield Ea
#  mutate(E_a = Ea_emp) %>% 
#  select(-Ea_emp)

model <- MCMCglmm(
  cbind(CT_min_scaled, T_opt_scaled, r_peak_scaled, E_a_scaled, r_scaled) ~ trait - 1,  # Multivariate response
  random = ~ us(trait):Clone,  # Random effects for Clone with unstructured covariance
  rcov = ~ us(trait):units,    # Residual covariance
  family = rep("gaussian", 5),  # Gaussian responses
  data = as.data.frame(TPC_summary_spread2),
  prior = prior,
  nitt = 13000, burnin = 3000, thin = 10  # Adjust as needed
)

summary(model)
plot(model)

# Extract genetic variance-covariance matrix (G)
G_matrix <- posterior.mode(model$VCV[, grepl("Clone", colnames(model$VCV))])
# Convert to a readable matrix format
G_matrix <- matrix(G_matrix, nrow = 5, ncol = 5)
rownames(G_matrix) <- colnames(G_matrix) <- c("CT_min", "T_opt", "r_peak", "E_a", "r")
G_matrix

# Extract residual variance-covariance matrix (R)
R_matrix <- posterior.mode(model$VCV[, grepl("units", colnames(model$VCV))])
# Convert to a readable matrix format
R_matrix <- matrix(R_matrix, nrow = 5, ncol = 5)
rownames(R_matrix) <- colnames(R_matrix) <- c("CT_min", "T_opt", "r_peak", "E_a", "r")
R_matrix

genetic_variances <- diag(G_matrix)
residual_variances <- diag(R_matrix)

# Calculate mean heritabilities
heritabilities <- genetic_variances / (genetic_variances + residual_variances)
names(heritabilities) <- c("CT_min", "T_opt", "r_peak", "E_a", "r")
heritabilities ## Suggests that CT_min and r_peak have more G than residual variance.

## Posterior distribution of heritabilities:
# Extract G and R columns from model$VCV
G_post <- model$VCV[, grepl("Clone", colnames(model$VCV))]
R_post <- model$VCV[, grepl("units", colnames(model$VCV))]

# Trait names in order
traits <- c("CT_min_scaled", "r_peak_scaled", "T_opt_scaled", "E_a_scaled")

# Match diagonals for each trait (same order in G and R)
herit_samples <- sapply(traits, function(trait) {
  G_diag <- G_post[, paste0("trait", trait, ":trait", trait, ".Clone")]
  R_diag <- R_post[, paste0("trait", trait, ":trait", trait, ".units")]
  G_diag / (G_diag + R_diag)
})

# Name columns for clarity
colnames(herit_samples) <- gsub("_scaled", "", traits)

# View summaries (posterior means and 95% credible intervals)
apply(herit_samples, 2, function(x) {
  c(mean = mean(x), quantile(x, probs = c(0.025, 0.975)))
})

# Convert matrix to long format
herit_df <- as.data.frame(herit_samples) %>%
  pivot_longer(cols = everything(), names_to = "Trait", values_to = "Heritability")

herit_fig <- ggplot(herit_df, aes(x = Heritability, fill = Trait, color = Trait)) +
  geom_density(alpha = 0.3, linewidth = 0.65) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)  # Add tick at 1 explicitly
  ) +
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "Heritability", y = "Density") +
  geom_hline(yintercept = 0, color = "gray", linewidth = 0.7) +
  theme(
    plot.background     = element_blank(),
    panel.background    = element_blank(),
    panel.border        = element_rect(color = "black", size = 1, fill = NA),
    axis.text           = element_text(size = 16),
    axis.title          = element_text(size = 16),
    axis.title.x.top    = element_blank(),
    axis.text.x.top     = element_blank(),
    axis.ticks.length.x.top = unit(0, "cm"),
    axis.ticks.length   = unit(-0.15, "cm"),
    plot.margin         = unit(c(1, 1, 1, 1), "cm"),
    strip.text          = element_text(size = 11),
    aspect.ratio        = 0.7,
    legend.key          = element_blank(),
    legend.position     = "none",
    legend.text         = element_text(size = 12, face = "italic")
  )

## Plot covariances:
# Traits to test for covariation with r
# Summarize posterior covariances with r
traits <- c("CT_min_scaled", "r_peak_scaled", "T_opt_scaled", "E_a_scaled")
cov_summary <- lapply(traits, function(trait) {
  vec <- G_post[, paste0("trait", trait, ":traitr_scaled.Clone")]
  tibble(
    Trait = gsub("_scaled", "", trait),
    mean = mean(vec),
    lower = quantile(vec, 0.025),
    upper = quantile(vec, 0.975)
  )
}) %>%
  bind_rows()

cov_simple <- ggplot(cov_summary, aes(y = Trait, x = mean)) +
  geom_vline(xintercept = 0, color = "gray40", linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_y_discrete(limits = rev(unique(cov_summary$Trait))) +
  labs(x = "Covariance with r", y = NULL) +
  theme(
    panel.background    = element_blank(),
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.border        = element_blank(),
    axis.line.x         = element_line(color = "black", linewidth = 0.4),
    axis.line.y         = element_line(color = "black", linewidth = 0.4),
    axis.ticks          = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length   = unit(-0.15, "cm"),  # negative = inward ticks
    axis.text.x         = element_text(size = 12),
    axis.text.y         = element_text(size = 12, face = "italic"),
    axis.title.x        = element_text(size = 14),
    axis.title.y        = element_blank(),
    plot.margin         = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  )

final_fig <- herit_fig + inset_element(cov_simple, left = 0.02, bottom = 0.4, right = 0.7, top = 0.95)

## Appendix table of shape parameter moments:

descript_df <- TPC_summary_spread_r %>% 
  summarize(CT_min_mean=mean(CT_min),
            CT_min_SE=sd(CT_min)/sqrt(20),
            E_a_mean=mean(E_a),
            E_a_SE=sd(E_a)/sqrt(20),
            T_opt_mean=mean(T_opt, na.rm=TRUE),
            T_opt_SE=sd(T_opt, na.rm=TRUE)/sqrt(20),
            r_peak_mean=mean(r_peak),
            r_peak_SE=sd(r_peak)/sqrt(20)
  )


############################################################################################
### FIGURE 4A, 4B, 4C, 4D ------------------

# GENETIC CORRELATION PLOT BETWEEN SHAPE PARAMETERS
  
# Here we only keep relevant parameters E_a, E_d, T_opt, CT_min
p.mat <- cor.mtest(TPC_summary_spread[,c(-1,-2,-7,-8, -9, -10)][-7,]) #We eliminate genotype B2192.III for this because T_top is NA
dev.new()
# generate correlation plot between all parameters
corrplot(cor(as.data.frame(TPC_summary_spread2)[,c(-1,-2,-7,-8, -9, -10)][-7,]),method="color", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat$p, sig.level = 0.05, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE) 


# QUANTIFYING SELECTION
TPC_eqn <- function(a, E_a, E_d, Th, Temperature, r_scale) {
  exp(a + (E_a / (8.6 * 10^-5)) * (1 / 298.15 - 1 / (Temperature + 273.15)) -
        log(1 + exp((E_d / (8.6 * 10^-5)) * (1 / Th - 1 / (Temperature + 273.15))))) - r_scale
}

# Filter out only necessary temperatures
TPC_predicted_a<-expand.grid(Clone=Clones, Temperature= c(13.5000000,16.5000000,19.5000000, 22.0000000,25.0000000,28.0000000,31.0000000,34.0000000,37.0000000)) %>%
  left_join(TPC_summary_spread) %>%
  #left_join(dplyr::select(TPC_fits, -TPC_fit)) %>%
  #left_join(distinct(dplyr::select(tpcs, Clone))) %>%
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
  filter(T_opt.x>32) %>% #filtering outliers (B2192 III)
  select(Temp_range, r, E_a.x, T_opt.x, CT_min.x, r_peak.x) %>%
  mutate(E_a=E_a.x, T_opt=T_opt.x, r_peak=r_peak.x, CT_min=CT_min.x, .keep = "unused") %>%
  filter(is.finite(T_opt) & !is.na(T_opt)) %>%
  relocate(E_a, r_peak, CT_min,T_opt, .after=r) 
Sel_fig <- cbind(Sel_fig[,-c(3:11)], data.frame(lapply(Sel_fig[,3:6],scale)))

TPC_summary_spread %>% print(n=22)

## Generate models for each parameter of interest
T_opt_mod <- lm(formula = r ~ (T_opt+I(T_opt^2))*Temp_range , data = Sel_fig)
E_a_mod <- lm(formula = r ~ (E_a+I(E_a^2))*Temp_range , data = Sel_fig)
CT_min_mod <- lm(formula = r ~ (CT_min+I(CT_min^2))*Temp_range , data = Sel_fig)
r_peak_mod <- lm(formula = r ~ (r_peak+I(r_peak^2))*Temp_range , data = Sel_fig)

summary(T_opt_mod)
summary(E_a_mod)
summary(CT_min_mod)
summary(r_peak_mod)

# Predict r for each parameter of interest across temperature interests
pred_func<-function(Var, Temp){
  x_vals<-seq(min(Sel_fig_gathered$value, na.rm=T), max(Sel_fig_gathered$value, na.rm=T), length.out=200)
  temp_mod<-eval(parse(text=paste0(Var, "_mod")))
  temp_df<-data.frame(Variable=x_vals, Temp_range=Temp)
  colnames(temp_df)[1]<-Var
  temp_predict<-predict(temp_mod, temp_df, se=T)
  return(data.frame(Variable=x_vals, Value=temp_predict$fit, variables=Var, Temp_range=Temp))
}

Sel_fig_gathered <- Sel_fig %>%
  pivot_longer(cols= E_a:T_opt, names_to='variables') %>%
  mutate(variables=factor(variables, levels=c("r_peak", "E_a", "CT_min", "T_opt")))

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
  geom_hline(yintercept=0, color="darkgray", linewidth=0.7, linetype=2) +
  geom_point(data=Sel_fig_gathered, aes(x=value, y=r, color=Temp_range), size=2, shape=16) +
  geom_line(data=pred_df, aes(Variable, Value, color=Temp_range), linewidth=1) +
  scale_color_manual(values=c("orangered", "skyblue", "orange")) +
  facet_wrap(. ~ variables, scales = "free", ncol=2) +
  scale_y_continuous(limits=c(-5,13)) +
  labs(x="Variable", y=Intrinsic~growth~rate~(r)~(d^-1)) +
  theme(plot.background=element_blank(), panel.background=element_blank(),
        panel.border=element_rect(color = "black", size=1, fill=NA),
        axis.text=element_text(size=16), axis.title=element_text(size=16),
        axis.title.x.top=element_blank(), axis.text.x.top=element_blank(),
        axis.ticks.length.x.top=unit(0, "cm"),
        axis.ticks.length=unit(-0.15, "cm"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        strip.text = element_text(size=11),
        aspect.ratio=0.9,
        legend.key=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12, face="italic")
  ) 



###################################################################################################
### FIGURE 3E, 3F, 3G ----------


# PREDICTED EVOLUTION
#G-matrix estimate
#creating reduced dataset for analyses
TPC_summary_spread_r<- TPC_summary_spread %>% select(Clone,CT_min, E_a, T_opt, r_peak)
TPC_summary_spread_r<-TPC_summary_spread_r[-c(7, 15),] #filtering outliers (B2192 III)
# RESULTS IN MAIN TEXT CONTAIN THESE OUTLIERS, THEY DO NOT ALTER THE RESULTS

#TPC_summary_spread_r<-TPC_summary_spread_r[-c(2),]

print(TPC_summary_spread_r, n=22)

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
    joined$W<-log(exp(joined$r))
    #joined$W<-joined$W/mean(joined$W)
    # joined[,-1]<-sweep(joined[,-1],MARGIN = 2,STATS = colMeans(joined[,-1]),FUN = "/")
    out<-lm(cbind(CT_min,E_a,T_opt,r_peak,W)~1, data=joined) %>%
      evolqg::BayesianCalculateMatrix(.,samples=1000)
    out$Ps*0.5
  }
names(Gwdist)<-c(13,19, 22, 25, 30, 32, 38)

traits<-c("r_peak","E_a", "CT_min", "T_opt")
colTraits<-RColorBrewer::brewer.pal(4, "Accent")

deltaz.plot2 <- deltaz.plot
#Figure 3G, plotting the evolutionary responses according to Price equation
deltaz.plot <- (
  llply(Gwdist, function(x)  {
    foreach(i = 1:1000, .combine = "rbind") %do% {
      dz <- x[i, traits, "W"]
      G  <- x[i, traits, traits]
      dz * sqrt(diag(G))
    }
  }) %>%
    melt %>%
    mutate(Var2 = factor(Var2, traits)) %>%
    filter(
      (Var2 == "r_peak" & value > -2.2 & value < 5) |
      (Var2 == "E_a"    & value > -0.1 & value < 0.12) |
      (Var2 == "CT_min" & value > -5 & value < 5) |
      (Var2 == "T_opt"  & value > -3 & value < 2.5)
    ) %>%
    ggplot(aes(value, L1)) +
    geom_vline(xintercept = 0, linetype = 2) +
    ggridges::geom_density_ridges(aes(fill = Var2), alpha = 0.4, show.legend = FALSE,
                                  quantile_lines = TRUE, quantiles = c(0.025, 0.975)) +
    facet_grid(. ~ Var2, scales = "free", labeller = label_parsed) +
    scale_fill_manual(name = "Variable", values = colTraits) +
    ylim(names(Gwdist)) +
    ylab(expression(paste(Temp, "(", C^o, ")"))) +
    xlab(expression(paste("Predicted evolutionary change (", Delta, "z)"))) +
    theme(
      axis.ticks.length = unit(-0.15, "cm"),
      strip.text.x = element_blank(),
      strip.background = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank()
    )
)

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
        axis.ticks.length=unit(-0.15, "cm"),
        panel.background = element_blank(),       # remove gray background
    	plot.background = element_blank(),        # remove plot background
    	panel.grid = element_blank()              # remove grid lines
)

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
  
# Custom theme matching your target aesthetic
# Theme to match aesthetic style
custom_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.ticks = element_line(size = 0.5),
    axis.line = element_line(size = 0.6, color = "black"),
    plot.margin = margin(5, 5, 5, 5),
    plot.title = element_text(face = "bold", size = 16, hjust = -0.1, vjust = 1.5)
  )

# Function to return rho and p-value string
get_rho_p_label <- function(x, y) {
  test <- cor.test(x, y, method = "pearson")
  rho <- unname(test$estimate)
  pval <- test$p.value
  
  rho_str <- paste0("rho = ", format(round(rho, 2), nsmall = 2))  
  paste(rho_str, sep = "\n")
}

# Color palette (if needed)
colTraits <- brewer.pal(4, "Accent")

# Remove outliers already identified
#TPC_summary_spread_r <- TPC_summary_spread_r %>% 
#						filter(Clone!="SB3539.1")

# Plot 1: E_a vs CT_min
# Plot 1: E_a vs CT_min
label1 <- get_rho_p_label(TPC_summary_spread_r$E_a, TPC_summary_spread_r$CT_min)
ea_ctmin.plot <- ggplot(TPC_summary_spread_r, aes(scale(E_a), scale(CT_min))) + 
  stat_ellipse(geom = "polygon", fill = "gray30", alpha = 0.2) +
  geom_point(aes(color = Clone), size = 3, show.legend = FALSE) +
  xlab(expression(E[a])) +
  ylab(expression(CT[min])) +
  annotate("text", x = Inf, y = -Inf, label = label1, hjust = 1.1, vjust = -0.5, size = 5) +
  ggtitle("e") +
  custom_theme +
  theme(axis.ticks.length = unit(-0.15, "cm"))

# Plot 2: T_opt vs CT_min
label2 <- get_rho_p_label(TPC_summary_spread_r$r_peak, TPC_summary_spread_r$CT_min)
topt_ctmin.plot <- ggplot(TPC_summary_spread_r, aes(scale(r_peak), scale(CT_min))) +
  stat_ellipse(geom = "polygon", fill = "gray30", alpha = 0.2) +
  geom_point(aes(color = Clone), size = 3, show.legend = FALSE) +
  xlab(expression(r[peak])) +
  ylab(expression(CT[min])) +
  annotate("text", x = Inf, y = -Inf, label = label2, hjust = 1.1, vjust = -0.5, size = 5) +
  ggtitle("f") +
  custom_theme +
  theme(axis.ticks.length = unit(-0.15, "cm"))

# Plot 3: E_a vs r_peak
label3 <- get_rho_p_label(TPC_summary_spread_r$E_a, TPC_summary_spread_r$r_peak)
ea_rpeak.plot <- ggplot(TPC_summary_spread_r, aes(scale(E_a), scale(r_peak))) +   
  stat_ellipse(geom = "polygon", fill = "gray30", alpha = 0.2) +
  geom_point(aes(color = Clone), size = 3, show.legend = FALSE) +
  xlab(expression(E[a])) +
  ylab(expression(r[peak])) +
  annotate("text", x = Inf, y = -Inf, label = label3, hjust = 1.1, vjust = -0.5, size = 5) +
  ggtitle("g") +
  custom_theme +
  theme(axis.ticks.length = unit(-0.15, "cm"))

Fig3_B <- (ea_rpeak.plot+ea_ctmin.plot+topt_ctmin.plot)/beta.plot/deltaz.plot


##########################################################################################
### FIGURE 4D, 4E, 4F -------------------------------------

# Data processing 
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
  mutate(CU = total - adjCount) %>% 
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

## Appendix Table
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

library("dplyr")
library("nls.multstart")
library("tidyr")
library("ggplot2")

## 1) Fit TPCs to the TPCs for clones AXS and CU416
tpc_exp_pred <- tpc_exp_pred %>% mutate( r_scale=2,log_r=log(r+r_scale))
tpc_exp_pred$Clone[tpc_exp_pred$Clone =="CU416"] <- "CU4106"

TPC_eqn<-function(a, E_a, E_d, Th, Temperature, r_scale){exp(a + (E_a/(8.6*10^-5))*(1/298.15-1/(Temperature+273.15)) - log(1+exp((E_d/(8.6*10^-5))*(1/Th-1/(Temperature+273.15)))))-r_scale}

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
         #T_opt=E_d*Th/(E_d+8.6e-5*Th*log(E_d/E_a-1))
         )

pars <- dplyr::select(TPC_fits_exp, -TPC_fit)  

TPC_predicted_exp<-expand.grid(Clone=c("AXS", "CU4106"), Temperature=seq(0, 50, length.out=500)) %>%
  left_join(dplyr::select(TPC_fits_exp, -TPC_fit)) %>%
  left_join(distinct(dplyr::select(tpc_exp_pred, Clone, r_scale))) %>%
  mutate(r=TPC_eqn(a, E_a, E_d, Th, Temperature, r_scale)) #%>%
#filter(r>=-1)

# plot all species' TPCs together
TPC_predicted_exp$Clone<-factor(TPC_predicted_exp$Clone)
head(TPC_predicted_exp)

#2) Use the TPCs to predict genotypic frequencies across temperatures

# Calculate mean_fitness and prepare data 
mean_fitness <- TPC_predicted_exp %>%
  dplyr::group_by(Temperature) %>%               
  dplyr::summarize(meanFit=mean(r), SD=sd(r), dissim=mean(dist(r))) 

TPC_predicted_exp <- merge(TPC_predicted_exp,mean_fitness,by="Temperature")
TPC_predicted_exp$rel_fitness <- (TPC_predicted_exp$r+10)/(TPC_predicted_exp$meanFit+10)
Temp <- unique(TPC_predicted_exp$Temperature)

# Plot
g_exp<-
  ggplot()+
  geom_hline(yintercept=0, color="gray", linewidth=0.7, linetype=2)+
  geom_line(data=TPC_predicted_exp, aes(Temperature, r, color=Clone), linewidth=1)+
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

## Define function for model
## Temperature has to be passed as Temp[nmbr] because it has to match existing temperatures in the TPC_predict dataset that contains the TPCs that will be used to make the average fitness predictions.
# Row sampler
sample_from_vector <- function(values, i, target_variance, n) {
  # If variance is zero, return the row in i
  if(target_variance==0){return(rep(values[i],n))}
  else{
    # Otherwise, extract the mean and calculate the variance based on neighbors
    mean_i <- values[i]
    current_variance <- mean(c((values[i - 1] - mean_i)^2, (values[i + 1] - mean_i)^2))
    # Calculate the scaling factor for the target variance
    scaling_factor <- sqrt(target_variance / current_variance)
    # Calculate weights for sampling based on proximity to mean
    deviations <- (values - mean_i)^2
    weights <- exp(-deviations / (2 * scaling_factor^2))
    weights <- weights / sum(weights) # Normalize weights
    # Sample n values with the computed weights
    sampled_indices <- sample(seq_along(values), size = n, replace = TRUE, prob = weights)
    # Return the sampled values
    return(values[sampled_indices])
  }
}

#### DEFINE MODEL
evo_mod <- function(temp,variance,time_steps){
  ## Select appropriate temperature data
  row <- which(Temp==temp)[1]
  if(row!=500){
    # Sample vector of random temperature fluctuations
    r_t <- sample_from_vector(Temp, i = row, target_variance = variance, n = time_steps)
  }else{
    # Because sampler cannot sample at i=500
    r_t <- sample_from_vector(Temp, i = 499, target_variance = variance, n = time_steps)	
  }
  
  ## Run model
  N_clones <- 2
  time <- seq(1,time_steps,1)
  mat_freq <- matrix(rep(0,N_clones*time_steps),nrow=N_clones,ncol=time_steps)
  mat_freq[,1] <- rep(1/N_clones,N_clones)	
  for(i in 1:(time_steps-1)){
    if(i==1){ #We add +4 to r to avoid dealing with negative vaues of r and how it impacts the calculation of rel.fit.
      # Find temperature after fluctuation
      fit <- TPC_predicted_exp %>%
        dplyr::filter(Temperature==r_t[i]) %>%
        dplyr::select(r)+2 #To ensure no negative values are possible
      # Find new frequencies
      new_freq <- diag(fit$r/sum((fit$r)*mat_freq[,1]),N_clones,N_clones)%*%rep(1/N_clones,N_clones) 
      mat_freq[,2] <- new_freq 
    }else{
      fit <- TPC_predicted_exp %>%
        dplyr::filter(Temperature==r_t[i]) %>%
        dplyr::select(r)+4 #To ensure no negative values are possible
      # Find new frequencies
      mat_freq[,i+1] <-diag(fit$r/sum((fit$r)*mat_freq[,i]),N_clones,N_clones)%*%mat_freq[,i]
    }
  }
  return(rowMeans(mat_freq[,1:time_steps]))
  #return(mat_freq[,1:time_steps]) 
  # If change for return(mat_freq) it is possible to plot the dynamics of the model against time
}

## Check that function is working
mat <- evo_mod(Temp[200],0.05,100)
dim(mat)
## This is good to see the dynamics of the model against time but function needs to be altered in last lines as explained in function body
matplot(seq(1,100),t(mat[,1:100]), type = 'l', lty=1) 
plot(c(1,1))

## Actual runs (takes a few seconds to run depending on how many temperatures)
runs <- sapply(Temp,evo_mod,0.01,100)
runs2 <- sapply(Temp,evo_mod,0.01,50)
runs3 <- sapply(Temp,evo_mod,0.01,20)

#Runs with no noise
runs4 <- sapply(Temp,evo_mod,0.00,100)
runs5 <- sapply(Temp,evo_mod,0.00,50)
runs6 <- sapply(Temp,evo_mod,0.00,20)

#Runs with more noise
runs7 <- sapply(Temp,evo_mod,0.02,100)
runs8 <- sapply(Temp,evo_mod,0.02,50)
runs9 <- sapply(Temp,evo_mod,0.02,20)
#matplot(Temp,t(runs), type = 'l', lty=1)

## Prep data for plotting
runs_clone <- as.data.frame(runs) %>% mutate(Clone=c("AXS", "CU4106")) %>% relocate("Clone")
colnames(runs_clone) <- c("Clone", Temp)
runs_clone2 <- as.data.frame(runs2) %>% mutate(Clone=c("AXS", "CU4106")) %>% relocate("Clone")
colnames(runs_clone2) <- c("Clone", Temp)
runs_clone3 <- as.data.frame(runs3) %>% mutate(Clone=c("AXS", "CU4106")) %>% relocate("Clone")
colnames(runs_clone3) <- c("Clone", Temp)

runs_clone4 <- as.data.frame(runs4) %>% mutate(Clone=c("AXS", "CU4106")) %>% relocate("Clone")
colnames(runs_clone4) <- c("Clone", Temp)
runs_clone5 <- as.data.frame(runs5) %>% mutate(Clone=c("AXS", "CU4106")) %>% relocate("Clone")
colnames(runs_clone5) <- c("Clone", Temp)
runs_clone6 <- as.data.frame(runs6) %>% mutate(Clone=c("AXS", "CU4106")) %>% relocate("Clone")
colnames(runs_clone6) <- c("Clone", Temp)

runs_clone7 <- as.data.frame(runs7) %>% mutate(Clone=c("AXS", "CU4106")) %>% relocate("Clone")
colnames(runs_clone7) <- c("Clone", Temp)
runs_clone8 <- as.data.frame(runs8) %>% mutate(Clone=c("AXS", "CU4106")) %>% relocate("Clone")
colnames(runs_clone8) <- c("Clone", Temp)
runs_clone9 <- as.data.frame(runs9) %>% mutate(Clone=c("AXS", "CU4106")) %>% relocate("Clone")
colnames(runs_clone9) <- c("Clone", Temp)


new_runs <- runs_clone %>%
  gather("Temp","Freq",2:501) %>%
  mutate(Temp=as.numeric(Temp)) 
new_runs2 <- runs_clone2 %>%
  gather("Temp","Freq",2:501) %>%
  mutate(Temp=as.numeric(Temp)) 
new_runs3 <- runs_clone3 %>%
  gather("Temp","Freq",2:501) %>%
  mutate(Temp=as.numeric(Temp)) 

new_runs4 <- runs_clone4 %>%
  gather("Temp","Freq",2:501) %>%
  mutate(Temp=as.numeric(Temp)) 
new_runs5 <- runs_clone5 %>%
  gather("Temp","Freq",2:501) %>%
  mutate(Temp=as.numeric(Temp)) 
new_runs6 <- runs_clone6 %>%
  gather("Temp","Freq",2:501) %>%
  mutate(Temp=as.numeric(Temp)) 

new_runs7 <- runs_clone7 %>%
  gather("Temp","Freq",2:501) %>%
  mutate(Temp=as.numeric(Temp)) 
new_runs8 <- runs_clone8 %>%
  gather("Temp","Freq",2:501) %>%
  mutate(Temp=as.numeric(Temp)) 
new_runs9 <- runs_clone9 %>%
  gather("Temp","Freq",2:501) %>%
  mutate(Temp=as.numeric(Temp)) 

## Plot

new_runs$Clone <- factor(new_runs$Clone, levels = rev(unique(new_runs$Clone)))
g_1 <- ggplot(new_runs, aes(x = Temp, y = Freq, fill = Clone)) + 
  geom_area(size = 0.1, colour = "white") + 
  scale_x_continuous(limits = c(10, 38)) + 
  scale_fill_manual(name = "Clone", values = c("AXS" = "#C3D48A", "CU4106" = "#de95c0")) + 
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    axis.title.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.length.x.top = unit(0, "cm"),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    aspect.ratio = 0.55,
    legend.position = "none"
  )

new_runs2$Clone <- factor(new_runs2$Clone, levels = rev(unique(new_runs2$Clone)))
g_2 <- ggplot(new_runs2, aes(x = Temp, y = Freq, fill = Clone)) + 
  geom_area(size = 0.1, colour = "white") + 
  scale_x_continuous(limits = c(10, 38)) + 
  scale_fill_manual(name = "Clone", values = c("AXS" = "#C3D48A", "CU4106" = "#de95c0")) + 
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    axis.title.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.length.x.top = unit(0, "cm"),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    aspect.ratio = 0.55,
    legend.position = "none"
  )

new_runs3$Clone <- factor(new_runs3$Clone, levels = rev(unique(new_runs3$Clone)))
g_3 <- ggplot(new_runs3, aes(x = Temp, y = Freq, fill = Clone)) + 
  geom_area(size = 0.1, colour = "white") + 
  scale_x_continuous(limits = c(10, 38)) + 
  scale_fill_manual(name = "Clone", values = c("AXS" = "#C3D48A", "CU4106" = "#de95c0")) + 
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    axis.title.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.length.x.top = unit(0, "cm"),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    aspect.ratio = 0.55,
    legend.position = "none"
  )

new_runs4$Clone <- factor(new_runs4$Clone, levels = rev(unique(new_runs4$Clone)))
g_4 <- ggplot(new_runs4, aes(x = Temp, y = Freq, fill = Clone)) + 
  geom_area(size = 0.1, colour = "white") + 
  scale_x_continuous(limits = c(10, 38)) + 
  scale_fill_manual(name = "Clone", values = c("AXS" = "#C3D48A", "CU4106" = "#de95c0")) + 
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    axis.title.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.length.x.top = unit(0, "cm"),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    aspect.ratio = 0.55,
    legend.position = "none"
  )

new_runs5$Clone <- factor(new_runs5$Clone, levels = rev(unique(new_runs5$Clone)))
g_5 <- ggplot(new_runs5, aes(x = Temp, y = Freq, fill = Clone)) + 
  geom_area(size = 0.1, colour = "white") + 
  scale_x_continuous(limits = c(10, 38)) + 
  scale_fill_manual(name = "Clone", values = c("AXS" = "#C3D48A", "CU4106" = "#de95c0")) + 
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    axis.title.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.length.x.top = unit(0, "cm"),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    aspect.ratio = 0.55,
    legend.position = "none"
  )

new_runs6$Clone <- factor(new_runs6$Clone, levels = rev(unique(new_runs6$Clone)))
g_6 <- ggplot(new_runs6, aes(x = Temp, y = Freq, fill = Clone)) + 
  geom_area(size = 0.1, colour = "white") + 
  scale_x_continuous(limits = c(10, 38)) + 
  scale_fill_manual(name = "Clone", values = c("AXS" = "#C3D48A", "CU4106" = "#de95c0")) + 
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    axis.title.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.length.x.top = unit(0, "cm"),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    aspect.ratio = 0.55,
    legend.position = "none"
  )

new_runs7$Clone <- factor(new_runs7$Clone, levels = rev(unique(new_runs7$Clone)))
g_7 <- ggplot(new_runs7, aes(x = Temp, y = Freq, fill = Clone)) + 
  geom_area(size = 0.1, colour = "white") + 
  scale_x_continuous(limits = c(10, 38)) + 
  scale_fill_manual(name = "Clone", values = c("AXS" = "#C3D48A", "CU4106" = "#de95c0")) + 
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    axis.title.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.length.x.top = unit(0, "cm"),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    aspect.ratio = 0.55,
    legend.position = "none"
  )

new_runs8$Clone <- factor(new_runs8$Clone, levels = rev(unique(new_runs8$Clone)))
g_8 <- ggplot(new_runs8, aes(x = Temp, y = Freq, fill = Clone)) + 
  geom_area(size = 0.1, colour = "white") + 
  scale_x_continuous(limits = c(10, 38)) + 
  scale_fill_manual(name = "Clone", values = c("AXS" = "#C3D48A", "CU4106" = "#de95c0")) + 
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    axis.title.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.length.x.top = unit(0, "cm"),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    aspect.ratio = 0.55,
    legend.position = "none"
  )

new_runs9$Clone <- factor(new_runs9$Clone, levels = rev(unique(new_runs9$Clone)))
g_9 <- ggplot(new_runs9, aes(x = Temp, y = Freq, fill = Clone)) + 
  geom_area(size = 0.1, colour = "white") + 
  scale_x_continuous(limits = c(10, 38)) + 
  scale_fill_manual(name = "Clone", values = c("AXS" = "#C3D48A", "CU4106" = "#de95c0")) + 
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    axis.title.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.length.x.top = unit(0, "cm"),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    aspect.ratio = 0.55,
    legend.position = "none"
  )


Fig4appendix_timeflux <- grid.arrange(
             g_6, g_5, g_4, 
             g_3, g_2, g_1, 
             g_9, g_8, g_7, 
             nrow = 3, ncol = 3,
             padding = unit(0.1, "cm"))

## 
# Find indexes for the experimental temperatures used
close_index <- function(ref_temp, Temp){
  which.min(abs(ref_temp - Temp))
}
index_list <- sapply(c(19,22,25,30,32,38),close_index, new_runs$Temp)

new_runs_select <- rbind(new_runs4[index_list,],new_runs4[index_list+1,]) %>% 
						mutate(Temp=as.factor(round(Temp))) 
new_runs_select$Freq[1]
# Apply control relative to 19ºC (i.e., clearly the quantitative pattern driven by more then the TPCs,
# so we subtract the difference between AXS at 19º observed and predicted across the board but present
# both in the figure )	
Freq_cont <- rbind(
	(new_runs_select$Freq[which(new_runs_select$Clone=="AXS")]-new_runs_select$Freq[1]+10^-6) %>% t() %>% t(),
	1-(new_runs_select$Freq[which(new_runs_select$Clone=="AXS")]-new_runs_select$Freq[1]+10^-6) %>% t() %>% t()
	)

new_runs_select$Freq_cont <- Freq_cont
  						
## Last figure of main text.
## Set order of factors first
new_runs_select$Clone <- factor(
  new_runs_select$Clone,
  levels = c("CU4106","AXS")
)

g_11 <- ggplot(new_runs_select, aes(x = factor(Temp),   # treat Temp as discrete
                                    y = Freq_cont,
                                    fill = Clone)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_errorbar(
    data = filter(new_runs_select, Clone == "AXS"),
    aes(x = factor(Temp), ymin = Freq, ymax = Freq),    
    width = .65,                                        
    colour = "black",
    linewidth = .7,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = "Clone",
    values = c("AXS" = "#C3D48A", "CU4106" = "#de95c0")
  ) +
  labs(x = "Temperature (ºC)", y = "Pred. Frequency") +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    aspect.ratio = 0.55,
    legend.position = "right",
    legend.text = element_text(size = 12)
  )


