##################################
# Proximal and Distal Model Code #
##################################

# This file performs the proximal and distal model analysis for a hybrid SMART-MRT. It uses the following inputs:
# proximal_w_resp: data for the proximal model including responders and nonresponders
# distal_nonresp: data for the distal model including only nonresponders
# Tables and figures are numbered to match those in the manuscript
# All output tables are formatted using Kable

# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "rootSolve", "kableExtra", "knitr", "geepack", "data.table")

# Source for simulated data
proximal_w_resp = read.csv("sim_data_for_proximal.csv")
distal_nonresp = read.csv("sim_data_for_distal_nonresponders_only.csv")

# edit names
names(proximal_w_resp) = c("id", "Biological Sex (Mean Centered)", "Baseline BMI (Mean Centered)", "Z1", "Responder", "Week Classified as Nonresponder", "Z2", "A", "Days Since Classified as Nonresponder", "Outcome", "Probability of Microrandomization")
names(distal_nonresp) = c("id", "Biological Sex (Mean Centered)", "Baseline BMI (Mean Centered)", "Mean A", "A (Mean Centered)", "Z1", "Z2", "Outcome")

# Source for estimator
source("Estimator_For_Proximal_Analysis.R")

#####
# Proximal Model
#####
# Treat ID as a nominal variable
proximal_w_resp$id = factor(proximal_w_resp$id)
# Add variable for interaction of Z1 and Z2
proximal_w_resp1 = proximal_w_resp %>%
  mutate(Z1_Z2 = Z1*Z2) 

# Control variable regression for proximal model
resp_model_control = geeglm(Outcome ~ Z1 + Z2 + `Week Classified as Nonresponder` + 
                               `Biological Sex (Mean Centered)` + `Days Since Classified as Nonresponder` + `Baseline BMI (Mean Centered)`, 
                             data = proximal_w_resp, id = id)
print(summary(resp_model_control))

# Control variable output table
resp_model_control_table = data.table(tidy(resp_model_control))
cilow = data.table(resp_model_control_table$estimate - 1.96 * resp_model_control_table$std.error)
cihigh = data.table(resp_model_control_table$estimate + 1.96 * resp_model_control_table$std.error)
resp_model_control_table = data.frame(cbind(resp_model_control_table[,2:3], cilow, cihigh, resp_model_control_table[,4:5]))
names(resp_model_control_table) = c("Estimate", "Robust SE", "95% CI LL", "95% CI UL", "Wald", "Pr>|W|")
rownames(resp_model_control_table) = c("Intercept", "$Z_{i1}$", "$Z_{i2}$", "Week Rerandomized", "Sex", "Days Since Rerandomization", "Baseline BMI (Mean Centered)")
resp_model_control_table = round(resp_model_control_table, digits = 4)
print(resp_model_control_table)
resp_model_control_table %>% 
  kbl(caption = "Table 8.1: Control Variables") %>%
  kable_classic(full_width = F, html_font = "Times New Roman") %>%
  footnote(general = "CI: confidence interval; LL: lower limit; UL: upper limit")

# Proximal model
proximal = binary_outcome_moderated_effect_resp(
  dta = proximal_w_resp1,
  control_var = c("Z1", "Z2", "Week Classified as Nonresponder", "Biological Sex (Mean Centered)", 
                  "Days Since Classified as Nonresponder", "Baseline BMI (Mean Centered)"),
  moderator = c("Z1", "Z2", "Z1_Z2"),
  id_var = "id",
  day_var = "Days Since Classified as Nonresponder",
  trt_var = "A",
  outcome_var = "Outcome",
  avail_var = NULL,
  prob_treatment = "Probability of Microrandomization",
  significance_level = 0.05)

# Proximal model output table
proximal_model = as.data.frame(proximal[1])
names(proximal_model) = c("Estimate", "Robust SE", "95% CI LL", "95% CI UL", "t", "Pr>|t|")
rownames(proximal_model) = c("$A_{it}$", 
                              "$A_{it} Z_{i1}$",
                              "$A_{it} Z_{i2}$", 
                              "$A_{it} Z_{i1} Z_{i2}$")
proximal_round = round(proximal_model, digits = 4)
print(proximal_round)
proximal_round %>% 
  kbl(caption = "Table 8.2: Proximal Model Results") %>%
  kable_classic(full_width = F, html_font = "Times New Roman") %>% 
  footnote(general = "CI: confidence interval; LL: lower limit; UL: upper limit")

#####
# Distal Model
#####

##### Distal Model : nonresponders only #####

# Distal model 
model_distal = geeglm(Outcome ~ Z1 * Z2 * `A (Mean Centered)` + 
                         `Biological Sex (Mean Centered)` + `Baseline BMI (Mean Centered)`, 
                       data = distal_nonresp, id = id)
print(summary(model_distal))

# Distal model output table
model_distal_table = data.table(tidy(model_distal))
cilow = data.table(model_distal_table$estimate - 1.96 * model_distal_table$std.error)
cihigh = data.table(model_distal_table$estimate + 1.96 * model_distal_table$std.error)
model_distal_table = data.frame(cbind(model_distal_table[,2:3], cilow, cihigh, model_distal_table[,4:5]))
names(model_distal_table) = c("Estimate", "Robust SE", "95% CI LL", "95% CI UL", "Wald", "Pr>|W|")
rownames(model_distal_table) = c("Intercept", "$Z_{i1}$", "$Z_{i2}$", "$\\bar{A_i}^{(2)}$", "Biological Sex (Mean Centered)", "Baseline BMI (Mean Centered)", "$Z_{i1} Z_{i2}$", "$Z_{i1} \\bar{A_i}^{(2)}$", "$Z_{i2}\\bar{A_i}^{(2)}$", "$Z_{i1}Z_{i2}\\bar{A_i}^{(2)}$")
model_distal_table = round(model_distal_table, digits = 4)
model_distal_table %>% 
  kbl(caption = "Table 9: Results (Parameter Estimates) for Distal Model") %>% 
  kable_classic(full_width = F, html_font = "Times New Roman") %>% 
  footnote(general = "CI: confidence interval; LL: lower limit; UL: upper limit; all covariates centered")

##### Distal Plot #####

beta = model_distal$coefficients;
data_to_plot = expand.grid(Z1=c(-1,+1),
                            Z2=c(-1,+1),
                            rate=c("mean","below","above"));
data_to_plot$regimen = paste("(",
                              data_to_plot$Z1,
                              ",",
                              data_to_plot$Z2,
                              ")",
                              sep="");
the_mean = round(mean(distal_nonresp$`Mean A`),4)
the_sd = round(sd(distal_nonresp$`A (Mean Centered)`),4);
data_to_plot$`A (Mean Centered)`[which(data_to_plot$rate=="mean")] = 0;
data_to_plot$`A (Mean Centered)`[which(data_to_plot$rate=="above")] = 0 + the_sd;
data_to_plot$`A (Mean Centered)`[which(data_to_plot$rate=="below")] = 0 - the_sd;
data_to_plot$Z1 = as.integer(data_to_plot$Z1);
data_to_plot$Z2 = as.integer(data_to_plot$Z2); 
beta = model_distal$coefficients;
data_to_plot$fitted_value = beta["(Intercept)"] + 
  beta["Z1"]*data_to_plot$Z1 + 
  beta["Z2"]*data_to_plot$Z2 + 
  beta["`A (Mean Centered)`"]*data_to_plot$`A (Mean Centered)` + 
  beta["Z1:Z2"]*data_to_plot$Z1*data_to_plot$Z2 + 
  beta["Z1:`A (Mean Centered)`"]*data_to_plot$Z1*data_to_plot$`A (Mean Centered)` + 
  beta["Z2:`A (Mean Centered)`"]*data_to_plot$Z2*data_to_plot$`A (Mean Centered)` + 
  beta["Z1:Z2:`A (Mean Centered)`"]*data_to_plot$Z1*data_to_plot$Z2*data_to_plot$`A (Mean Centered)` ;
print(data_to_plot);

plot1 = ggplot(data_to_plot, 
               aes(x = factor(regimen),y=fitted_value)) +
  geom_point(aes(y = fitted_value,
                 colour=as.factor(rate),
                 shape=as.factor(rate)),
             size = 2) + 
  scale_colour_manual(name = 'Rate of message delivery', 
                      values =c('red', 'blue', "green"),
                      labels = c(paste('Mean (',the_mean,')'),
                                 paste('1 SD Below Mean (',the_mean-the_sd,')'), 
                                 paste('1 SD Above Mean (',the_mean+the_sd,')'))) +
  scale_shape_manual(name = 'Rate of message delivery', 
                     values =c(15, 17, 16),
                     labels = c(paste('Mean (',the_mean,')'),
                                paste('1 SD Below Mean (',the_mean-the_sd,')'), 
                                paste('1 SD Above Mean (',the_mean+the_sd,')'))) +
  labs(title = str_wrap("Figure 4: Estimated weight loss by the initial options, subsequent options and rate of message delivery", 60), 
       x = "Sequence of initial and subsequent options", 
       y = "Estimated weight loss") ;
plot(plot1)
