library(mrgsolve)
log((0.778^2)+1)
log((0.261^2)+1)
log((0.801^2)+1)

VF2 <- 0.366
eta_VF2 <- 0.4733007
PHI_2i <- log(VF2/(1-VF2));
VF2_2i <- exp(PHI_2i + eta_VF2 ) /(1 + exp(PHI_2i + eta_VF2));

model <- '
$PARAM  
KA = 0.464,    // Absorption rate (1/h)  
ALAG2 = 2.30,  // Lag time 2nd dose (h)  
CL = 19.1,     // Clearance (L/h)  
V3 = 79.4,     // Central volume (L)  
Q = 12.0,      // Intercompartmental clearance (L/h)  
V4 = 199,      // Peripheral volume (L)  
Ftot = 0.0793, // Bioavailability  
VF2 = 0.366,   // Fraction to 2nd depot  

$OMEGA @name IIV @labels eta_VF2 eta_CL eta_Ftot 
0.4733007  // IIV_VF2 
0.06590103  // IIV_CL  
0.495672  // IIV_Ftot  

$SIGMA  
0.01 
0.00001  

$CMT DEPOT1 DEPOT2 CENTRAL PERIPH 

$MAIN
double CLi = CL * exp(eta_CL); 
double V3i = V3; 
double Qi = Q;
double V4i = V4;
double KAi = KA;
double PHI_2i = log(VF2/(1-VF2));
double VF2_2i = exp(PHI_2i + eta_VF2 ) /(1 + exp(PHI_2i + eta_VF2));
double F2i = VF2_2i * Ftot * exp(eta_Ftot);
double VF1i = (1 - VF2_2i);
double F1i = VF1i * Ftot * exp(eta_Ftot); 
F_DEPOT1 = F1i;
F_DEPOT2 = F2i;
ALAG_DEPOT2 = ALAG2;

$ODE  
dxdt_DEPOT1 = -KAi * DEPOT1;  
dxdt_DEPOT2 = -KAi * DEPOT2;  
dxdt_CENTRAL = KAi*DEPOT1 + KAi*DEPOT2 - (CLi/V3i)*CENTRAL - (Qi/V3i)*CENTRAL + (Qi/V4i)*PERIPH;  
dxdt_PERIPH = (Qi/V3i)*CENTRAL - (Qi/V4i)*PERIPH;

$TABLE
double CP = CENTRAL / V3i;
double DV = (CP * (1 + EPS(1)))+ EPS(2);

$CAPTURE CP DV CLi V3i Qi V4i F1i F2i;
'

mod <- mcode("final_model", model)




set.seed(4040)
dose_events1 <- ev(cmt = 1, ID = 1:2500, amt = 5000, ii = 12, addl = 0)
dose_events2 <- ev(cmt = 2, ID = 1:2500, amt = 5000, ii = 12, addl = 0)
dose_events <- dose_events1 + dose_events2
dose1 <- as_tibble(dose_events) %>%
  select(ID, time, amt, ii, addl, evid, cmt) %>%
  mutate(dose_group = "5mg ") %>%
  arrange(ID, time)

set.seed(4041)
dose_events3 <- ev(cmt = 1, ID = 2501:5000, amt = 2000, ii = 12, addl = 0)
dose_events4 <- ev(cmt = 2, ID = 2501:5000, amt = 2000, ii = 12, addl = 0)
dose_events <- dose_events3 + dose_events4
dose2 <- as_tibble(dose_events) %>%
  select(ID, time, amt, ii, addl, evid, cmt) %>%
  mutate(dose_group = "2mg ") %>%
  arrange(ID, time)

dose <- rbind(dose1,dose2)
#model simulation
set.seed(2)

out <- mod %>%
  data_set(dose) %>%
  Req(CP, DV, CLi, V3i, Qi, V4i, F1i, F2i) %>%
  mrgsim(end = 12, delta = 0.5)


simul_rsv <- as_tibble(out) 

aucs <- simul_rsv %>%
  group_by(ID) %>%
  mutate(AUCt = trapz(time, DV))

 dadosaa <- aucs %>% filter(ID == c("6","2","3","1","5"))
 ggplot(dadosaa, aes(x = time, y = DV, group = ID)) +
   geom_line(linewidth = 1.1, position = position_dodge(width=0.05))

simul_rsv3 <- aucs %>% filter(time==0.5 | time==1 |time==1.5 | time==2|time==3 | time==4| time== 6| time== 8)
pivot_conc <- simul_rsv3 %>% 
  select(ID, time, DV) %>% 
  pivot_wider(id_cols = ID, names_from = time, names_prefix = "conc_", values_from = DV)

base_prediction_auc <- simul_rsv3 %>% 
  filter(time==0.5) %>% 
  select( -time, -CP) %>% left_join(pivot_conc)

#Remove <1 >99 percentiles
centile_1_auc <- quantile(base_prediction_auc$AUCt, 0.01)
centile_99_auc <- quantile(base_prediction_auc$AUCt, 0.99)

base_prediction_auc1 <- base_prediction_auc %>%
  filter(AUCt >= centile_1_auc, AUCt <= centile_99_auc)

 fwrite(base_prediction_auc1, file = "rsvsim202112h.csv")

setwd("C:/Users/Matheus/Desktop/R")
data <- as.data.frame(fread("rsvsim2021.csv"))
df1 <- as.data.frame(fread("rsvalef.csv"))
# df1 <- df1 %>% filter(!time %in% c("16", "24"))
df1 <- df1 %>%
  group_by(pacientes, phase) %>%
  mutate(AUCt = trapz(time, conc))

# df1 <- df1 %>% mutate(time = case_when(
#   time == 0 ~ "C0",
#   time == 0.5 ~ "C05",
#   time == 1 ~ "C1",
#   time == 1.5 ~ "C15",
#   time == 2 ~ "C2",
#   time == 3 ~ "C3",
#   time == 4 ~ "C4",
#   time == 6 ~ "C6",
#   time == 8 ~ "C8",
#   time == 12 ~ "C12",
#   TRUE ~ as.character(time)
# ))

df1 <- df1 %>% mutate(time = case_when(
  time == 0 ~ "C0",
  time == 0.5 ~ "C05",
  time == 1 ~ "C1",
  time == 1.5 ~ "C15",
  time == 2 ~ "C15",
  time == 3 ~ "C3",
  time == 4 ~ "C4",
  time == 6 ~ "C6",
  time == 8 ~ "C8",
  time == 12 ~ "C12",
  TRUE ~ as.character(time)
))

pivot_df <- df1 %>% pivot_wider(id_cols = c(pacientes,phase), names_from = time, values_from = conc)

df1 <- df1 %>% 
  filter(time=="C05") %>% 
  select( -time, -conc) %>% left_join(pivot_df)
df1 <- as.data.frame(df1)


df2 <- as.data.frame(fread("rsvleandro.csv"))

df2 <- df2 %>% filter(!time %in% c("48")) %>% 
  group_by(ID, group,phase) %>%
  mutate(AUCt = trapz(time, conc))
df2 <- df2 %>% filter(phase == "Phase 2")

df2 <- df2 %>% mutate(time = case_when(
  time == 0 ~ "C0",
  time == 0.5 ~ "C05",
  time == 1 ~ "C1",
  time == 1.5 ~ "C15",
  time == 2 ~ "C2",
  time == 3 ~ "C3",
  time == 4 ~ "C4",
  time == 6 ~ "C6",
  time == 8 ~ "C8",
  time == 12 ~ "C12",
  TRUE ~ as.character(time)
))

pivot_df <- df2 %>% pivot_wider(id_cols = c(ID, group,phase), names_from = time, values_from = conc)

df2 <- df2 %>% 
  filter(time=="C05") %>% 
  select( -time, -conc) %>% left_join(pivot_df)
df2 <- as.data.frame(df2)

# df3 <- df2 %>% filter(group == "Group 1")

data <- data %>% rename(C05 = conc_0.5,C1 = conc_1,C15 = conc_1.5, C2 = conc_2, C3 = conc_3,
         C4= conc_4,C6 = conc_6,C8 = conc_8)

n_cores <- detectCores() - 2
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)

results <- test_combinations(
  data = data,
  # time_predictors = c("C05","C1","C15", "C2","C3", "C4","C6","C8"),
   time_predictors = c("C15","C4"),
  fixed_predictors = c(),
  response = "AUCt",
  external_dfs = list(
    "alef" = df1,
    "leandrof2" = df2
  ),
  combo_sizes = c(2),
  include_fixed_predictors = F,
  include_ratio_predictors = F,
  include_diff_predictors = F,
  return_fits = T,
  return_workflows = T,
  # algorithms = c( "xgb","svm","rf","glmnet"))
  algorithms = c( "rf"))
final_results <- bind_rows(results$all_tabs)
fwrite(wafinal_results, "resultados_completosrsv24h2arrumado.csv")

# wafinal_results24 <- as.data.frame(fread("resultados_completosrsv.csv"))

final_data <- data %>% select(C15,C4,AUCt)

# Training/testing split with stratification
set.seed(123)
final_split <- initial_split(final_data, strata = AUCt, prop = 0.75)
final_train <- training(final_split)
final_test <- testing(final_split)

set.seed(123)

rsv_split <- initial_split(data, 
                           strata = AUCt, 
                           prop= 0.75)
rsv_ml_train  <- training(rsv_split )
rsv_ml_test  <- testing(rsv_split )

# Create preprocessing recipe
final_rec <- recipe(AUCt ~ ., data = final_train) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_zv(all_numeric_predictors())

# Importance plot
# rsv_ml_train <- rsv_ml_train %>% select(C05, C2, C4, AUCt, PREG, ATV, FED, LOWFAT, EFV, WT)
final_rec2 <- recipe(AUCt ~ ., data = rsv_ml_train)
rsv_ml_rec_prep <-  prep(final_rec2)

final_rf_spec <- rand_forest(mode = "regression",
                             mtry = tune(),
                             trees = tune(),
                             min_n = tune()) %>%
  set_engine("ranger", nthread = 10)
final_rf_spec <- finalize_model(final_rf_spec,results[["best_params_list"]][["size_2"]][["C15_C4rf"]])

final_rf_spec %>% set_engine("ranger", importance = "permutation") %>%
  fit(AUCt ~ ., data = juice(rsv_ml_rec_prep)) %>%
  vip::vip(geom = "col") + theme_bw() +  theme(axis.title.x = element_text(size = 20),
                                               axis.text.x = element_text(size = 17),
                                               axis.text.y = element_text(size = 17),
                                               axis.title.y = element_text(size = 20),
                                               legend.text = element_text(size = 17),
                                               plot.title = element_text(size = 16))

final_wf <- results$workflows$size_2$C15_C4rf

#Save model
final_wf_rf <- results$fits$size_2$C15_C4rf
# saveRDS(final_wf_rf, "model_rfrsv2021.rds")

explainer_external <- explain_tidymodels(
  model = final_wf_rf, 
  data = final_train, 
  y = final_train$AUCt,
  label = "rf")
# saveRDS(explainer_external, file ="explainer_externalrfrsv2021.rds" )

# Predictions on sets
# final_predictions <- final_wf_rf %>% predict(final_train) %>% bind_cols(final_train)

#resample
set.seed(456)
folds <- vfold_cv(final_train, strata = AUCt)

#10 fold CV
rf_rs<- fit_resamples(object = final_wf, 
                       resamples = folds, 
                       control = control_resamples(verbose=T, save_pred = T))


rf_pred_obs <- rf_rs%>% collect_predictions()
CV10F <- calculate_metrics(rf_pred_obs)
CV10F$val <- "CV10F"

Train_pred_obs <- rf_pred_obs %>% 
  left_join(rsv_ml_train, by="AUCt") 

Figure_2a <- Train_pred_obs %>%
  ggplot(mapping = aes(x = .pred, y = AUCt)) + 
  geom_point() +
  geom_smooth(method=lm, color="black") + 
  labs(x="Reference AUC (ng*h/mL)", y= "Predicted rosuvastatin AUC (ng*h/mL)",
       color="", title = "A") + 
  theme_bw()+
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17),
        plot.title = element_text(size = 16))
Figure_2a

final_res <- final_wf %>%
  last_fit(final_split,type = "conf_int")

#Performances biais imprecision test
final_res_predictions <- final_res %>% collect_predictions()
ab <- calculate_metrics(final_res_predictions)
ab$val <- "Test"


Test_pred_obs <- final_res_predictions %>% left_join(rsv_ml_test, by="AUCt")

Figure_2b <- Test_pred_obs %>%
  ggplot(mapping = aes(x = AUCt, y = .pred)) +
  geom_point() +
  geom_smooth(method=lm,color = "black")+
  labs(x = "Reference AUC (ng*h/mL)", y = "Predicted rosuvastatin AUC (ng*h/mL)", 
       color ="", title ="B") +  
  theme_bw() +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17),
        plot.title = element_text(size = 16))
Figure_2b


#Alef external validation
predictions <- predict(final_wf_rf, new_data = df1)
predictions$.pred
df1$AUC_pred <- predictions$.pred
df1$.pred <- predictions$.pred


df1 <- df1 %>% mutate(phase = case_when(
  phase == "TT" ~ "Third trimester",
  phase == "PP" ~ "Postpartum"))
df1$phase <- as.factor(df1$phase)
df1 %>%
  ggplot(mapping = aes(x = AUCt, y = AUC_pred, colour = phase)) +
  geom_point() +
  geom_smooth(method=lm,color = "black")+
  labs(x = "Reference AUC (ng*h/mL)", y = "Predicted rosuvastatin AUC (ng*h/mL)", 
       color="", title = "B") + 
  theme_bw()+
  scale_colour_manual(values = c("black","red")) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17),
        plot.title = element_text(size = 16)) 


#Plots


ae2 <- calculate_metrics(df1)
ae2$val <- "Alef"

#Leandro External validation, metrics and plots (Phase 2)

predictions <- predict(final_wf_rf, new_data = df2)
df2$AUC_pred <- predictions$.pred
df2$.pred <- predictions$.pred

be2 <- calculate_metrics(df2)
be2$val <- "Leandro"

df2$group <- as.factor(df2$group)
df2 <- df2 %>% mutate(group = case_when(
  group == "Group 1" ~ "F1/F2",
  group == "Group 2" ~ "F3/F4"))

#Plots
df2 %>%
  ggplot(mapping = aes(x = AUCt, y = AUC_pred, colour = group)) +
  geom_point() +
  geom_smooth(method=lm,color = "black")+
  labs(x = "Reference AUC (ng*h/mL)", y = "Predicted rosuvastatin AUC (ng*h/mL)", 
       color ="METAVIR score", title ="A") +  
  theme_bw() +
  scale_colour_manual(values = c("black","red")) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17),
        plot.title = element_text(size = 16))

#Final metrics
tab <- rbind(CV10F,ab,ae2,be2)

#Correlation pairs plot
datagraf <- final_data 

# Id numeric and binary variables
bin_vars <- c()  #  bináry
num_vars <- c("AUCt", "C15", "C4")  # numeric

# ggpairs customization
custom_plot <- function(data, mapping, ...) {
  x_name <- as.character(mapping$x)[2]
  y_name <- as.character(mapping$y)[2]
  
  if ((x_name %in% bin_vars & y_name %in% num_vars) || (x_name %in% num_vars & y_name %in% bin_vars)) {
    ggplot(data, mapping) + geom_boxplot()
  } else if (x_name %in% bin_vars & y_name %in% bin_vars) {
    data %>%
      count(!!sym(x_name), !!sym(y_name)) %>%
      ggplot(aes(x = !!sym(x_name), y = n, fill = as.factor(!!sym(y_name)))) +
      geom_col(position = "dodge") +
      scale_fill_brewer(palette = "Set1") +
      theme_minimal() +
      labs(fill = y_name)
  } else {
    ggplot(data, mapping) + geom_point(alpha = 0.5) + geom_smooth(method = "lm", se = FALSE)
  }
}

# Plot
p1 <- ggpairs(datagraf, 
              columns = c(bin_vars, num_vars),
              lower = list(continuous = wrap("points"), combo = wrap(custom_plot), discrete = wrap(custom_plot)),
              upper = list(continuous = wrap("cor")))

print(p1)

####
rsv_ml_train$grupo <- "grupo"
rsv_ml_test$grupo <- "grupo"
df1$grupo <- "grupo"
df2$grupo <- "grupo"


model <- "$PROBLEM ROSUVASTATIN
$INPUT ID TIME TAD AMT RATE DOSE DV CMT MDV EVID STUDY
$DATA ../DATASET/DATASET_ROSUVASTATIN_V01.csv IGNORE=@
$SUBROUTINES ADVAN6 TOL=6
$MODEL
NCOMPARTMENTS = 4
COMP = (DEPOT1)
COMP = (DEPOT2)
COMP = (CENTRAL, DEFOBS)
COMP = (PERIPH1)

$PK
KA1 = THETA(1)
ALAG2 = THETA(2)
CL = THETA(3)*EXP(ETA(2)) ; Elimination from central
V3 = THETA(4) ; Central volume
Q = THETA(5) ; Q
V4 = THETA(6) ; Peripheral volume
S3 = V3/1000
FTOT = THETA(7)
VF2 = THETA(8)

PHI_2 = LOG(VF2/(1- VF2))
VF2_2 = EXP(PHI_2 + ETA(1))/(1 + EXP(PHI_2 + ETA(1)))
F2 = VF2_2*FTOT ; Bioavailability of Depot 2

VF1 = (1- VF2_2)
F1 = VF1*FTOT*EXP(ETA(3)) ; Bioavailability of Depot 1

$DES
K30 = CL/V3
K34 = Q/V3
K43 = Q/V4

DADT(1) =- KA1*A(1)
DADT(2) =- KA1*A(2)
DADT(3) = KA1*A(1) + KA1*A(2)- K30*A(3)- K34*A(3) + K43*A(4)
DADT(4) = K34*A(3)- K43*A(4)

$ERROR
IPRED = F
DEL = 0
IF(IPRED.EQ.0) DEL = 0.0001

W = F
IRES = DV- IPRED
IWRES = IRES/(W + DEL)
Y = IPRED + W*EPS(1) + EPS(2)

$THETA
(0, 0.3) ; KA1
(0, 4) ; ALAG2
(0, 55) ; CL
(0, 70) ; V3
(0, 10) ; Q
(0, 300) ; V4
(0, 0.2, 1) ; FTOT
(0, 0.5) ; VF2

$OMEGA
0.1 ; IIV VF2
0.01 ; IIV CL
0.01 ; IIV FTOT

$SIGMA
0.01 ; Prop RE
0.00006 FIX ; Add RE

$EST METHOD=1 INTER MAXEVAL=9999 NOABORT PRINT=
$COV PRINT=E
$TABLE ID TIME TAD ETAS(1:LAST) DV IPRED IWRES CWRES EVID F1 F2 ALAG2 NOPRINT ONEHEADER FILE=sdtab0007b
"

