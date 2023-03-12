# seroreversion rate
reversion = 32/109
revresion_test = prop.test(32, 109, correct = FALSE)
reversion_ci = c(revresion_test$conf.int[1], revresion_test$conf.int[2]) # 0.2162870 0.3849234

# sensitivity and specificity
se_ci = c(0.837, 0.921)
se = 0.879

sp_ci = c(0.961, 1)
sp = 0.9975

#################
# Montreal-Laval 
#################

# 95 est and CI
ml = prop.test(215, 1578, correct = FALSE)
ml$estimate # 0.1362484 
c(ml$conf.int[1], ml$conf.int[2]) # 0.1202033 0.1540603

# correct for seroreversion
p1_ml = 90/3061
p1_ml_test = prop.test(90, 3061, correct = FALSE)
p1_ml_ci = c(p1_ml_test$conf.int[1], p1_ml_test$conf.int[2])
ml$estimate + p1_ml * reversion # 0.1448802 
c(ml$conf.int[1] + p1_ml_ci[1] * reversion_ci[1], ml$conf.int[2] + p1_ml_ci[2] * reversion_ci[2]) # 0.1253903 0.1679182

# correct for both seroreversion and test-kit performance
(ml$estimate - 1 + sp) / (se - 1 + sp) + reversion * p1_ml # 0.1612256
c((ml$estimate - 1 + sp_ci[1]) / (se_ci[2] - 1 + sp_ci[1]) + p1_ml_ci[1] * reversion_ci[1], 
  (ml$estimate - 1 + sp_ci[2]) / (se_ci[1] - 1 + sp_ci[2]) + p1_ml_ci[2] * reversion_ci[2]) # 0.1154460 0.1766398

#################
# Surrounding 
#################

# 95 est and CI
sr = prop.test(128, 1422, correct = FALSE)
sr$estimate # 0.09001406 
c(sr$conf.int[1], sr$conf.int[2]) # 0.07622221 0.10601507

# correct for seroreversion
p1_sr = 48/1925
p1_sr_test = prop.test(48, 1925, correct = FALSE)
p1_sr_ci = c(p1_sr_test$conf.int[1], p1_sr_test$conf.int[2])
sr$estimate + p1_sr * reversion # 0.09733445 
c(sr$conf.int[1] + p1_sr_ci[1] * reversion_ci[1], sr$conf.int[2] + p1_sr_ci[2] * reversion_ci[2]) # 0.08030107 0.11868052

# correct for both seroreversion and test-kit performance
(sr$estimate - 1 + sp) / (se - 1 + sp) + reversion * p1_sr # 0.1071653 
c((sr$estimate - 1 + sp_ci[1]) / (se_ci[2] - 1 + sp_ci[1]) + p1_sr_ci[1] * reversion_ci[1], 
  (sr$estimate - 1 + sp_ci[2]) / (se_ci[1] - 1 + sp_ci[2]) + p1_sr_ci[2] * reversion_ci[2]) # 0.06191794 0.12020914 

#################
# Other
#################

# 95 est and CI
ot = prop.test(372, 4304, correct = FALSE)
ot$estimate # 0.08643123
c(ot$conf.int[1], ot$conf.int[2]) # 0.07840072 0.09519932

# correct for seroreversion
p1_ot = 35/2705
p1_ot_test = prop.test(35, 2705, correct = FALSE)
p1_ot_ci = c(p1_ot_test$conf.int[1], p1_ot_test$conf.int[2])
ot$estimate + p1_ot * reversion # 0.09022983
c(ot$conf.int[1] + p1_ot_ci[1] * reversion_ci[1], ot$conf.int[2] + p1_ot_ci[2] * reversion_ci[2]) # 0.08041614 0.10210530

# correct for both seroreversion and test-kit performance
(ot$estimate - 1 + sp) / (se - 1 + sp) + reversion * p1_ot # 0.09955585 
c((ot$estimate - 1 + sp_ci[1]) / (se_ci[2] - 1 + sp_ci[1]) + p1_ot_ci[1] * reversion_ci[1], 
  (ot$estimate - 1 + sp_ci[2]) / (se_ci[1] - 1 + sp_ci[2]) + p1_ot_ci[2] * reversion_ci[2]) # 0.05579232 0.11016909

#################
# Total
#################

# 95 est and CI
tot = prop.test(715, 7304, correct = FALSE)
tot$estimate # 0.09789157
c(tot$conf.int[1], tot$conf.int[2]) # 0.09128639 0.10491949

# correct for seroreversion
p1_tot = 173/7691
p1_tot_test = prop.test(173, 7691, correct = FALSE)
p1_tot_ci = c(p1_tot_test$conf.int[1], p1_tot_test$conf.int[2])
tot$estimate + p1_tot * reversion # 0.1044953
c(tot$conf.int[1] + p1_tot_ci[1] * reversion_ci[1], tot$conf.int[2] + p1_tot_ci[2] * reversion_ci[2]) # 0.09548463 0.11494825

# correct for both seroreversion and test-kit performance
(tot$estimate - 1 + sp) / (se - 1 + sp) + reversion * p1_tot # 0.1154361 
c((tot$estimate - 1 + sp_ci[1]) / (se_ci[2] - 1 + sp_ci[1]) + p1_tot_ci[1] * reversion_ci[1], 
  (tot$estimate - 1 + sp_ci[2]) / (se_ci[1] - 1 + sp_ci[2]) + p1_tot_ci[2] * reversion_ci[2]) # 0.07096872 0.12698404