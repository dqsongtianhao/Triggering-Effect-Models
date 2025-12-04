setwd("~/paper1/Table1")
library(dplyr)

n_list <- c(100, 400, 1000)

IQR <- function(x) (quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T)) / 1.349

bind_rows(apply(expand.grid(n_list, 1:5), 1, function(x) {load(paste0("n",x[1],"_",x[2],".RData")); get(paste0("n",x[1],"_",x[2]))[1:1000, ]})) %>%
  group_by(n_sub) %>%
  summarise(conv = n() / 5000,
            n_event = mean(n_event0 / n_sub)) %>%
  round(3) %>%
  View()
            
# Method 1
# report both empirical se and IQR-adjusted se to avoid extreme values
temp <- bind_rows(apply(expand.grid(n_list, 1:5), 1, function(x) {load(paste0("n",x[1],"_",x[2],".RData")); get(paste0("n",x[1],"_",x[2]))[1:1000, ]})) %>%
  filter(conv1 == 0, !is.na(gamma1_sd), !is.na(alpha11_sd), !is.na(alpha12_sd), !is.na(alpha13_sd), !is.na(beta11_sd), !is.na(beta12_sd), !is.na(beta13_sd)) %>%
  group_by(n_sub) %>%
  summarise(conv = n() / 5000,
            across(.cols = starts_with(c("gamma1", "alpha1", "beta1", "cH1")) & ends_with(as.character(0:9)), .fns = list(mean = mean), .names = "{.col}_{.fn}"),
            across(.cols = starts_with(c("gamma1", "alpha1", "beta1")) & ends_with(as.character(0:9)), .fns = list(median = median), .names = "{.col}_{.fn}"),
            gamma1_cov = mean(abs(gamma1 - 0.5)/gamma1_sd < qnorm(0.975)),
            alpha11_cov = mean(abs(alpha11 - 0.5)/alpha11_sd < qnorm(0.975)),
            alpha12_cov = mean(abs(alpha12 - 0.5)/alpha12_sd < qnorm(0.975)),
            alpha13_cov = mean(abs(alpha13 - 0.5)/alpha13_sd < qnorm(0.975)),
            beta11_cov = mean(abs(beta11 - 0.5)/beta11_sd < qnorm(0.975)),
            beta12_cov = mean(abs(beta12 - 0.5)/beta12_sd < qnorm(0.975)),
            beta13_cov = mean(abs(beta13 - 0.5)/beta13_sd < qnorm(0.975)),
            across(.cols = starts_with(c("gamma1", "alpha1", "beta1")) & ends_with(as.character(0:9)), .fns = list(sd = sd), .names = "{.col}_emp_se"),
            across(.cols = starts_with(c("gamma1", "alpha1", "beta1")) & ends_with(as.character(0:9)), .fns = list(sd = IQR), .names = "{.col}_IQR_se"),
            across(.cols = starts_with(c("gamma1", "alpha1", "beta1")) & ends_with("sd"), .fns = list(mean = mean), .names = "{.col}_model_se")) %>%
  round(3)

do.call("paste", c((temp %>% ungroup() %>% select(starts_with("gamma1") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha11") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha12") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha13") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta11") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta12") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta13") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("cH1") & !ends_with("se"))), sep = " & ", collapse = " & "))

do.call("paste", c((temp %>% ungroup() %>% select(starts_with("gamma1") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha11") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha12") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha13") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta11") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta12") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta13") & ends_with("se"))), sep = " & ", collapse = " & "))


# Method 2
# report both empirical se and IQR-adjusted se to avoid extreme values
temp <- bind_rows(apply(expand.grid(n_list, 1:5), 1, function(x) {load(paste0("n",x[1],"_",x[2],".RData")); get(paste0("n",x[1],"_",x[2]))[1:1000, ]})) %>%
  filter(conv2 == 0, !is.na(gamma2_sd), !is.na(alpha21_sd), !is.na(alpha22_sd), !is.na(alpha23_sd), !is.na(alpha24_sd), !is.na(beta2_sd)) %>%
  group_by(n_sub) %>%
  summarise(conv = n() / 5000,
            across(.cols = starts_with(c("gamma2", "alpha2", "beta2", "cH2")) & ends_with(as.character(0:9)), .fns = list(mean = mean), .names = "{.col}_{.fn}"),
            across(.cols = starts_with(c("gamma2", "alpha2", "beta2")) & ends_with(as.character(0:9)), .fns = list(median = median), .names = "{.col}_{.fn}"),
            gamma2_cov = mean(abs(gamma2 - 0.5)/gamma2_sd < qnorm(0.975)),
            alpha21_cov = mean(abs(alpha21 - 0.5)/alpha21_sd < qnorm(0.975)),
            alpha22_cov = mean(abs(alpha22 - 0.5)/alpha22_sd < qnorm(0.975)),
            alpha23_cov = mean(abs(alpha23 - 0.5)/alpha23_sd < qnorm(0.975)),
            alpha24_cov = mean(abs(alpha24 - 0)/alpha24_sd < qnorm(0.975)),
            beta2_cov = mean(abs(beta2 - 0.5)/beta2_sd < qnorm(0.975)),
            across(.cols = starts_with(c("gamma2", "alpha2", "beta2")) & ends_with(as.character(0:9)), .fns = list(sd = sd), .names = "{.col}_emp_se"),
            across(.cols = starts_with(c("gamma2", "alpha2", "beta2")) & ends_with(as.character(0:9)), .fns = list(sd = IQR), .names = "{.col}_IQR_se"),
            across(.cols = starts_with(c("gamma2", "alpha2", "beta2")) & ends_with("sd"), .fns = list(mean = mean), .names = "{.col}_model_se")) %>%
  round(3)

do.call("paste", c((temp %>% ungroup() %>% select(starts_with("gamma2") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha21") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha22") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha23") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha24") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta2") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("cH2") & !ends_with("se"))), sep = " & ", collapse = " & "))

do.call("paste", c((temp %>% ungroup() %>% select(starts_with("gamma2") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha21") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha22") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha23") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha24") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta2") & ends_with("se"))), sep = " & ", collapse = " & "))
