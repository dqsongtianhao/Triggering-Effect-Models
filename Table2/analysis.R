setwd("~/paper1/Table2")
library(dplyr)

n_list <- c(100, 400, 1000)

IQR <- function(x) (quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T)) / 1.349

bind_rows(apply(expand.grid(n_list, 1:5), 1, function(x) {load(paste0("n",x[1],"_",x[2],".RData")); get(paste0("n",x[1],"_",x[2]))[1:1000, ]})) %>%
  group_by(n_sub) %>%
  summarise(conv = n() / 5000,
            n_event = mean(n_event0 / n_sub)) %>%
  round(3) %>%
  View()

temp <- bind_rows(apply(expand.grid(n_list, 1:5), 1, function(x) {load(paste0("n",x[1],"_",x[2],".RData")); get(paste0("n",x[1],"_",x[2]))[1:1000, ]})) %>%
  filter(conv == 0) %>%
  group_by(n_sub) %>%
  summarise(conv = n() / 5000,
            across(.cols = starts_with(c("cH")), .fns = list(mean = mean), .names = "{.col}_{.fn}"),
            across(.cols = starts_with(c("gamma", "alpha", "beta")) & !ends_with("sd"), .fns = list(mean = mean, median = median), .names = "{.col}_{.fn}"),
            gamma_cov = mean(abs(gamma - 0.5)/gamma_sd < qnorm(0.975)),
            alpha1_cov = mean(abs(alpha1 - 0.5)/alpha1_sd < qnorm(0.975)),
            alpha2_cov = mean(abs(alpha2 - 0.5)/alpha2_sd < qnorm(0.975)),
            alpha3_cov = mean(abs(alpha3 - 0.5)/alpha3_sd < qnorm(0.975)),
            beta1_cov = mean(abs(beta1 - 0.5)/beta1_sd < qnorm(0.975)),
            beta2_cov = mean(abs(beta2 - 0.5)/beta2_sd < qnorm(0.975)),
            beta3_cov = mean(abs(beta3 - 0.5)/beta3_sd < qnorm(0.975)),
            across(.cols = starts_with(c("alpha", "beta")) & ends_with(as.character(0:9)) | ends_with("gamma"), .fns = list(sd = sd), .names = "{.col}_emp_se"),
            across(.cols = starts_with(c("alpha", "beta")) & ends_with(as.character(0:9)) | ends_with("gamma"), .fns = list(sd = IQR), .names = "{.col}_IQR_se"),
            across(.cols = starts_with(c("gamma", "alpha", "beta")) & ends_with("sd"), .fns = list(mean = mean), .names = "{.col}_model_se")) %>%
  round(3)

do.call("paste", c((temp %>% ungroup() %>% select(starts_with("gamma") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha1") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha2") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha3") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta1") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta2") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta3") & !ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("cH") & !ends_with("se"))), sep = " & ", collapse = " & "))

do.call("paste", c((temp %>% ungroup() %>% select(starts_with("gamma") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha1") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha2") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("alpha3") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta1") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta2") & ends_with("se"))), sep = " & ", collapse = " & "))
do.call("paste", c((temp %>% ungroup() %>% select(starts_with("beta3") & ends_with("se"))), sep = " & ", collapse = " & "))