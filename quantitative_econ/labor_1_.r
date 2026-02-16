#Start Date:14/02/2025
#Author: Jordi Torres
#PS Macro 1.


#1. LOAD/MANIPULATE DATA.
library(data.table)

#Note, update paths manually to run the code:
data <- fread("/Users/jorditorresvallverdu/Library/Mobile Documents/com~apple~CloudDocs/tse/year2/term2/macro_labor/ln.data.1.AllData.txt")
series <- fread("/Users/jorditorresvallverdu/Library/Mobile Documents/com~apple~CloudDocs/tse/year2/term2/macro_labor/ln.series.txt")

#This is for the output, updtate too.
setwd("/Users/jorditorresvallverdu/Documents/GitHub/TSE-2025-2026/quantitative_econ")


library(data.table)

#Keep only monthly observations 
data <- data[period != "M13"]
data <- data[period != "Q01"]
data <- data[period != "Q02"]
data <- data[period != "Q03"]
data <- data[period != "Q04"]
data <- data[period != "A01"]


library(dplyr)
library(tidyr)

keep_ids <- c(
  "LNS12000002","LNS13000002","LNS17100002","LNS17200002",
  "LNS12000001","LNS13000001","LNS17100001","LNS17200001"
)

data_clean <- data %>%
  filter(substr(period,1,1) == "M") %>%
  filter(series_id %in% keep_ids) %>%
  mutate(
    value = ifelse(value == ".", NA, value),
    value = as.numeric(value),
    month = as.numeric(substr(period,2,3)),
    date = as.Date(paste(year, month, "01", sep="-"))
  )


wide_data <- data_clean %>%
  select(date, series_id, value) %>%
  pivot_wider(names_from = series_id, values_from = value)

wide_data <- wide_data %>%
  rename(
    E_w  = LNS12000002,
    U_w  = LNS13000002,
    UE_w = LNS17100002,
    EU_w = LNS17200002,
    E_m  = LNS12000001,
    U_m  = LNS13000001,
    UE_m = LNS17100001,
    EU_m = LNS17200001
  )

#restrict the dataset
  wide_data <- wide_data %>%
  filter(!is.na(E_w) & !is.na(E_m)) %>%
  filter(date >= as.Date("1990-01-01") & date<= as.Date("2020-01-01"))


#2. EXERCISE 2

wide_data <- wide_data %>%
  arrange(date) %>%
  mutate(
    F_w = UE_w / lag(U_w),
    S_w = EU_w / lag(E_w),
    F_m = UE_m / lag(U_m),
    S_m = EU_m / lag(E_m)
  )


wide_data <- wide_data %>%
  mutate(
    # Women
    f_w = -log(1 - F_w),
    s_w = -log(1 - S_w),
    u_w = U_w / (U_w + E_w),

    # Men
    f_m = -log(1 - F_m),
    s_m = -log(1 - S_m),
    u_m = U_m / (U_m + E_m)
  )

tab_clean <- data.frame(
  Statistic = c("Mean F", "SD F",
                "Mean S", "SD S",
                "Mean f", "SD f",
                "Mean s", "SD s"),
  Women = c(mean(wide_data$F_w, na.rm=TRUE),
            sd(wide_data$F_w, na.rm=TRUE),
            mean(wide_data$S_w, na.rm=TRUE),
            sd(wide_data$S_w, na.rm=TRUE),
            mean(wide_data$f_w, na.rm=TRUE),
            sd(wide_data$f_w, na.rm=TRUE),
            mean(wide_data$s_w, na.rm=TRUE),
            sd(wide_data$s_w, na.rm=TRUE)),
  Men = c(mean(wide_data$F_m, na.rm=TRUE),
          sd(wide_data$F_m, na.rm=TRUE),
          mean(wide_data$S_m, na.rm=TRUE),
          sd(wide_data$S_m, na.rm=TRUE),
          mean(wide_data$f_m, na.rm=TRUE),
          sd(wide_data$f_m, na.rm=TRUE),
          mean(wide_data$s_m, na.rm=TRUE),
          sd(wide_data$s_m, na.rm=TRUE))
)

library(stargazer)

stargazer(tab_clean,
          summary = FALSE,
          rownames = FALSE,
          digits = 3,
          title = "Job Finding and Separation Rates (Monthly, 1990–2019)",
          label = "tab:rates",
          type = "latex",
          out = "table_rates.tex"
)



library(ggplot2)

df_f <- wide_data %>%
  select(date, f_w, f_m) %>%
  pivot_longer(-date, names_to="group", values_to="f") %>%
  mutate(group = recode(group, f_w="Women", f_m="Men"))

p_f <-ggplot(df_f, aes(x=date, y=f, color=group)) +
  geom_line() +
  labs(title="Job-finding rate (hazard) over time", x=NULL, y="f_t") +
  theme_minimal()

  ggsave("job_finding_rate.png", plot = p_f, width = 7, height = 5)


df_s <- wide_data %>%
  select(date, s_w, s_m) %>%
  pivot_longer(-date, names_to="group", values_to="s") %>%
  mutate(group = recode(group, s_w="Women", s_m="Men"))

p_s <- ggplot(df_s, aes(x=date, y=s, color=group)) +
  geom_line() +
  labs(title="Separation rate (hazard) over time", x=NULL, y="s_t") +
  theme_minimal()


  ggsave("job_separation_rate.png", plot = p_s, width = 7, height = 5)


#3. EXERCISE 3-> SHIMER decomposition

fbar_w <- mean(wide_data$f_w, na.rm=TRUE)
sbar_w <- mean(wide_data$s_w, na.rm=TRUE)

fbar_m <- mean(wide_data$f_m, na.rm=TRUE)
sbar_m <- mean(wide_data$s_m, na.rm=TRUE)

wide_data <- wide_data %>%
  mutate(
    u_w_f = sbar_w / (sbar_w + f_w),
    u_w_s = s_w / (s_w + fbar_w),

    u_m_f = sbar_m / (sbar_m + f_m),
    u_m_s = s_m / (s_m + fbar_m)
  )

  wide_data <- wide_data %>%
  mutate(
    u_w_ss = s_w / (s_w + f_w),
    u_m_ss = s_m / (s_m + f_m)
  )


summary(wide_data[,c("u_w","u_w_ss")])



df_w <- wide_data %>%
  select(date, u_w_ss, u_w_f, u_w_s) %>%
  pivot_longer(-date, names_to="series", values_to="value")

df_g <- ggplot(df_w, aes(x=date, y=value, color=series)) +
  geom_line() +
  labs(title="Shimer Decomposition – Women (Steady-State)",
       y="Unemployment Rate",
       x=NULL) +
  theme_minimal()

  ggsave("shimer_w.png", plot = df_g, width = 7, height = 5)


df_m <- wide_data %>%
  select(date, u_m_ss, u_m_f, u_m_s) %>%
  pivot_longer(-date, names_to="series", values_to="value")

df_m_plot <- ggplot(df_m, aes(x=date, y=value, color=series)) +
  geom_line() +
  labs(title="Shimer Decomposition – Men (Steady-State)",
       y="Unemployment Rate",
       x=NULL) +
  theme_minimal()


  ggsave("shimer_m.png", plot = df_m_plot, width = 7, height = 5)


tab_decomp <- tibble(
  Group = c("Women","Men"),
  SD_u_ss = c(sd(wide_data$u_w_ss, na.rm=TRUE),
              sd(wide_data$u_m_ss, na.rm=TRUE)),
  SD_u_f  = c(sd(wide_data$u_w_f, na.rm=TRUE),
              sd(wide_data$u_m_f, na.rm=TRUE)),
  SD_u_s  = c(sd(wide_data$u_w_s, na.rm=TRUE),
              sd(wide_data$u_m_s, na.rm=TRUE))
) %>%
  mutate(across(-Group, ~ round(.x, 3)))


stargazer(tab_decomp,
          summary = FALSE,
          rownames = FALSE,
          digits = 3,
          title = "Shimer Decomposition: Volatility of Steady-State Unemployment",
          label = "tab:decomp",
          type = "latex",
          out = "table_decomposition.tex")


#4. EXERCISE 4-> Filtering.
##HP filtering
#install.packages("mFilter")

library(mFilter)

#hp
wide_data_hp <- wide_data[
  complete.cases(wide_data[, c("f_w","s_w","u_w_ss","u_w",
                               "f_m","s_m","u_m_ss","u_m")]), ]

hp_f_w <- hpfilter(log(wide_data_hp$f_w), freq = 129600)
hp_s_w <- hpfilter(log(wide_data_hp$s_w), freq = 129600)
hp_u_w <- hpfilter(log(wide_data_hp$u_w_ss), freq = 129600)
hp_u_w_actual <- hpfilter(log(wide_data_hp$u_w), freq = 129600)

hp_f_m <- hpfilter(log(wide_data_hp$f_m), freq = 129600)
hp_s_m <- hpfilter(log(wide_data_hp$s_m), freq = 129600)
hp_u_m <- hpfilter(log(wide_data_hp$u_m_ss), freq = 129600)
hp_u_m_actual <- hpfilter(log(wide_data_hp$u_m), freq = 129600)

cyc_hp <- data.frame(
  c_f_w = hp_f_w$cycle,
  c_s_w = hp_s_w$cycle,
  c_u_w = hp_u_w$cycle,
  c_u_w_actual = hp_u_w_actual$cycle,
  c_f_m = hp_f_m$cycle,
  c_s_m = hp_s_m$cycle,
  c_u_m = hp_u_m$cycle,
  c_u_m_actual = hp_u_m_actual$cycle
)
cyc_hp <- na.omit(cyc_hp)

#linear

t <- 1:nrow(wide_data_hp)

cyc_lin <- data.frame(
  c_f_w = resid(lm(log(wide_data_hp$f_w) ~ t)),
  c_s_w = resid(lm(log(wide_data_hp$s_w) ~ t)),
  c_u_w = resid(lm(log(wide_data_hp$u_w_ss) ~ t)),
  c_u_w_actual = resid(lm(log(wide_data_hp$u_w) ~ t)),
  c_f_m = resid(lm(log(wide_data_hp$f_m) ~ t)),
  c_s_m = resid(lm(log(wide_data_hp$s_m) ~ t)),
  c_u_m = resid(lm(log(wide_data_hp$u_m_ss) ~ t)),
  c_u_m_actual = resid(lm(log(wide_data_hp$u_m) ~ t))
)
cyc_lin <- na.omit(cyc_lin)


#Tabs
tab_sd <- data.frame(
  Method = rep(c("HP","Linear"), each=2),
  Group  = rep(c("Women","Men"), 2),
  SD_f = c(sd(cyc_hp$c_f_w), sd(cyc_hp$c_f_m),
           sd(cyc_lin$c_f_w), sd(cyc_lin$c_f_m)),
  SD_s = c(sd(cyc_hp$c_s_w), sd(cyc_hp$c_s_m),
           sd(cyc_lin$c_s_w), sd(cyc_lin$c_s_m)),
  SD_u = c(sd(cyc_hp$c_u_w), sd(cyc_hp$c_u_m),
           sd(cyc_lin$c_u_w), sd(cyc_lin$c_u_m)),
  SD_u_actual = c(sd(cyc_hp$c_u_w_actual), sd(cyc_hp$c_u_m_actual),
                  sd(cyc_lin$c_u_w_actual), sd(cyc_lin$c_u_m_actual))
)

tab_cor <- data.frame(
  Method = rep(c("HP","Linear"), each=2),
  Group  = rep(c("Women","Men"), 2),
  Corr_f_u = c(cor(cyc_hp$c_f_w, cyc_hp$c_u_w),
               cor(cyc_hp$c_f_m, cyc_hp$c_u_m),
               cor(cyc_lin$c_f_w, cyc_lin$c_u_w),
               cor(cyc_lin$c_f_m, cyc_lin$c_u_m)),
  Corr_s_u = c(cor(cyc_hp$c_s_w, cyc_hp$c_u_w),
               cor(cyc_hp$c_s_m, cyc_hp$c_u_m),
               cor(cyc_lin$c_s_w, cyc_lin$c_u_w),
               cor(cyc_lin$c_s_m, cyc_lin$c_u_m)),
  Corr_f_u_actual = c(cor(cyc_hp$c_f_w, cyc_hp$c_u_w_actual),
                      cor(cyc_hp$c_f_m, cyc_hp$c_u_m_actual),
                      cor(cyc_lin$c_f_w, cyc_lin$c_u_w_actual),
                      cor(cyc_lin$c_f_m, cyc_lin$c_u_m_actual)),
  Corr_s_u_actual = c(cor(cyc_hp$c_s_w, cyc_hp$c_u_w_actual),
                      cor(cyc_hp$c_s_m, cyc_hp$c_u_m_actual),
                      cor(cyc_lin$c_s_w, cyc_lin$c_u_w_actual),
                      cor(cyc_lin$c_s_m, cyc_lin$c_u_m_actual))
)

#round numeric columns 
tab_sd[, c("SD_f","SD_s","SD_u","SD_u_actual")] <- 
  round(tab_sd[, c("SD_f","SD_s","SD_u","SD_u_actual")], 3)

tab_cor[, c("Corr_f_u","Corr_s_u",
            "Corr_f_u_actual","Corr_s_u_actual")] <- 
  round(tab_cor[, c("Corr_f_u","Corr_s_u",
                    "Corr_f_u_actual","Corr_s_u_actual")], 3)

#export
stargazer(tab_sd, summary=FALSE, rownames=FALSE, digits=3,
          title="Cyclical Volatility (SD of Detrended Log Series)",
          label="tab:sd_cycles",
          type="latex", out="table_q4_sd.tex")

stargazer(tab_cor, summary=FALSE, rownames=FALSE, digits=3,
          title="Cyclical Co-movement (Correlations with Unemployment Cycle)",
          label="tab:cor_cycles",
          type="latex", out="table_q4_cor.tex")

##########EOF