library(ggplot2)
# rmax <- 30
out.df <- matrix(NA, ncol = 2, nrow = 0)
# a <- 0.01
r <- seq(1.8, 4, by = 0.001) # seq(1.9, 2.975, by = 0.0005)

n <- 2000

for (z in 1:length(r)) {
  xl <- vector()
  xl[1] <- 0.1000
  for (i in 2:n) {
    # xl[i] <- xl[i - 1] * r[z] * exp(-a * xl[i - 1])
    xl[i] <- r[z] * xl[i - 1] * (1 - xl[i - 1])
  }
  uval <- unique(xl[(n - 500):n])
  ### Here is where we can save the output for ggplot
  out.df <- rbind(out.df, cbind(rep(r[z], length(uval)), uval))
}
out.df <- as.data.frame(out.df)
colnames(out.df) <- c("r", "N")
ggplot(out.df, aes(x = r, y = N)) +
  geom_point(size = 0.5)

###############################
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(mathjaxr)
library(latex2exp)

# 生长率值
lambda <- seq(0.25, 4, 0.25)
plot.list <- list()

df.2 <- data.frame()

for (lam in lambda) {
  # 状态变量向量
  N <- numeric(30)
  # N的初始值
  N[1] <- 0.01
  for (t in 2:30) {
    N[t] <- lam * N[t - 1] * (1 - N[t - 1])
  }
  df <- data.frame(
    lambda = rep(lam, 30),
    t = 1:30,
    N = N
  )
  df.2 <- rbind(df.2, df)
}

appender <- function(string) {
  TeX(paste("$\\lambda = $", string))
}

df.2 |>
  transform(lambda = as.character(lambda)) |>
  # mutate(lambda1 = expression(paste("lambda = ", lambda)))
  ggplot(aes(x = t, y = N)) +
  geom_line() +
  facet_wrap(~lambda,
    labeller = as_labeller(appender,
      default = label_parsed
    ),
    scales = "free_y"
  )


lambda <- seq(1, 4, 0.01)

df_2 <- data.frame()

for (lam in lambda) {
  N <- numeric(100)
  N[1] <- 0.01
  for (t in 2:100) {
    N[t] <- lam * N[t - 1] * (1 - N[t - 1])
  }
  df <- data.frame(
    lambda = rep(lam, 50),
    t = 51:100,
    N = N[51:100]
  )
  df_2 <- rbind(df_2, df)
}

ggplot() +
  geom_point(data = df_2, aes(x = lambda, y = N), size = 1) +
  theme_bw()


#################
# https://genchanghsu.github.io/2021_Fall_Introduction_to_Theoretical_Ecology/week-4.html

library(tidyverse)

### (1) Set the parameters
r <- 1.8
K <- 500
N0 <- 10
time <- 100

### (2) Define the discrete logistic growth equation
log_fun <- function(r, N, K) {
  N + r * N * (1 - N / K)
}

### (3) Use for loop to iterate over the time sequence
pop_size <- numeric(time)
pop_size[1] <- N0

for (i in 2:time) {
  pop_size[i] <- log_fun(r = r, N = pop_size[i - 1], K = K)
}

pop_data <- pop_size %>%
  as.data.frame() %>%
  rename(., pop_size = `.`) %>%
  mutate(time = 0:(time - 1)) %>%
  relocate(time)

head(pop_data)

### Population trajectory
ggplot(pop_data, aes(x = time, y = pop_size)) +
  geom_point() +
  geom_line() +
  geom_hline(
    yintercept = K,
    color = "red",
    size = 1.2,
    linetype = "dashed"
  ) +
  geom_text(
    x = time * 1.02,
    y = K + 50,
    label = "italic(K)",
    color = "red",
    size = 6.5, parse = T
  ) +
  labs(
    y = expression(italic(N)),
    title = paste0(
      "Discrete logistic growth", "\n",
      "(r = ", r, ", K = ", K, ", N0 = ", N0, ")"
    )
  ) +
  scale_x_continuous(limits = c(0, time * 1.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max(pop_size) * 1.1), expand = c(0, 0)) +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))


### Cobweb plot/logistic map
cobweb_data <- data.frame(
  Nt = rep(pop_size[-time], each = 2),
  Nt1 = c(0, rep(pop_size[-1], each = 2)[-length(rep(pop_size[-1],
    each = 2
  ))])
)

logistic_map <- data.frame(Nt = seq(0, (r + 1) / r * K, by = 0.1)) %>%
  mutate(Nt1 = Nt + r * Nt * (1 - Nt / K))

ggplot() +
  geom_line(
    data = logistic_map, aes(x = Nt, y = Nt1),
    color = "green", size = 1.2
  ) +
  geom_path(
    data = cobweb_data, aes(x = Nt, y = Nt1),
    color = "blue", size = 0.5
  ) +
  geom_abline(
    slope = 1, intercept = 0,
    color = "red", size = 1
  ) +
  labs(
    x = expression(italic(N[t])),
    y = expression(italic(N[t + 1])),
    title = paste0(
      "Cobweb plot/logistic map", "\n",
      "(r = ", r, ", K = ", K, ", N0 = ", N0, ")"
    )
  ) +
  scale_x_continuous(
    limits = c(0, (r + 1) / r * K * 1.05),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, max(pop_size) * 1.1),
    expand = c(0, 0)
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )
