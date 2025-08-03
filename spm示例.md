

```R
# 安装必要包
if (!require("cmdstanr")) install.packages("cmdstanr")
if (!require("loo")) install.packages("loo")
if (!require("tidyverse")) install.packages("tidyverse")
library(cmdstanr)
library(loo)
library(tidyverse)

# 生成模拟数据 ----------------------------------------------------------------
set.seed(123)
n_year <- 30
true_params <- list(
  Schaefer = list(r = 0.3, K = 10000, q = 0.001),
  Fox = list(r = 0.25, K = 8000, q = 0.0012),
  Pella = list(r = 0.35, K = 12000, q = 0.0008, m = 1.5) # m=1时为Schaefer模型
)

generate_data <- function(model_type) {
  B <- numeric(n_year)
  B[1] <- true_params[[model_type]]$K
  Catch <- round(runif(n_year, 1000, 3000))
  
  for(t in 2:n_year) {
    if(model_type == "Schaefer") {
      B[t] <- B[t-1] + true_params$Schaefer$r*B[t-1]*(1 - B[t-1]/true_params$Schaefer$K) - Catch[t-1]
    } else if(model_type == "Fox") {
      B[t] <- B[t-1] + true_params$Fox$r*B[t-1]*(1 - log(B[t-1])/log(true_params$Fox$K)) - Catch[t-1]
    } else {
      m <- true_params$Pella$m
      B[t] <- B[t-1] + true_params$Pella$r/m*B[t-1]*(1 - (B[t-1]/true_params$Pella$K)^m) - Catch[t-1]
    }
    B[t] <- max(B[t], 100)
  }
  data.frame(Year = 1:n_year, Biomass = B, Catch = Catch)
}

data_schaefer <- generate_data("Schaefer")
data_fox <- generate_data("Fox")
data_pella <- generate_data("Pella")

# 定义Stan模型 ----------------------------------------------------------------
stan_code <- '
data {
  int<lower=0> N;
  vector[N] Catch;
  vector[N] Biomass;
}
parameters {
  real<lower=0> r;
  real<lower=0> K;
  real<lower=0> q;
  real<lower=0> sigma_proc;
  vector[N] B_pred;
}
model {
  r ~ normal(0.3, 0.1);
  K ~ normal(10000, 2000);
  q ~ normal(0.001, 0.0005);
  sigma_proc ~ normal(0, 0.1);
  B_pred[1] ~ normal(Biomass[1], sigma_proc);
  for (t in 2:N) {
    B_pred[t] ~ normal(B_pred[t-1] + r * B_pred[t-1] * (1 - B_pred[t-1] / K) - Catch[t-1], sigma_proc);
  }
  Biomass ~ normal(B_pred, sigma_proc);
}
generated quantities {
  real MSY = r * K / 4;
  real B_current = B_pred[N];
  real catch_limit = r * B_current * (1 - B_current / K);
}
'

# 编译Stan模型
stan_file <- write_stan_file(stan_code)
mod <- cmdstan_model(stan_file)

# 准备数据
stan_data <- list(N = n_year, Catch = data_schaefer$Catch, Biomass = data_schaefer$Biomass)

# 进行贝叶斯参数估计
fit <- mod$sample(data = stan_data, seed = 123, chains = 4, parallel_chains = 4, iter_sampling = 2000, iter_warmup = 1000)

# 提取结果
fit_summary <- fit$summary()
print(fit_summary)

# 提取MSY和可捕量
msy <- fit$draws("MSY")
catch_limit <- fit$draws("catch_limit")

cat("最大可持续产量 (MSY):", mean(msy), "\n")
cat("当前生物量下的可捕量:", mean(catch_limit), "\n")

# 交叉验证 ----------------------------------------------------------------
loo_fit <- loo(fit)
print(loo_fit)
```

要使用R语言编写一个拟合三种剩余产量模型（Surplus Production Models）的代码，并使用贝叶斯方法（通过`cmdstan`）进行参数估算，同时使用交叉验证方法评估模型的拟合效果，可以按照以下步骤进行：

### 1. 安装和加载必要的R包

首先，确保你已经安装了`cmdstanr`和`rstan`包，以及其他必要的包。

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install.packages("rstan")
install.packages("loo")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")

library(cmdstanr)
library(rstan)
library(loo)
library(dplyr)
library(tidyr)
library(ggplot2)
```

### 2. 准备数据

假设你有一个包含年份、捕捞量（Catch）和生物量指数（Biomass Index）的数据集。

```r
# 示例数据
data <- data.frame(
  Year = 1980:2020,
  Catch = rlnorm(41, meanlog = log(100), sdlog = 0.2),
  Biomass_Index = rlnorm(41, meanlog = log(500), sdlog = 0.1)
)
```

### 3. 定义Stan模型

我们将定义三个不同的剩余产量模型，每个模型都包含过程误差项。

#### 模型1: Schaefer模型

```stan
// schaefer_model.stan
data {
  int<lower=0> N; // 数据点数量
  vector[N] Catch; // 捕捞量
  vector[N] Biomass_Index; // 生物量指数
}

parameters {
  real<lower=0> r; // 内禀增长率
  real<lower=0> K; // 环境承载力
  real<lower=0> sigma; // 过程误差
  vector[N] B; // 生物量
}

model {
  // 先验
  r ~ normal(0.5, 0.2);
  K ~ normal(1000, 500);
  sigma ~ normal(0, 1);
  
  // 初始生物量
  B[1] ~ normal(K, sigma);
  
  // 过程模型
  for (t in 2:N) {
    B[t] ~ normal(B[t-1] + r * B[t-1] * (1 - B[t-1] / K) - Catch[t-1], sigma);
  }
  
  // 观测模型
  Biomass_Index ~ normal(B, sigma);
}
```

#### 模型2: Fox模型

```stan
// fox_model.stan
data {
  int<lower=0> N; // 数据点数量
  vector[N] Catch; // 捕捞量
  vector[N] Biomass_Index; // 生物量指数
}

parameters {
  real<lower=0> r; // 内禀增长率
  real<lower=0> K; // 环境承载力
  real<lower=0> sigma; // 过程误差
  vector[N] B; // 生物量
}

model {
  // 先验
  r ~ normal(0.5, 0.2);
  K ~ normal(1000, 500);
  sigma ~ normal(0, 1);
  
  // 初始生物量
  B[1] ~ normal(K, sigma);
  
  // 过程模型
  for (t in 2:N) {
    B[t] ~ normal(B[t-1] * exp(r * (1 - log(B[t-1]) / log(K))) - Catch[t-1], sigma);
  }
  
  // 观测模型
  Biomass_Index ~ normal(B, sigma);
}
```

#### 模型3: Pella-Tomlinson模型

```stan
// pella_tomlinson_model.stan
data {
  int<lower=0> N; // 数据点数量
  vector[N] Catch; // 捕捞量
  vector[N] Biomass_Index; // 生物量指数
}

parameters {
  real<lower=0> r; // 内禀增长率
  real<lower=0> K; // 环境承载力
  real<lower=0> m; // 形状参数
  real<lower=0> sigma; // 过程误差
  vector[N] B; // 生物量
}

model {
  // 先验
  r ~ normal(0.5, 0.2);
  K ~ normal(1000, 500);
  m ~ normal(1, 0.5);
  sigma ~ normal(0, 1);
  
  // 初始生物量
  B[1] ~ normal(K, sigma);
  
  // 过程模型
  for (t in 2:N) {
    B[t] ~ normal(B[t-1] + r / m * B[t-1] * (1 - (B[t-1] / K)^m) - Catch[t-1], sigma);
  }
  
  // 观测模型
  Biomass_Index ~ normal(B, sigma);
}
```

### 4. 编译Stan模型

```r
# 编译Stan模型
schaefer_model <- cmdstan_model("schaefer_model.stan")
fox_model <- cmdstan_model("fox_model.stan")
pella_tomlinson_model <- cmdstan_model("pella_tomlinson_model.stan")
```

### 5. 拟合模型

```r
# 准备数据
stan_data <- list(
  N = nrow(data),
  Catch = data$Catch,
  Biomass_Index = data$Biomass_Index
)

# 拟合Schaefer模型
fit_schaefer <- schaefer_model$sample(
  data = stan_data,
  chains = 4,
  iter_sampling = 2000,
  iter_warmup = 1000,
  seed = 123
)

# 拟合Fox模型
fit_fox <- fox_model$sample(
  data = stan_data,
  chains = 4,
  iter_sampling = 2000,
  iter_warmup = 1000,
  seed = 123
)

# 拟合Pella-Tomlinson模型
fit_pella_tomlinson <- pella_tomlinson_model$sample(
  data = stan_data,
  chains = 4,
  iter_sampling = 2000,
  iter_warmup = 1000,
  seed = 123
)
```

### 6. 交叉验证和模型比较

使用`loo`包进行交叉验证和模型比较。

```r
# 提取对数似然
log_lik_schaefer <- fit_schaefer$draws("log_lik")
log_lik_fox <- fit_fox$draws("log_lik")
log_lik_pella_tomlinson <- fit_pella_tomlinson$draws("log_lik")

# 计算LOO
loo_schaefer <- loo(log_lik_schaefer)
loo_fox <- loo(log_lik_fox)
loo_pella_tomlinson <- loo(log_lik_pella_tomlinson)

# 比较模型
compare(loo_schaefer, loo_fox, loo_pella_tomlinson)
```

### 7. 结果可视化

```r
# 提取MSY和可捕量
msy_schaefer <- fit_schaefer$summary("MSY")
msy_fox <- fit_fox$summary("MSY")
msy_pella_tomlinson <- fit_pella_tomlinson$summary("MSY")

# 可视化结果
ggplot() +
  geom_point(aes(x = data$Year, y = data$Biomass_Index), color = "black") +
  geom_line(aes(x = data$Year, y = fit_schaefer$summary("B")$mean), color = "blue") +
  geom_line(aes(x = data$Year, y = fit_fox$summary("B")$mean), color = "red") +
  geom_line(aes(x = data$Year, y = fit_pella_tomlinson$summary("B")$mean), color = "green") +
  labs(title = "Biomass Index vs. Year", x = "Year", y = "Biomass Index")
```

### 8. 结论

通过交叉验证和模型比较，你可以选择最合适的剩余产量模型来估算MSY和可捕量。`compare`函数将帮助你确定哪个模型的拟合效果最好。

### 注意事项

1. **数据质量**：确保你的数据质量高，特别是在生物量指数和捕捞量的测量上。
2. **模型选择**：根据实际情况选择合适的模型，必要时可以调整模型的复杂度。
3. **计算资源**：贝叶斯方法可能需要较多的计算资源，特别是在数据量较大或模型较复杂时。

希望这些步骤能帮助你完成剩余产量模型的拟合和评估。如果你有更多问题，欢迎继续提问！



```r
# 安装并加载必要的包
install.packages("TMB")
library(TMB)

# 定义剩余产量模型的TMB代码
cpp_code <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(C); // 捕捞量
  DATA_VECTOR(B); // 生物量
  PARAMETER(r);   // 内禀增长率
  PARAMETER(K);   // 环境承载量
  PARAMETER(log_sigma); // 误差项的对数标准差

  Type sigma = exp(log_sigma);
  int n = C.size();
  vector<Type> B_pred(n);
  Type nll = 0.0; // 负对数似然

  B_pred(0) = B(0);
  for(int i = 1; i < n; i++) {
    B_pred(i) = B_pred(i-1) + r * B_pred(i-1) * (1 - B_pred(i-1) / K) - C(i-1);
    nll -= dnorm(B(i), B_pred(i), sigma, true);
  }

  return nll;
}
'

# 编译TMB代码
write(cpp_code, file = "schaefer_model.cpp")
compile("schaefer_model.cpp")
dyn.load(dynlib("schaefer_model"))

# 示例数据
C <- c(100, 150, 200, 250, 300, 350, 400, 450, 500, 550) # 捕捞量
B <- c(1000, 900, 800, 700, 600, 500, 400, 300, 200, 100) # 生物量

# 数据和参数的列表
data <- list(C = C, B = B)
parameters <- list(r = 0.5, K = 1000, log_sigma = log(0.1))

# 构建模型
obj <- MakeADFun(data, parameters, DLL = "schaefer_model")

# 优化模型
opt <- nlminb(obj$par, obj$fn, obj$gr)

# 提取参数估计值
r_est <- opt$par["r"]
K_est <- opt$par["K"]

# 计算MSY
MSY <- r_est * K_est / 4

# 计算可捕量
B_current <- B[length(B)]
catch_limit <- r_est * B_current * (1 - B_current / K_est)

# 输出结果
cat("估计的内禀增长率 (r):", r_est, "\n")
cat("估计的环境承载量 (K):", K_est, "\n")
cat("最大可持续产量 (MSY):", MSY, "\n")
cat("当前生物量下的可捕量:", catch_limit, "\n")
```

```r
# 安装必要包
if (!require("TMB")) install.packages("TMB")
if (!require("tidyverse")) install.packages("tidyverse")
library(TMB)
library(tidyverse)

# 生成模拟数据 ----------------------------------------------------------------
set.seed(123)
n_year <- 30
true_params <- list(
  Schaefer = list(r = 0.3, K = 10000, q = 0.001),
  Fox = list(r = 0.25, K = 8000, q = 0.0012),
  Pella = list(r = 0.35, K = 12000, q = 0.0008, m = 1.5) # m=1时为Schaefer模型
)

generate_data <- function(model_type) {
  B <- numeric(n_year)
  B[1] <- true_params[[model_type]]$K
  Catch <- round(runif(n_year, 1000, 3000))
  
  for(t in 2:n_year) {
    if(model_type == "Schaefer") {
      B[t] <- B[t-1] + true_params$Schaefer$r*B[t-1]*(1 - B[t-1]/true_params$Schaefer$K) - Catch[t-1]
    } else if(model_type == "Fox") {
      B[t] <- B[t-1] + true_params$Fox$r*B[t-1]*(1 - log(B[t-1])/log(true_params$Fox$K)) - Catch[t-1]
    } else {
      m <- true_params$Pella$m
      B[t] <- B[t-1] + true_params$Pella$r/m*B[t-1]*(1 - (B[t-1]/true_params$Pella$K)^m) - Catch[t-1]
    }
    B[t] <- max(B[t], 100)
  }
  
  CPUE <- B * true_params[[model_type]]$q * rnorm(n_year, 1, 0.1)
  data.frame(Year = 1:n_year, Catch = c(NA, Catch[1:(n_year-1)]), CPUE = CPUE)
}

# 编译TMB模板 ----------------------------------------------------------------
compile("SPM.cpp") # 需要创建对应的C++文件
dyn.load(dynlib("SPM"))

# 定义三种模型的TMB实现 --------------------------------------------------------
fit_spm <- function(data, model_type) {
  # 准备TMB数据
  tmb_data <- list(
    model_type = switch(model_type,
                       "Schaefer" = 1,
                       "Fox" = 2,
                       "Pella" = 3),
    Catch = data$Catch,
    CPUE = data$CPUE,
    n = nrow(data)
  )
  
  # 设置初始参数
  parameters <- switch(model_type,
                      "Schaefer" = list(
                        log_r = log(0.2),
                        log_K = log(8000),
                        log_q = log(0.001)),
                      "Fox" = list(
                        log_r = log(0.15),
                        log_K = log(7000),
                        log_q = log(0.0015)),
                      "Pella" = list(
                        log_r = log(0.25),
                        log_K = log(10000),
                        log_q = log(0.001),
                        log_m = log(1.2)))
  
  # 创建TMB对象
  obj <- MakeADFun(
    data = tmb_data,
    parameters = parameters,
    DLL = "SPM",
    silent = TRUE
  )
  
  # 参数优化
  opt <- nlminb(obj$par, obj$fn, obj$gr,
               lower = switch(model_type,
                            "Schaefer" = c(-5, 8, -8),
                            "Fox" = c(-5, 8, -8),
                            "Pella" = c(-5, 8, -8, -2)),
               upper = switch(model_type,
                             "Schaefer" = c(2, 15, 2),
                             "Fox" = c(2, 15, 2),
                             "Pella" = c(2, 15, 2, 3)))
  
  # 提取结果
  report <- sdreport(obj)
  
  # 计算管理参数
  MSY <- switch(model_type,
               "Schaefer" = exp(opt$par[1])*exp(opt$par[2])/4,
               "Fox" = exp(opt$par[1])*exp(opt$par[2])/exp(1),
               "Pella" = {
                 r <- exp(opt$par[1])
                 K <- exp(opt$par[2])
                 m <- exp(opt$par[4])
                 r*K/(m^(1/(m-1)))
               })
  
  list(
    params = report$par.fixed,
    MSY = MSY,
    AIC = 2*length(opt$par) + 2*opt$objective,
    pred = obj$report()$CPUE_pred
  )
}

# 运行三种模型 ----------------------------------------------------------------
models <- c("Schaefer", "Fox", "Pella")
results <- list()

for(model in models) {
  data <- generate_data(model)
  results[[model]] <- fit_spm(data, model)
  results[[model]]$data <- data
}

# 模型比较分析 ----------------------------------------------------------------
model_compare <- map_df(models, ~{
  tibble(
    Model = .x,
    AIC = results[[.x]]$AIC,
    RMSE = sqrt(mean((results[[.x]]$data$CPUE - results[[.x]]$pred)^2)),
    R2 = cor(results[[.x]]$data$CPUE, results[[.x]]$pred)^2
  )
})

# 可视化结果 ------------------------------------------------------------------
plot_data <- map_df(models, ~{
  tibble(
    Year = results[[.x]]$data$Year,
    Observed = results[[.x]]$data$CPUE,
    Predicted = results[[.x]]$pred,
    Model = .x
  )
})

ggplot(plot_data, aes(x = Year)) +
  geom_line(aes(y = Observed, color = "观测值")) +
  geom_line(aes(y = Predicted, color = "预测值"), linetype = 2) +
  facet_wrap(~Model, scales = "free_y") +
  labs(title = "不同剩余产量模型拟合效果比较",
       y = "CPUE",
       color = "") +
  theme_bw()

# 输出结果 ------------------------------------------------------------------
cat("模型比较结果:\n")
print(model_compare)

cat("\n管理参数估计:\n")
walk(models, ~{
  cat("\n", .x, "模型:\n")
  cat("MSY:", round(results[[.x]]$MSY, 1), "\n")
  cat("参数估计:\n")
  print(exp(results[[.x]]$params))
})

```

