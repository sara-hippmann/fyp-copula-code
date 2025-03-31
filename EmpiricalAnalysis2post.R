library(plotly)
library(GGally)
library(copula)
library(MASS)
library(VineCopula)
library(dplyr)
library(readr)
library(scatterplot3d)
library(rafalib)
library(TSP)
library(fGarch)
library(network)
library(ggplot2)
library(lmtest)
library(pscl)
library(rugarch)
library(tseries)

data <- read_csv("/Users/sarahippmann/Desktop/thesis resources/thesisdata2post.csv")
data <- data %>%
  filter(data$`SPX Index` != 0, data$`LUATTRUU Index` != 0, data$`BCOMGC Index` != 0, data$`BCOMCO Index` != 0, data$`FNER Index` != 0,
         data$`BGCI Index` != 0)
spx_pct_change <- diff(log(data$`SPX Index`))
bond_pct_change <- diff(log(data$`LUATTRUU Index`))
gold_pct_change <- diff(log(data$`BCOMGC Index`))
oil_pct_change <- diff(log(data$`BCOMCO Index`))
re_pct_change <- diff(log(data$`FNER Index`))
crypto_pct_change <- diff(log(data$`BGCI Index`))

#SPX marginals
spx_spec <- ugarchspec(mean.model=list(armaOrder=c(2,2), include.mean=TRUE),
                       variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                       distribution.model = "std")
spx_fit <- ugarchfit(spx_spec, data=spx_pct_change)
spx_residuals <- residuals(spx_fit, standardize=TRUE)
adf.test(spx_residuals)
spx_kde <- density(spx_residuals, bw="SJ")
plot(spx_kde, 
     main = "Kernel Density Estimation",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(spx_residuals)
dx <- diff(spx_kde$x)
spx_cdf <- cumsum(spx_kde$y[-1]*dx)
spx_cdf <- spx_cdf / max(spx_cdf)
spx_transformed <- approx(spx_kde$x[-1], spx_cdf, xout=spx_residuals, rule=2)$y
head(spx_transformed)
par(mfrow=c(1,2))
hist(spx_residuals, main = "Original Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(spx_kde, col = "blue", lwd = 2)
hist(spx_transformed, main="Transformed SPX (Unif)", xlab="Value", col="lightgreen", freq=FALSE)

#BND marginals
bnd_spec <- ugarchspec(mean.model=list(armaOrder=c(0,2)),
                       variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                       distribution.model = "std")
bnd_fit <- ugarchfit(bnd_spec, data=bond_pct_change)
bnd_residuals <- residuals(bnd_fit, standardize=TRUE)
adf.test(bnd_residuals)
bnd_kde <- density(bnd_residuals, bw="SJ")
plot(bnd_kde, 
     main = "Kernel Density Estimation",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(bnd_residuals)
dx <- diff(bnd_kde$x)
bnd_cdf <- cumsum(bnd_kde$y[-1]*dx)
bnd_cdf <- bnd_cdf / max(bnd_cdf)
bnd_transformed <- approx(bnd_kde$x[-1], bnd_cdf, xout=bnd_residuals, rule=2)$y
head(bnd_transformed)
par(mfrow=c(1,2))
hist(bnd_residuals, main = "Original Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(bnd_kde, col = "blue", lwd = 2)
hist(bnd_transformed, main="Transformed BND (Unif)", xlab="Value", col="lightgreen", freq=FALSE)

#GOLD marginals
gold_spec <- ugarchspec(mean.model=list(armaOrder=c(0,0)),
                        variance.model=list(model="sGARCH", garchOrder=c(2,1)),
                        distribution.model = "std")
gold_fit <- ugarchfit(gold_spec, data=gold_pct_change)
gold_residuals <- residuals(gold_fit, standardize=TRUE)
adf.test(gold_residuals)
gold_kde <- density(gold_residuals, bw="SJ")
plot(gold_kde, 
     main = "Kernel Density Estimation",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(gold_residuals)
dx <- diff(gold_kde$x)
gold_cdf <- cumsum(gold_kde$y[-1]*dx)
gold_cdf <- gold_cdf / max(gold_cdf)
gold_transformed <- approx(gold_kde$x[-1], gold_cdf, xout=gold_residuals, rule=2)$y
head(gold_transformed)
par(mfrow=c(1,2))
hist(gold_residuals, main = "Original Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(gold_kde, col = "blue", lwd = 2)
hist(gold_transformed, main="Transformed GOLD (Unif)", xlab="Value", col="lightgreen", freq=FALSE)
print(gold_kde$bw)

#OIL marginals
oil_spec <- ugarchspec(mean.model=list(armaOrder=c(0,0)),
                       variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                       distribution.model = "std")
oil_fit <- ugarchfit(oil_spec, data=oil_pct_change)
oil_residuals <- residuals(oil_fit, standardize=TRUE)
adf.test(oil_residuals)
oil_kde <- density(oil_residuals, bw="SJ")
plot(oil_kde, 
     main = "Kernel Density Estimation",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(oil_residuals)
dx <- diff(oil_kde$x)
oil_cdf <- cumsum(oil_kde$y[-1]*dx)
oil_cdf <- oil_cdf / max(oil_cdf)
oil_transformed <- approx(oil_kde$x[-1], oil_cdf, xout=oil_residuals, rule=2)$y
head(oil_transformed)
par(mfrow=c(1,2))
hist(oil_residuals, main = "Original Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(oil_kde, col = "blue", lwd = 2)
hist(oil_transformed, main="Transformed OIL (Unif)", xlab="Value", col="lightgreen", freq=FALSE)

#RE marginals
re_spec <- ugarchspec(mean.model=list(armaOrder=c(0,2)),
                      variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                      distribution.model = "std")
re_fit <- ugarchfit(re_spec, data=re_pct_change)
re_residuals <- residuals(re_fit, standardize=TRUE)
adf.test(re_residuals)
re_kde <- density(re_residuals, bw="SJ")
plot(re_kde, 
     main = "Kernel Density Estimation",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(re_residuals)
dx <- diff(re_kde$x)
re_cdf <- cumsum(re_kde$y[-1]*dx)
re_cdf <- re_cdf / max(re_cdf)
re_transformed <- approx(re_kde$x[-1], re_cdf, xout=re_residuals, rule=2)$y
head(re_transformed)
par(mfrow=c(1,2))
hist(re_residuals, main = "Original Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(re_kde, col = "blue", lwd = 2)
hist(re_transformed, main="Transformed RE (Unif)", xlab="Value", col="lightgreen", freq=FALSE)
print(re_kde$bw)

#CRYPTO marginals
crypto_spec <- ugarchspec(mean.model=list(armaOrder=c(0,0,0)),
                          variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                          distribution.model = "std")
crypto_fit <- ugarchfit(crypto_spec, data=crypto_pct_change)
crypto_residuals <- residuals(crypto_fit, standardize=TRUE)
adf.test(crypto_residuals)
crypto_kde <- density(crypto_residuals, bw="SJ")
plot(re_kde, 
     main = "Kernel Density Estimation",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(crypto_residuals)
dx <- diff(crypto_kde$x)
crypto_cdf <- cumsum(crypto_kde$y[-1]*dx)
crypto_cdf <- crypto_cdf / max(crypto_cdf)
crypto_transformed <- approx(crypto_kde$x[-1], crypto_cdf, xout=crypto_residuals, rule=2)$y
head(crypto_transformed)
par(mfrow=c(1,2))
hist(crypto_residuals, main = "Original Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(crypto_kde, col = "blue", lwd = 2)
hist(crypto_transformed, main="Transformed CRYPTO (Unif)", xlab="Value", col="lightgreen", freq=FALSE)
print(crypto_kde$bw)

#pairs plots and relationship investigations
df = data.frame(spx_residuals, bnd_residuals, gold_residuals, oil_residuals, re_residuals, crypto_residuals)
pairs(df)

plot <- ggpairs(
  df,
  upper = list(continuous = wrap("cor", size = 4, color = "red")),  # Correlation coefficients in the upper triangle
  lower = list(continuous = wrap("points", alpha = 0.5)),  # Density contour in the lower triangle
  diag = list(continuous = "density")    # Histogram in the diagonal
)
print(plot)

options(digits=2)
round(cor(df, method="kendall"),digits=2)

filled.contour(kde2d(spx_residuals, bnd_residuals), 
               main = "Filled Contour Plot of Bivariate Data",
               xlab = "SPX", 
               ylab = "BND",
               color.palette = heat.colors)
spx_bnd = data.frame(spx_residuals, bnd_residuals)
ggplot(spx_bnd, aes(x = spx_residuals, y = bnd_residuals)) +
  geom_density_2d(color = "blue") +
  labs(title = "Contour Plot of Bivariate Data", x = "X-axis", y = "Y-axis") +
  theme_minimal()

spx_bnd_kde <- kde2d(spx_residuals, bnd_residuals)
plot_ly(x = spx_bnd_kde$x, y = spx_bnd_kde$y, z = spx_bnd_kde$z, type = "surface") %>%
  layout(title = "SPX BND residual contour plot",
         scene = list(xaxis = list(title = "SPX"),
                      yaxis = list(title = "BND"),
                      zaxis = list(title = "Density")))

spx_gold_kde <- kde2d(spx_residuals, gold_residuals)
plot_ly(x = spx_gold_kde$x, y = spx_gold_kde$y, z = spx_gold_kde$z, type = "surface") %>%
  layout(title = "SPX GOLD residual contour plot",
         scene = list(xaxis = list(title = "SPX"),
                      yaxis = list(title = "GOLD"),
                      zaxis = list(title = "Density")))

spx_oil_kde <- kde2d(spx_residuals, oil_residuals)
plot_ly(x = spx_oil_kde$x, y = spx_oil_kde$y, z = spx_oil_kde$z, type = "surface") %>%
  layout(title = "SPX OIL residual contour plot",
         scene = list(xaxis = list(title = "SPX"),
                      yaxis = list(title = "OIL"),
                      zaxis = list(title = "Density")))

spx_re_kde <- kde2d(spx_residuals, re_residuals)
plot_ly(x = spx_re_kde$x, y = spx_re_kde$y, z = spx_re_kde$z, type = "surface") %>%
  layout(title = "SPX RE residual contour plot",
         scene = list(xaxis = list(title = "SPX"),
                      yaxis = list(title = "RE"),
                      zaxis = list(title = "Density")))

spx_crypto_kde <- kde2d(spx_residuals, crypto_residuals)
plot_ly(x = spx_crypto_kde$x, y = spx_crypto_kde$y, z = spx_crypto_kde$z, type = "surface") %>%
  layout(title = "SPX CRYPTO residual contour plot",
         scene = list(xaxis = list(title = "SPX"),
                      yaxis = list(title = "CRYPTO"),
                      zaxis = list(title = "Density")))

bnd_gold_kde <- kde2d(bnd_residuals, gold_residuals)
plot_ly(x = bnd_gold_kde$x, y = bnd_gold_kde$y, z = bnd_gold_kde$z, type = "surface") %>%
  layout(title = "BND GOLD residual contour plot",
         scene = list(xaxis = list(title = "BND"),
                      yaxis = list(title = "GOLD"),
                      zaxis = list(title = "Density")))

bnd_oil_kde <- kde2d(bnd_residuals, oil_residuals)
plot_ly(x = bnd_oil_kde$x, y = bnd_oil_kde$y, z = bnd_oil_kde$z, type = "surface") %>%
  layout(title = "BND OIL residual contour plot",
         scene = list(xaxis = list(title = "BND"),
                      yaxis = list(title = "OIL"),
                      zaxis = list(title = "Density")))

bnd_re_kde <- kde2d(bnd_residuals, re_residuals)
plot_ly(x = bnd_re_kde$x, y = bnd_re_kde$y, z = bnd_re_kde$z, type = "surface") %>%
  layout(title = "BND RE residual contour plot",
         scene = list(xaxis = list(title = "BND"),
                      yaxis = list(title = "RE"),
                      zaxis = list(title = "Density")))

bnd_crypto_kde <- kde2d(bnd_residuals, crypto_residuals)
plot_ly(x = bnd_crypto_kde$x, y = bnd_crypto_kde$y, z = bnd_crypto_kde$z, type = "surface") %>%
  layout(title = "BND CRYPTO residual contour plot",
         scene = list(xaxis = list(title = "BND"),
                      yaxis = list(title = "CRYPTO"),
                      zaxis = list(title = "Density")))

gold_oil_kde <- kde2d(gold_residuals, oil_residuals)
plot_ly(x = gold_oil_kde$x, y = gold_oil_kde$y, z = gold_oil_kde$z, type = "surface") %>%
  layout(title = "GOLD OIL residual contour plot",
         scene = list(xaxis = list(title = "GOLD"),
                      yaxis = list(title = "OIL"),
                      zaxis = list(title = "Density")))

gold_re_kde <- kde2d(gold_residuals, re_residuals)
plot_ly(x = gold_re_kde$x, y = gold_re_kde$y, z = gold_re_kde$z, type = "surface") %>%
  layout(title = "GOLD RE residual contour plot",
         scene = list(xaxis = list(title = "GOLD"),
                      yaxis = list(title = "RE"),
                      zaxis = list(title = "Density")))

gold_crypto_kde <- kde2d(gold_residuals, crypto_residuals)
plot_ly(x = gold_crypto_kde$x, y = gold_crypto_kde$y, z = gold_crypto_kde$z, type = "surface") %>%
  layout(title = "GOLD CRYPTO residual contour plot",
         scene = list(xaxis = list(title = "GOLD"),
                      yaxis = list(title = "CRYPTO"),
                      zaxis = list(title = "Density")))

oil_re_kde <- kde2d(oil_residuals, re_residuals)
plot_ly(x = oil_re_kde$x, y = oil_re_kde$y, z = oil_re_kde$z, type = "surface") %>%
  layout(title = "OIL RE residual contour plot",
         scene = list(xaxis = list(title = "OIL"),
                      yaxis = list(title = "RE"),
                      zaxis = list(title = "Density")))

oil_crypto_kde <- kde2d(oil_residuals, crypto_residuals)
plot_ly(x = oil_crypto_kde$x, y = oil_crypto_kde$y, z = oil_crypto_kde$z, type = "surface") %>%
  layout(title = "OIL CRYPTO residual contour plot",
         scene = list(xaxis = list(title = "OIL"),
                      yaxis = list(title = "CRYPTO"),
                      zaxis = list(title = "Density")))

re_crypto_kde <- kde2d(re_residuals, crypto_residuals)
plot_ly(x = re_crypto_kde$x, y = re_crypto_kde$y, z = re_crypto_kde$z, type = "surface") %>%
  layout(title = "OIL CRYPTO residual contour plot",
         scene = list(xaxis = list(title = "RE"),
                      yaxis = list(title = "CRYPTO"),
                      zaxis = list(title = "Density")))

#GAUSSIAN MULTIVARIATE COPULA
norm_copula <- normalCopula(dim=6)
norm_fit <- fitCopula(norm_copula, cbind(spx_transformed, bnd_transformed, gold_transformed, oil_transformed, re_transformed, crypto_transformed), method="ml")
samples <- rCopula(500, norm_fit@copula)
pairs(as.data.frame(samples),
      main="Pairwise Scatterplots of Fitted Multivariate Gaussian Copula",
      pch=1, col="blue")

norm_loglik <- logLik(norm_fit)
norm_params <- norm_fit@copula@dimension + 1
norm_obs <- nrow(cbind(spx_transformed, bnd_transformed, gold_transformed, oil_transformed, re_transformed, crypto_transformed))
norm_aic <- -2 * norm_loglik + 2 * norm_params
print(norm_aic)
norm_bic <- -2 * norm_loglik + norm_params * log(norm_obs)
print(norm_bic)


#STUDENT T MULTIVARIATE COPULA
student_t_copula <- tCopula(df=7, dim=6)
student_fit <- fitCopula(student_t_copula, cbind(spx_transformed, bnd_transformed, gold_transformed, oil_transformed, re_transformed, crypto_transformed), method="ml", optim.method="BFGS")
summary(student_fit)

student_loglik <- logLik(student_fit)
student_params <- student_fit@copula@dimension + 1
student_obs <- nrow(cbind(spx_transformed, bnd_transformed, gold_transformed, oil_transformed, re_transformed, crypto_transformed))
student_aic <- -2 * student_loglik + 2 * student_params
print(student_aic)
student_bic <- -2 * student_loglik + student_params * log(student_obs)
print(student_bic)

#function for getting the model evaluation criteria
model_eval <-  function(fit=fit.vine, digits=3){
  df <- sum(abs(fit$par) > 0)+sum(fit$par2 > 0)
  out <- round(c(fit$AIC, fit$BIC, fit$logLik, df), digits)
  names(out) <- c("AIC", "BIC", "logLik", "par")
  out
}

#RVINE COPULA
vine_df = data.frame(spx_transformed, bnd_transformed, gold_transformed, oil_transformed, re_transformed, crypto_transformed)
fit.rv.sum=RVineStructureSelect(vine_df, selectioncrit = "AIC",
                                indeptest = FALSE, level=0.05, type="RVine")
fit.rv.sum.mle=RVineMLE(vine_df, fit.rv.sum)

rv_sum_matrix <- fit.rv.sum.mle$RVM
rv_sum_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(rv_sum_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(rv_sum_matrix)
print(fit.rv.sum$AIC)
print(fit.rv.sum$BIC)
print(fit.rv.sum$logLik)

fit.rv.ind=RVineStructureSelect(vine_df, selectioncrit = "AIC",
                                indeptest = TRUE, level=0.05, type="RVine")
fit.rv.ind.mle=RVineMLE(vine_df, fit.rv.ind)

rv_ind_matrix <- fit.rv.ind.mle$RVM
rv_ind_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(rv_ind_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(rv_ind_matrix)
rv_ind_loglik <- logLik(fit.rv.ind.mle)

#RVINE COPULA (T AND GAUSS)
fit.trv.sum=RVineStructureSelect(vine_df, familyset=c(1,2), selectioncrit = "AIC",
                                 indeptest = FALSE, level=0.05, type="RVine")
fit.trv.sum.mle=RVineMLE(vine_df, fit.trv.sum)

trv_sum_matrix <- fit.trv.sum.mle$RVM
trv_sum_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(trv_sum_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(trv_sum_matrix)
trv_sum_loglik <- logLik(fit.trv.sum.mle)

fit.trv.ind=RVineStructureSelect(vine_df, familyset=c(1,2), selectioncrit = "AIC",
                                 indeptest = TRUE, level=0.05, type="RVine")
fit.trv.ind.mle=RVineMLE(vine_df, fit.trv.ind)

trv_ind_matrix <- fit.trv.ind.mle$RVM
trv_ind_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(trv_ind_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(trv_ind_matrix)
trv_ind_loglik <- logLik(fit.trv.ind.mle)

#RVINE COPULA (GAUSS)
fit.grv.sum=RVineStructureSelect(vine_df, familyset=1, selectioncrit = "AIC",
                                 indeptest = FALSE, level=0.05, type="RVine")
fit.grv.sum.mle=RVineMLE(vine_df, fit.grv.sum)

grv_sum_matrix <- fit.trv.sum.mle$RVM
grv_sum_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(grv_sum_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(grv_sum_matrix)
grv_sum_loglik <- logLik(fit.grv.sum.mle)

fit.grv.ind=RVineStructureSelect(vine_df, familyset=1, selectioncrit = "AIC",
                                 indeptest = TRUE, level=0.05, type="RVine")
fit.grv.ind.mle=RVineMLE(vine_df, fit.grv.ind)

grv_ind_matrix <- fit.grv.ind.mle$RVM
grv_ind_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(grv_ind_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(grv_ind_matrix)
grv_ind_loglik <- logLik(fit.grv.ind.mle)

#CVINE COPULA
fit.cv.sum=RVineStructureSelect(vine_df, selectioncrit = "AIC",
                                indeptest = FALSE, level=0.05, type="CVine")
fit.cv.sum.mle=RVineMLE(vine_df, fit.cv.sum)

cv_sum_matrix <- fit.cv.sum.mle$RVM
cv_sum_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(cv_sum_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(cv_sum_matrix)
cv_sum_loglik <- logLik(fit.cv.sum.mle)

fit.cv.ind=RVineStructureSelect(vine_df, selectioncrit = "AIC",
                                indeptest = TRUE, level=0.05, type="CVine")
fit.cv.ind.mle=RVineMLE(vine_df, fit.cv.ind)

cv_ind_matrix <- fit.cv.ind.mle$RVM
cv_ind_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(cv_ind_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(cv_ind_matrix)
cv_ind_loglik <- logLik(fit.cv.ind.mle)

#CVINE COPULA (T AND GAUSS)
fit.tcv.sum=RVineStructureSelect(vine_df, familyset=c(1,2), selectioncrit = "AIC",
                                 indeptest = FALSE, level=0.05, type="CVine")
fit.tcv.sum.mle=RVineMLE(vine_df, fit.tcv.sum)

tcv_sum_matrix <- fit.tcv.sum.mle$RVM
tcv_sum_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(tcv_sum_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(tcv_sum_matrix)
tcv_sum_loglik <- logLik(fit.tcv.sum.mle)

fit.tcv.ind=RVineStructureSelect(vine_df, familyset=c(1,2), selectioncrit = "AIC",
                                 indeptest = TRUE, level=0.05, type="CVine")
fit.tcv.ind.mle=RVineMLE(vine_df, fit.tcv.ind)

tcv_ind_matrix <- fit.tcv.ind.mle$RVM
tcv_ind_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(tcv_ind_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(tcv_ind_matrix)
tcv_ind_loglik <- logLik(fit.tcv.ind.mle)

#CVINE COPULA (GAUSS)
fit.gcv.sum=RVineStructureSelect(vine_df, familyset=1, selectioncrit = "AIC",
                                 indeptest = FALSE, level=0.05, type="CVine")
fit.gcv.sum.mle=RVineMLE(vine_df, fit.gcv.sum)

gcv_sum_matrix <- fit.gcv.sum.mle$RVM
gcv_sum_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(gcv_sum_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(gcv_sum_matrix)
gcv_sum_loglik <- logLik(fit.gcv.sum.mle)

fit.gcv.ind=RVineStructureSelect(vine_df, familyset=1, selectioncrit = "AIC",
                                 indeptest = TRUE, level=0.05, type="CVine")
fit.gcv.ind.mle=RVineMLE(vine_df, fit.grv.ind)

gcv_ind_matrix <- fit.gcv.ind.mle$RVM
gcv_ind_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(gcv_ind_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(gcv_ind_matrix)
gcv_ind_loglik <- logLik(fit.gcv.ind.mle)

#DVINE COPULA
d = dim(vine_df)[2]
M = 1- abs(TauMatrix(vine_df))
hamilton = insert_dummy(TSP(M),label="cut")
sol = solve_TSP(hamilton,method="repetitive_nn")
order = cut_tour(sol,"cut")
DVM= D2RVine(order,family=rep(0,d*(d-1)/2),par=rep(0,d*(d-1)/2))
fit.dv.sum=RVineCopSelect(data=vine_df,familyset=NA,indeptest=FALSE,
                          level=0.05,Matrix=DVM$Matrix,selectioncrit="AIC")
fit.dv.sum.mle=RVineMLE(vine_df,fit.dv.sum)

dv_sum_matrix <- fit.dv.sum.mle$RVM
dv_sum_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(dv_sum_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(dv_sum_matrix)

fit.dv.ind=RVineCopSelect(data=vine_df,familyset=NA,indeptest=TRUE,
                          level=0.05,Matrix=DVM$Matrix,selectioncrit="AIC")
fit.dv.ind.mle=RVineMLE(vine_df,fit.dv.ind)

dv_ind_matrix <- fit.dv.ind.mle$RVM
dv_ind_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(dv_ind_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(dv_ind_matrix)

#DVINE COPULA (T AND GAUSS)
d = dim(vine_df)[2]
M = 1- abs(TauMatrix(vine_df))
hamilton = insert_dummy(TSP(M),label="cut")
sol = solve_TSP(hamilton,method="repetitive_nn")
order = cut_tour(sol,"cut")
DVM= D2RVine(order,family=rep(0,d*(d-1)/2),par=rep(0,d*(d-1)/2))
fit.tdv.sum=RVineCopSelect(data=vine_df,familyset=c(1,2),indeptest=FALSE,
                           level=0.05,Matrix=DVM$Matrix,selectioncrit="AIC")
fit.tdv.sum.mle=RVineMLE(vine_df,fit.tdv.sum)

tdv_sum_matrix <- fit.tdv.sum.mle$RVM
tdv_sum_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(tdv_sum_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(tdv_sum_matrix)

fit.tdv.ind=RVineCopSelect(data=vine_df,familyset=c(1,2),indeptest=TRUE,
                           level=0.05,Matrix=DVM$Matrix,selectioncrit="AIC")
fit.tdv.ind.mle=RVineMLE(vine_df,fit.tdv.ind)

tdv_ind_matrix <- fit.tdv.ind.mle$RVM
tdv_ind_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(tdv_ind_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(tdv_ind_matrix)

#DVINE COPULA (GAUSS)
d = dim(vine_df)[2]
M = 1- abs(TauMatrix(vine_df))
hamilton = insert_dummy(TSP(M),label="cut")
sol = solve_TSP(hamilton,method="repetitive_nn")
order = cut_tour(sol,"cut")
DVM= D2RVine(order,family=rep(0,d*(d-1)/2),par=rep(0,d*(d-1)/2))
fit.gdv.sum=RVineCopSelect(data=vine_df,familyset=1,indeptest=FALSE,
                           level=0.05,Matrix=DVM$Matrix,selectioncrit="AIC")
fit.gdv.sum.mle=RVineMLE(vine_df,fit.gdv.sum)

gdv_sum_matrix <- fit.gdv.sum.mle$RVM
gdv_sum_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(gdv_sum_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(gdv_sum_matrix)

fit.gdv.ind=RVineCopSelect(data=vine_df,familyset=1,indeptest=TRUE,
                           level=0.05,Matrix=DVM$Matrix,selectioncrit="AIC")
fit.gdv.ind.mle=RVineMLE(vine_df,fit.gdv.ind)

gdv_ind_matrix <- fit.gdv.ind.mle$RVM
gdv_ind_matrix[[5]] <- c("SPX", "BND", "GOLD", "OIL", "RE", "CRYPTO")
plot(gdv_ind_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(gdv_ind_matrix)

#MODEL COMPARISON
out.table <- rbind(
  model_eval(fit=fit.rv.sum),
  model_eval(fit=fit.rv.ind),
  model_eval(fit=fit.trv.sum),
  model_eval(fit=fit.trv.ind),
  model_eval(fit=fit.grv.sum),
  model_eval(fit=fit.grv.ind),
  model_eval(fit=fit.cv.sum),
  model_eval(fit=fit.cv.ind),
  model_eval(fit=fit.tcv.sum),
  model_eval(fit=fit.tcv.ind),
  model_eval(fit=fit.gcv.sum),
  model_eval(fit=fit.gcv.ind),
  model_eval(fit=fit.dv.sum),
  model_eval(fit=fit.dv.ind),
  model_eval(fit=fit.tdv.sum),
  model_eval(fit=fit.tdv.ind),
  model_eval(fit=fit.gdv.sum),
  model_eval(fit=fit.gdv.ind)
)
row.names(out.table) <- c(
  "Rvine-seq-all",
  "RVine-mle-all",
  "Rvine-seq-t",
  "Rvine-mle-t",
  "Rvine-seq-g",
  "Rvine-mle-g",
  "Cvine-seq-all",
  "CVine-mle-all",
  "Cvine-seq-t",
  "Cvine-mle-t",
  "Cvine-seq-g",
  "Cvine-mle-g",
  "Dvine-seq-all",
  "DVine-mle-all",
  "Dvine-seq-t",
  "Dvine-mle-t",
  "Dvine-seq-g",
  "Dvine-mle-g"
)
out.table

#VUONG TESTS - TESTING
RVineVuongTest(vine_df, fit.rv.sum.mle$RVM, fit.dv.sum.mle$RVM)$statistic

#VUONG TESTS - SEQUENTIAL MLE
vuong.p<-c(
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.cv.sum.mle$RVM)$p.value,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.dv.sum.mle$RVM)$p.value,
  RVineVuongTest(vine_df,fit.cv.sum.mle$RVM,fit.dv.sum.mle$RVM)$p.value,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.grv.sum.mle$RVM)$p.value)
vuong.p.aic<-c(
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.cv.sum.mle$RVM)$p.value.Akaike,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.dv.sum.mle$RVM)$p.value.Akaike,
  RVineVuongTest(vine_df,fit.cv.sum.mle$RVM,fit.dv.sum.mle$RVM)$p.value.Akaike,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.grv.sum.mle$RVM)$p.value.Akaike)
vuong.p.bic<-c(
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.cv.sum.mle$RVM)$p.value.Schwarz,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.dv.sum.mle$RVM)$p.value.Schwarz,
  RVineVuongTest(vine_df,fit.cv.sum.mle$RVM,fit.dv.sum.mle$RVM)$p.value.Schwarz,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.grv.sum.mle$RVM)$p.value.Schwarz)
vuong.stat<-c(
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.cv.sum.mle$RVM)$statistic,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.dv.sum.mle$RVM)$statistic,
  RVineVuongTest(vine_df,fit.cv.sum.mle$RVM,fit.dv.sum.mle$RVM)$statistic,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.grv.sum.mle$RVM)$statistic)
vuong.stat.aic<-c(
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.cv.sum.mle$RVM)$statistic.Akaike,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.dv.sum.mle$RVM)$statistic.Akaike,
  RVineVuongTest(vine_df,fit.cv.sum.mle$RVM,fit.dv.sum.mle$RVM)$statistic.Akaike,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.grv.sum.mle$RVM)$statistic.Akaike)
vuong.stat.bic<-c(
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.cv.sum.mle$RVM)$statistic.Schwarz,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.dv.sum.mle$RVM)$statistic.Schwarz,
  RVineVuongTest(vine_df,fit.cv.sum.mle$RVM,fit.dv.sum.mle$RVM)$statistic.Schwarz,
  RVineVuongTest(vine_df,fit.rv.sum.mle$RVM,fit.grv.sum.mle$RVM)$statistic.Schwarz)
vuong.table.rcdv<-round(cbind(vuong.stat,vuong.p,vuong.stat.aic,vuong.p.aic,
                              vuong.stat.bic, vuong.p.bic),digits=2)
rownames(vuong.table.rcdv)<-c("rv all vs dv all","rv all vs cv all",
                              "cv all vs dv all", "rv all vs gv")
colnames(vuong.table.rcdv)<-c("stat","p.value","stat-aic","p-value",
                              "stat-bic","p-value")
vuong.table.rcdv