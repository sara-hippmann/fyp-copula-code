library(copula)
library(readr)
library(plotly)
library(GGally)
library(MASS)
library(VineCopula)
library(dplyr)
library(scatterplot3d)
library(rafalib)
library(TSP)
library(fGarch)
library(network)
library(forecast)
library(tseries)
library(FinTS)
library(rugarch)

cryptodata <- read_csv("/Users/sarahippmann/Desktop/thesis resources/cryptodata.csv")

#renaming and formatting series
time <- cryptodata$Dates
btc_pct <- diff(log(cryptodata$`BTC Index`))
eth_pct <- diff(log(cryptodata$`ETH Index`))
bbgalaxy_pct <- diff(log(cryptodata$`BGCI Index`))
fstok10_pct <- diff(log(cryptodata$`FSTOK10 Index`))
fstok40_pct <- diff(log(cryptodata$`FSTOK40 Index`))
blockchain_pct <- diff(log(cryptodata$`SOBLCRTA Index`))
hodl5_pct <- diff(log(cryptodata$`HODL5 Index`))

#obtaining standardised residuals
btc_sarima <- auto.arima(btc_pct, seasonal = TRUE)
summary(btc_sarima)
btc_spec <- ugarchspec(mean.model=list(armaOrder=c(4,0), include.mean=TRUE),
                       variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                       distribution.model = "std")
btc_fit <- ugarchfit(btc_spec, data=btc_pct)
btc_residuals <- residuals(btc_fit, standardize=TRUE)
adf.test(btc_residuals)
kpss.test(btc_residuals)
Box.test(btc_residuals, lag=20, type="Ljung-Box") #test for autocorrelation
ArchTest(btc_residuals, lags=20) # test for ARCH effects
btc_kde <- density(btc_residuals, bw="SJ")
plot(btc_kde, 
     main = "BITCOIN Kernel Density Estimation",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(btc_residuals)
dx <- diff(btc_kde$x)
btc_cdf <- cumsum(btc_kde$y[-1]*dx)
btc_cdf <- btc_cdf / max(btc_cdf)
btc_transformed <- approx(btc_kde$x[-1], btc_cdf, xout=btc_residuals, rule=2)$y
head(btc_transformed)
par(mfrow=c(1,2))
hist(btc_residuals, main = "Bitcoin Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(btc_kde, col = "blue", lwd = 2)
legend("topleft",                    # Position of the legend (you can adjust this)
       legend = c("Bitcoin Data", "KDE"),  # Labels
       fill = c("lightblue", "blue"),  # Colors matching the plot elements
       border = "black",             # Border color of the legend
       bty = "n")    
hist(btc_transformed, main="Transformed Bitcoin (Unif)", xlab="Value", col="lightgreen", freq=FALSE)                # Remove border around legend (optional)
print(btc_kde$bw)

eth_sarima <- auto.arima(eth_pct, seasonal = TRUE)
summary(eth_sarima)
eth_spec <- ugarchspec(mean.model=list(armaOrder=c(2,1), include.mean=TRUE),
                       variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                       distribution.model = "std")
eth_fit <- ugarchfit(eth_spec, data=eth_pct)
eth_residuals <- residuals(eth_fit, standardize=TRUE)
adf.test(eth_residuals)
kpss.test(eth_residuals)
Box.test(eth_residuals, lag=20, type="Ljung-Box") #test for autocorrelation
ArchTest(eth_residuals, lags=20) # test for ARCH effects
eth_kde <- density(eth_residuals, bw="SJ")
plot(eth_kde, 
     main = "ETHEREUM Kernel Density Estimation",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(eth_residuals)
dx <- diff(eth_kde$x)
eth_cdf <- cumsum(eth_kde$y[-1]*dx)
eth_cdf <- eth_cdf / max(eth_cdf)
eth_transformed <- approx(eth_kde$x[-1], eth_cdf, xout=eth_residuals, rule=2)$y
head(eth_transformed)
par(mfrow=c(1,2))
hist(eth_residuals, main = "Original Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(eth_kde, col = "blue", lwd = 2)
hist(eth_transformed, main="Transformed ETHEREUM (Unif)", xlab="Value", col="lightgreen", freq=FALSE)
print(eth_kde$bw)

bgci_sarima <- auto.arima(bbgalaxy_pct, seasonal = TRUE)
summary(bgci_sarima)
bgci_spec <- ugarchspec(mean.model=list(armaOrder=c(2,1), include.mean=TRUE),
                       variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                       distribution.model = "std")
bgci_fit <- ugarchfit(bgci_spec, data=bbgalaxy_pct)
bgci_residuals <- residuals(bgci_fit, standardize=TRUE)
adf.test(bgci_residuals)
kpss.test(bgci_residuals)
Box.test(bgci_residuals, lag=20, type="Ljung-Box") #test for autocorrelation
ArchTest(bgci_residuals, lags=20) # test for ARCH effects
bgci_kde <- density(bgci_residuals, bw="SJ")
plot(bgci_kde, 
     main = "BLOOMBERG CRYPTO KDE",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(bgci_residuals)
dx <- diff(bgci_kde$x)
bgci_cdf <- cumsum(bgci_kde$y[-1]*dx)
bgci_cdf <- bgci_cdf / max(bgci_cdf)
bgci_transformed <- approx(bgci_kde$x[-1], bgci_cdf, xout=bgci_residuals, rule=2)$y
head(bgci_transformed)
par(mfrow=c(1,2))
hist(bgci_residuals, main = "Original Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(bgci_kde, col = "blue", lwd = 2)
hist(bgci_transformed, main="Transformed BLOOMBERG CRYPTO (Unif)", xlab="Value", col="lightgreen", freq=FALSE)
print(bgci_kde$bw)

sobl_sarima <- auto.arima(blockchain_pct, seasonal = TRUE)
summary(sobl_sarima)
sobl_spec <- ugarchspec(mean.model=list(armaOrder=c(4,2), include.mean=TRUE),
                           variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                           distribution.model = "std")
sobl_fit <- ugarchfit(sobl_spec, data=blockchain_pct)
sobl_residuals <- residuals(sobl_fit, standardize=TRUE)
adf.test(sobl_residuals)
kpss.test(sobl_residuals)
Box.test(sobl_residuals, lag=20, type="Ljung-Box") #test for autocorrelation
ArchTest(sobl_residuals, lags=20) # test for ARCH effects
sobl_kde <- density(sobl_residuals, bw="SJ")
plot(sobl_kde, 
     main = "BLOCKCHAIN KDE",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(sobl_residuals)
dx <- diff(sobl_kde$x)
sobl_cdf <- cumsum(sobl_kde$y[-1]*dx)
sobl_cdf <- sobl_cdf / max(sobl_cdf)
sobl_transformed <- approx(sobl_kde$x[-1], sobl_cdf, xout=sobl_residuals, rule=2)$y
head(sobl_transformed)
par(mfrow=c(1,2))
hist(sobl_residuals, main = "Original Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(sobl_kde, col = "blue", lwd = 2)
hist(sobl_transformed, main="Transformed BLOCKCHAIN (Unif)", xlab="Value", col="lightgreen", freq=FALSE)
print(sobl_kde$bw)


hodl_sarima <- auto.arima(hodl5_pct, seasonal = TRUE)
summary(hodl_sarima)
hodl_spec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                        variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                        distribution.model = "std")
hodl_fit <- ugarchfit(hodl_spec, data=hodl5_pct)
hodl_residuals <- residuals(hodl_fit, standardize=TRUE)
adf.test(hodl_residuals)
kpss.test(hodl_residuals)
Box.test(hodl_residuals, lag=20, type="Ljung-Box") #test for autocorrelation
ArchTest(hodl_residuals, lags=20) # test for ARCH effects
hodl_kde <- density(hodl_residuals, bw="SJ")
plot(hodl_kde, 
     main = "TOP 5 COINS KDE",
     xlab = "Value", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
rug(hodl_residuals)
dx <- diff(hodl_kde$x)
hodl_cdf <- cumsum(hodl_kde$y[-1]*dx)
hodl_cdf <- hodl_cdf / max(hodl_cdf)
hodl_transformed <- approx(hodl_kde$x[-1], hodl_cdf, xout=hodl_residuals, rule=2)$y
head(hodl_transformed)
par(mfrow=c(1,2))
hist(hodl_residuals, main = "Original Data", xlab = "Value", col = "lightblue", freq = FALSE)
lines(hodl_kde, col = "blue", lwd = 2)
hist(hodl_transformed, main="Transformed TOP5 COINS (Unif)", xlab="Value", col="lightgreen", freq=FALSE)
print(hodl_kde$bw)

#pairs plots and relationship investigations
df = data.frame(btc_residuals, eth_residuals, bgci_residuals, sobl_residuals, hodl_residuals)
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

#Fitting of CVine copula to data
crypto_df = data.frame(btc_transformed, eth_transformed, bgci_transformed, sobl_transformed, hodl_transformed)

#no independence testing CVine
fit.cryptocv.sum=RVineStructureSelect(crypto_df, selectioncrit = "AIC",
                                     indeptest = FALSE, level=0.05, type="CVine")
fit.cryptocv.sum.mle=RVineMLE(crypto_df, fit.cryptocv.sum)

cryptocv_sum_matrix <- fit.cryptocv.sum.mle$RVM
cryptocv_sum_matrix[[5]] <- c("BTC", "ETH", "BGCI", "SOBL", "HODL5")
plot(cryptocv_sum_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(cryptocv_sum_matrix)

#independence testing CVine
fit.cryptocv.ind=RVineStructureSelect(crypto_df, selectioncrit = "AIC",
                                     indeptest = TRUE, level=0.05, type="CVine")
fit.cryptocv.ind.mle=RVineMLE(crypto_df, fit.cryptocv.ind)

cryptocv_ind_matrix <- fit.cryptocv.ind.mle$RVM
cryptocv_ind_matrix[[5]] <- c("BTC", "ETH", "BGCI", "SOBL", "HODL5")
plot(cryptocv_ind_matrix, tree="ALL", edge.labels="family-tau", type=1)
contour(cryptocv_ind_matrix)

