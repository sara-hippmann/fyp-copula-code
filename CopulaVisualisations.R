install.packages("plot3D")
library(copula)
library(VineCopula)
library(reshape2)
library(ggplot2)
library(plotly)
library(dplyr)
library(MASS)
library(plot3D)

#FRECHET-HOEFFDING BOUND COPULA CONTOUR PLOTS
n <- 100 
u <- seq(0, 1, length.out = n)
v <- seq(0, 1, length.out = n)
grid <- expand.grid(u = u, v = v)
W <- pmin(grid$u, grid$v)             # Upper bound: min(u, v)
M <- pmax(grid$u + grid$v - 1, 0)     # Lower bound: max(u + v - 1, 0)
Pi <- grid$u * grid$v                 # Independence copula: u * v
W_matrix <- matrix(W, nrow = n, ncol = n)
M_matrix <- matrix(M, nrow = n, ncol = n)
Pi_matrix <- matrix(Pi, nrow = n, ncol = n)
par(mfrow = c(1, 3)) 
contour(u, v, W_matrix, main = "Upper Bound (W)", xlab = "u", ylab = "v", col = "#74abe3")
contour(u, v, M_matrix, main = "Lower Bound (M)", xlab = "u", ylab = "v", col = "#74abe3")
contour(u, v, Pi_matrix, main = "Independence Copula (Π)", xlab = "u", ylab = "v", col = "#74abe3")

# FRECHET-HOEFFDING BOUNDS 3D PLOTS
frechet_hoe = function(u, v) {
  pmin(u, v)
}
u_vals <- seq(0, 1, length.out = 50)
v_vals <- seq(0, 1, length.out = 50)
grid <- expand.grid(u = u_vals, v = v_vals)
grid$z <- mapply(frechet_hoe, grid$u, grid$v)
persp(x = u_vals, y = v_vals, z = matrix(grid$z, nrow = length(u_vals), byrow = TRUE),
      xlab = "U", ylab = "V", zlab = "M(u,v)",
      main = "Upper Fréchet-Hoeffding Bound",
      col = "#74abe3", theta = 0, phi = 30)
frechet_hoe_lower = function(u, v) {
  pmax(0, u + v - 1)
}
u_vals <- seq(0, 1, length.out = 50)
v_vals <- seq(0, 1, length.out = 50)
grid <- expand.grid(u = u_vals, v = v_vals)
grid$z <- mapply(frechet_hoe_lower, grid$u, grid$v)
persp(x = u_vals, y = v_vals, z = matrix(grid$z, nrow = length(u_vals), byrow = TRUE),
      xlab = "U", ylab = "V", zlab = "W(u,v)",
      main = "Lower Fréchet-Hoeffding Bound",
      col = "#74abe3", theta = 0, phi = 30)
independence_copula = function(u, v) {
  u * v
}
u_vals <- seq(0, 1, length.out = 50)
v_vals <- seq(0, 1, length.out = 50)
grid <- expand.grid(u = u_vals, v = v_vals)
grid$z <- mapply(independence_copula, grid$u, grid$v)
persp(x = u_vals, y = v_vals, z = matrix(grid$z, nrow = length(u_vals), byrow = TRUE),
      xlab = "U", ylab = "V", zlab = "Pi(u,v)",
      main = "Independence Copula",
      col = "#74abe3", theta = 30, phi = 30)

#GAUSSIAN COPULA
par(mfrow = c(1, 2)) 
norm_cop <- normalCopula(param = 0.5, dim = 2)
contour(norm_cop, main = "Contour Plot of Gaussian Copula", 
        FUN = dCopula, xlab="U", ylab="V")

u <- seq(0, 1, length.out = 50)
v <- seq(0, 1, length.out = 50)
grid <- expand.grid(u = u, v = v)
grid_matrix <- as.matrix(grid)
z <- dCopula(grid_matrix, copula = norm_cop)
z_matrix <- matrix(z, nrow = length(u), ncol = length(v))
persp(u, v, z_matrix, theta = 30, phi = 30, col = "#74abe3",
      xlab = "U", ylab = "V", zlab = "Density", main = "Density Plot of Gaussian Copula")


#STUDENT T COPULA
par(mfrow = c(1, 2)) 
t_cop <- tCopula(param = 0.5, dim = 2, df = 4)
contour(t_cop, main = "Contour Plot of Student t Copula", 
        FUN = dCopula, xlab="U", ylab="V")

u <- seq(0, 1, length.out = 50)
v <- seq(0, 1, length.out = 50)
grid <- expand.grid(u = u, v = v)
grid_matrix <- as.matrix(grid)
z <- dCopula(grid_matrix, copula = t_cop)
z_matrix <- matrix(z, nrow = length(u), ncol = length(v))
persp(u, v, z_matrix, theta = 30, phi = 30, col = "#74abe3",
      xlab = "U", ylab = "V", zlab = "Density", main = "Density Plot of Student t Copula")

#GUMBEL COPULA
par(mfrow = c(1, 2)) 
gumbel_cop <- gumbelCopula(param = 2, dim = 2) 
contour(gumbel_cop, main = "Contour Plot of Gumbel Copula", 
        FUN = dCopula, xlab="U", ylab="V")

u <- seq(0, 1, length.out = 50)
v <- seq(0, 1, length.out = 50)
grid <- expand.grid(u = u, v = v)
grid_matrix <- as.matrix(grid)
z <- dCopula(grid_matrix, copula = gumbel_cop)
z_matrix <- matrix(z, nrow = length(u), ncol = length(v))
persp(u, v, z_matrix, theta = 30, phi = 30, col = "#74abe3",
      xlab = "U", ylab = "V", zlab = "Density", main = "Density Plot of Gumbel Copula")

#JOE COPULA
par(mfrow = c(1, 2)) 
joe_cop <- joeCopula(param = 2, dim = 2)
contour(joe_cop, main = "Contour Plot of Joe Copula", 
        FUN = dCopula, xlab ="U", ylab="V")

u <- seq(0, 1, length.out = 50)
v <- seq(0, 1, length.out = 50)
grid <- expand.grid(u = u, v = v)
grid_matrix <- as.matrix(grid)
z <- dCopula(grid_matrix, copula = joe_cop)
z_matrix <- matrix(z, nrow = length(u), ncol = length(v))
persp(u, v, z_matrix, theta = 30, phi = 30, col = "#74abe3",
      xlab = "U", ylab = "V", zlab = "Density", main = "Density Plot of Joe Copula")

#CLAYTON COPULA
par(mfrow = c(1, 2)) 
clayton_cop <- claytonCopula(param = 2, dim = 2)
contour(clayton_cop, main = "Contour Plot of Clayton Copula", 
        FUN = dCopula, xlab="U", ylab="V")

u <- seq(0, 1, length.out = 50)
v <- seq(0, 1, length.out = 50)
grid <- expand.grid(u = u, v = v)
grid_matrix <- as.matrix(grid)
z <- dCopula(grid_matrix, copula = clayton_cop)
z_matrix <- matrix(z, nrow = length(u), ncol = length(v))
persp(u, v, z_matrix, theta = 30, phi = 30, col = "#74abe3",
      xlab = "U", ylab = "V", zlab = "Density", main = "Density Plot of Clayton Copula")

#FRANK COPULA
par(mfrow=c(1,2))
frank_cop <- frankCopula(param = 2, dim = 2)
contour(frank_cop, main = "Contour Plot of Frank Copula", 
        FUN = dCopula, xlab="U", ylab="V")

u <- seq(0, 1, length.out = 50)
v <- seq(0, 1, length.out = 50)
grid <- expand.grid(u = u, v = v)
grid_matrix <- as.matrix(grid)
z <- dCopula(grid_matrix, copula = frank_cop)
z_matrix <- matrix(z, nrow = length(u), ncol = length(v))
persp(u, v, z_matrix, theta = 30, phi = 30, col = "#74abe3",
      xlab = "U", ylab = "V", zlab = "Density", main = "Density Plot of Frank Copula")
