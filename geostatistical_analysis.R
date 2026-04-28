# =============================================================================
# Geostatistical Analysis & Spatial Statistics in R
# Author: Nafisa Ahmad
# Institution: La Trobe University — Master of Data Science (AI Specialization)
# Year: 2025
# =============================================================================
# Analyses covered:
#   1. Variogram modelling & ordinary kriging (head dataset)
#   2. Model validation — LOOCV & 5-fold cross-validation
#   3. Monitoring network optimisation
#   4. Spatial point process analysis (hyytiala dataset)
#   5. Spatio-temporal Poisson process simulation
# =============================================================================

# --- Install packages if needed ----------------------------------------------
# install.packages(c("geoR", "gstat", "spatstat", "sp", "lattice"))

# --- Load libraries ----------------------------------------------------------
library(spatstat)
library(gstat)
library(geoR)
library(sp)
library(lattice)


# =============================================================================
# PART 1: KRIGING & VARIOGRAM MODELLING
# =============================================================================

# Load data
data(head)
summary(head)

# --- Exploratory spatial plot ------------------------------------------------
plot(head$coords,
     main = "Spatial Locations of Observation Points",
     xlab = "X-coordinate",
     ylab = "Y-coordinate",
     pch  = 19,
     col  = "black")

# --- Sample variogram --------------------------------------------------------
sample_vario <- variog(head)
plot(sample_vario, main = "Sample Variogram")

# --- Fit Gaussian and Exponential models -------------------------------------
gaussian_model <- variofit(sample_vario,
                           ini.cov.pars = c(200000, 1),
                           cov.model    = "gaussian",
                           nugget       = 1)

exp_model <- variofit(sample_vario,
                      ini.cov.pars = c(200000, 1),
                      cov.model    = "exponential",
                      nugget       = 1)

# Compare models visually
plot(sample_vario, main = "Sample Variogram with Fitted Models")
lines(gaussian_model, col = "blue4", lwd = 2)
lines(exp_model,      col = "red4",  lwd = 2)
legend("bottomright",
       c("Gaussian", "Exponential"),
       col = c("blue4", "red4"),
       lwd = 2)

# AIC comparison
cat("Gaussian model AIC:", gaussian_model$value, "\n")
cat("Exponential model AIC:", exp_model$value, "\n")
# Gaussian model selected — substantially lower AIC

# --- Directional variograms --------------------------------------------------
dir_vario_15 <- variog(head, direction = pi/12, tolerance = pi/6)
dir_vario_30 <- variog(head, direction = pi/6,  tolerance = pi/6)

plot(dir_vario_15, type = "l", col = "blue", lwd = 2,
     main = "Directional Variograms: 15° vs 30°")
lines(dir_vario_30, col = "red", lwd = 2)
legend("topleft",
       c("15° direction", "30° direction"),
       col = c("blue", "red"),
       lwd = 2)

multi_vario <- variog4(head, direction = c(0, pi/4, pi/2, 3*pi/4))
plot(multi_vario, main = "Multidirectional Variograms (0°, 45°, 90°, 135°)")


# =============================================================================
# PART 2: ORDINARY KRIGING & SPATIAL PREDICTION
# =============================================================================

chosen_model <- exp_model

# --- Predict at spatial median -----------------------------------------------
coords   <- head$coords
median_x <- median(coords[, 1])
median_y <- median(coords[, 2])

cat("Median x-coordinate:", round(median_x, 2), "\n")
cat("Median y-coordinate:", round(median_y, 2), "\n")

plot(head$coords,
     main = "Observation Locations with Median Centre",
     xlab = "X-coordinate", ylab = "Y-coordinate",
     pch = 19, col = "black")
points(median_x, median_y, pch = 19, col = "red", cex = 1.5)
legend("topleft",
       legend = c("Observation Points", "Median Centre"),
       pch    = c(19, 19),
       col    = c("black", "red"),
       pt.cex = c(1, 1.5))

pred_points      <- data.frame(x = median_x, y = median_y)
krig_predictions <- krige.conv(head,
                               locations = pred_points,
                               krige     = krige.control(obj.model = exp_model))

cat("Predicted value at median:", round(krig_predictions$predict, 4), "\n")
cat("Prediction variance:",       round(krig_predictions$krige.var, 4), "\n")

# --- Grid prediction ---------------------------------------------------------
grid_x    <- seq(min(head$coords[,1]), max(head$coords[,1]), length.out = 50)
grid_y    <- seq(min(head$coords[,2]), max(head$coords[,2]), length.out = 50)
grid_locs <- as.matrix(expand.grid(x = grid_x, y = grid_y))

krig_grid <- krige.conv(head,
                        locations = grid_locs,
                        krige     = krige.control(obj.model = gaussian_model))

pred_matrix <- matrix(krig_grid$predict,
                      nrow  = length(grid_x),
                      ncol  = length(grid_y),
                      byrow = TRUE)

image(x = grid_x, y = grid_y, z = pred_matrix,
      main = "Kriging Prediction Surface",
      xlab = "X-coordinate", ylab = "Y-coordinate")
contour(x = grid_x, y = grid_y, z = pred_matrix, add = TRUE)
points(head$coords, pch = 19, cex = 0.5)


# =============================================================================
# PART 3: MODEL VALIDATION
# =============================================================================

# --- LOOCV -------------------------------------------------------------------
sp_data <- data.frame(x = head$coords[,1],
                      y = head$coords[,2],
                      z = head$data)
coordinates(sp_data) <- ~x+y

vgm_model <- vgm(psill  = exp_model$cov.pars[1],
                 model  = "Exp",
                 range  = exp_model$cov.pars[2],
                 nugget = exp_model$nugget)

cv_result <- krige.cv(z ~ 1, sp_data, vgm_model)
spplot(cv_result, "residual",
       col.regions = bpy.colors(),
       main        = "LOOCV Residuals")

# --- 5-Fold Cross-Validation -------------------------------------------------
set.seed(123)
n      <- length(head$data)
folds  <- sample(rep(1:5, length.out = n))
cv_5fold <- numeric(n)

for (i in 1:5) {
  test_idx <- which(folds == i)

  train_data <- list(coords = head$coords[-test_idx, ],
                     data   = head$data[-test_idx])
  class(train_data) <- "geodata"

  train_vario <- variog(train_data)
  train_model <- variofit(train_vario,
                          ini.cov.pars = exp_model$cov.pars,
                          cov.model    = "exponential",
                          nugget       = exp_model$nugget)

  test_data <- list(coords = head$coords[test_idx, , drop = FALSE],
                    data   = head$data[test_idx])

  test_pred <- krige.conv(geodata   = train_data,
                          locations = test_data$coords,
                          krige     = krige.control(obj.model = train_model))

  cv_5fold[test_idx] <- test_data$data - test_pred$predict
}

rmse_5fold <- sqrt(mean(cv_5fold^2))
mae_5fold  <- mean(abs(cv_5fold))
bias_5fold <- mean(cv_5fold)

cat("5-Fold CV Results:\n")
cat("RMSE:", round(rmse_5fold, 4), "\n")
cat("MAE: ", round(mae_5fold,  4), "\n")
cat("Bias:", round(bias_5fold, 4), "\n")

plot(head$coords,
     pch  = 19,
     col  = ifelse(abs(cv_5fold) > quantile(abs(cv_5fold), 0.9), "red", "black"),
     main = "5-Fold CV Residuals (Top 10% Errors in Red)",
     xlab = "X-coordinate", ylab = "Y-coordinate")


# =============================================================================
# PART 4: MONITORING NETWORK OPTIMISATION
# =============================================================================

mean_krig_vars <- numeric(nrow(head$coords))

max_var_point <- which.min(
  sapply(1:nrow(head$coords), function(i) {
    temp_data <- list(coords = head$coords[-i,], data = head$data[-i])
    class(temp_data) <- "geodata"
    mean_krig_vars[i] <<- mean(
      krige.conv(temp_data,
                 locations = grid_locs,
                 krige     = krige.control(obj.model = chosen_model))$krige.var
    )
  })
)

cat("Most valuable point to retain: Point", max_var_point, "\n")
cat("Coordinates:", round(head$coords[max_var_point,], 2), "\n")
cat("Observed value:", head$data[max_var_point], "\n")

# Sorted variance plot
sorted_vars <- sort(mean_krig_vars)
plot(sorted_vars,
     type = "b", pch = 10, col = "blue",
     main = "Sorted Mean Kriging Variance by Point Removal",
     xlab = "Rank (1 = Most Valuable)", ylab = "Mean Kriging Variance")
abline(h = median(sorted_vars),          col = "red",       lty = 2)
abline(h = quantile(sorted_vars, 0.9),   col = "darkgreen", lty = 3)
legend("bottomright",
       legend = c("Variance", "Median", "90th percentile"),
       col    = c("blue", "red", "darkgreen"),
       lty    = c(1, 2, 3))

# Network map
n_points  <- length(mean_krig_vars)
n_select  <- ceiling(0.1 * n_points)
bottom_10 <- order(mean_krig_vars)[1:n_select]
top_10    <- order(mean_krig_vars, decreasing = TRUE)[1:n_select]

plot(head$coords,
     main = "Network Optimisation — Top & Bottom 10%",
     xlab = "X-coordinate", ylab = "Y-coordinate",
     pch = 1, col = "gray40")
points(head$coords[bottom_10,], pch = 2, col = "blue", cex = 1.5)
points(head$coords[top_10,],   pch = 0, col = "red",  cex = 1.5)
legend("topleft",
       legend = c("Regular", "High value — retain", "Low value — removable"),
       col    = c("gray40", "blue", "red"),
       pch    = c(1, 2, 0), pt.cex = 1.5)


# =============================================================================
# PART 5: SPATIAL POINT PROCESS ANALYSIS
# =============================================================================

data("hyytiala")

# Overview plot
plot(hyytiala,
     main = "Hyytiala Forest — Tree Species Locations",
     cols = c("red", "blue", "darkgrey", "brown4"),
     pch  = c(1, 2, 3, 4))

# Kernel density estimate
plot(density(hyytiala), main = "Kernel Density Estimate")
contour(density(hyytiala), add = TRUE)

# CSR quadrat test
quad_test <- quadrat.test(hyytiala, 3, 3)
plot(quad_test, main = "Quadrat Test for CSR")
print(quad_test)

# Species-level intensity
aspen_cells <- subset(hyytiala, marks == "aspen")
birch_cells <- subset(hyytiala, marks == "birch")
pine_cells  <- subset(hyytiala, marks == "pine")
rowan_cells <- subset(hyytiala, marks == "rowan")

par(mfrow = c(2, 2))
plot(density(aspen_cells), main = "Aspen Intensity")
plot(density(birch_cells), main = "Birch Intensity")
plot(density(pine_cells),  main = "Pine Intensity")
plot(density(rowan_cells), main = "Rowan Intensity")
par(mfrow = c(1, 1))

# Cross-type pair correlation function
pcf_cross <- pcfcross(hyytiala, i = "birch", j = "rowan")
plot(pcf_cross, main = "Cross-Type Pair Correlation: Birch vs Rowan")

# F and G summary statistics
unmarked <- unmark(hyytiala)
f_func   <- Fest(unmarked)
g_func   <- Gest(unmarked)

par(mfrow = c(1, 2))
plot(f_func, main = "F-function (Empty Space)")
plot(g_func, main = "G-function (Nearest Neighbour)")
par(mfrow = c(1, 1))


# =============================================================================
# PART 6: SPATIO-TEMPORAL POISSON PROCESS SIMULATION
# =============================================================================

# Thinning algorithm — intensity: lambda(x,y,t) = 100(x+y)^2 * (5-t)^2
simulate_spatiotemporal <- function(t_max = 1, x_max = 1, y_max = 1) {
  lambda     <- function(x, y, t) { 100 * (x + y)^2 * (5 - t)^2 }
  lambda_max <- lambda(x_max, y_max, 5)
  events     <- data.frame(x = numeric(0), y = numeric(0), t = numeric(0))

  for (i in 1:1000) {
    x_cand <- runif(1, 0, x_max)
    y_cand <- runif(1, 0, y_max)
    t_cand <- runif(1, 0, t_max)

    if (runif(1) < lambda(x_cand, y_cand, t_cand) / lambda_max) {
      events <- rbind(events, data.frame(x = x_cand, y = y_cand, t = t_cand))
    }
  }
  events[order(events$t), ]
}

set.seed(123)
events <- simulate_spatiotemporal(t_max = 3, x_max = 2, y_max = 2)

# Spatial distribution
plot(events$x, events$y,
     pch = 1, col = "blue",
     main = "Simulated Event Locations",
     xlab = "X", ylab = "Y")

# Temporal dynamics
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(events$t, events$x^2 + events$y^2,
     pch = 19, col = "purple", cex = 0.6,
     main = "Time vs Spatial Component",
     xlab = "Time", ylab = "Spatial Intensity")
plot(sort(events$t), 1:nrow(events),
     type = "s", col = "darkblue", lwd = 1.5,
     main = "Cumulative Events Over Time",
     xlab = "Time", ylab = "Count")
par(mfrow = c(1, 1))

# Empirical vs theoretical intensity
t_max <- 3; x_max <- 2; y_max <- 2
spatial_integral <- 100 * (x_max^3/3 + x_max^2*y_max + x_max*y_max^2)
lambda_theory    <- function(t) spatial_integral * (5 - t)

time_breaks        <- seq(0, t_max, length = 21)
time_counts        <- hist(events$t, breaks = time_breaks, plot = FALSE)$counts
temporal_intensity <- time_counts / diff(time_breaks)[1]

plot(time_breaks[-1] - diff(time_breaks)/2, temporal_intensity,
     type = "b", pch = 19, col = "black",
     main = "Temporal Intensity: Empirical vs Theoretical",
     xlab = "Time", ylab = "Intensity (events/unit time)",
     ylim = c(0, max(temporal_intensity, lambda_theory(0)) * 1.1))
curve(lambda_theory(x), col = "red", lwd = 2, add = TRUE)
legend("topright",
       c("Empirical", "Theoretical"),
       col = c("black", "red"),
       lwd = c(1, 2), pch = c(19, NA))

# Subwindow simulations
lambda_spatial <- function(x, y) { 400 * (x + y)^2 }
win <- owin(c(0.5, 1), c(0, 1))

set.seed(123); sim1 <- rpoispp(lambda_spatial, win = win)
set.seed(456); sim2 <- rpoispp(lambda_spatial, win = win)

par(mfrow = c(1, 3), mar = c(2, 2, 3, 1))
plot(sim1, main = "Simulation 1", cols = "blue")
plot(sim2, main = "Simulation 2", cols = "red")
plot(events$x, events$y,
     pch = 19, cex = 0.3, col = "gray",
     main = "Original Process",
     xlim = c(0, 2), ylim = c(0, 2))
rect(1, 0, 2, 1, border = "red", lwd = 2)
par(mfrow = c(1, 1))

# =============================================================================
# End of script
# Author: Nafisa Ahmad
# GitHub: https://github.com/namoumita
# LinkedIn: https://www.linkedin.com/in/nafisa-ahmad-957474199
# =============================================================================
