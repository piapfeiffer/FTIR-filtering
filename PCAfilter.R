# Pre-Processing for FTIR spectra
# filtering of non-informative variables
library(ggplot2)

PCAfilter <- function(X, threshold = 0.95){
  # performs PCA on matrix X and reconstructs from components needed to explain 
  # the variance given by the threshold
  
  X_scaled <- scale(X)
  center <- attr(X_scaled, "scaled:center")
  scale <- attr(X_scaled, "scaled:scale")
  
  # PCA on scaled and centered data
  res.PCA <- prcomp(X_scaled, center = FALSE, scale.=FALSE)
  
  # reconstruction
  k <- min(which(cumsum(var_explained(res.PCA))>threshold))
  Xhat <- res.PCA$x[, 1:k]%*%t(res.PCA$rotation[, 1:k])
  Xhat <- as.matrix(t((t( Xhat) * scale)) + rep(center, each = nrow(v)))

  # computation of reconstruction error
  RE <- matrix(data = NA, nrow= nrow(X), ncol = ncol(X))
  colnames(RE) <- colnames(X)
  for (i in 1:nrow(X)) {
    comp <- rbind(X[i,], Xhat[i,])
    res <- t(comp[1,] - comp[2,])^2
    RE[i, ] <- res
  }
  MRE <- apply(RE, 2, mean)
  
  # smoothing step
  MRE_smooth <- ksmooth(1:length(MRE), MRE, "normal", bandwidth = 20)
  MRE <- MRE_smooth$y
  threshold_non_inf <- quantile(MRE, prob = 0.95, na.rm = TRUE)
  index_noise <- which(abs(MRE)>threshold_non_inf)
  
  # diagnostic plot
  plot_df <- data.frame(cbind(wn = as.numeric(colnames(X_orig)), mre = MRE))
  diag_plot <- ggplot(data = plot_df, aes(x = wn, y = mre)) +
    theme_minimal() +
    theme(panel.grid.major = element_line(),panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA)) +
    geom_line() +
    geom_hline(yintercept = threshold, linetype="dashed", color="red") +
    scale_x_reverse() +
    xlab(expression(paste("Wavenumber in ", bold(cm ^-1)))) +
    ylab("Mean reconstruction error")

  show(diag_plot)
  return(list(index = index_noise, X = Xhat))
}
