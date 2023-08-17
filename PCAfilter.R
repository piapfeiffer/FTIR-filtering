# Pre-Processing for FTIR spectra
# filtering of non-informative variables
library(ggplot2)
library(dplyr)

var_explained <- function(prcomp_obj) {
  a <- prcomp_obj$sdev^2
  a / sum(a)
}

PCAfilter <- function(X, threshold = 0.95){
  # performs PCA on matrix X and reconstructs from components needed to explain 
  # the variance given by the threshold
  
  X_scaled <- scale(X)
  center <- attr(X_scaled, "scaled:center")
  scale <- attr(X_scaled, "scaled:scale")
  
  const_cols <- apply(X_scaled %>% select(starts_with("x")),2,mad)>0
  X_scaled <- X_scaled[, const_cols]
  center <- center[const_cols]
  scale <- scale[const_cols]
  attr(X_scaled, "scaled:center") <- center
  attr(X_scaled, "scaled:scale") <- scale
  
  # PCA on scaled and centered data
  res.PCA <- prcomp(X_scaled, center = FALSE, scale.=FALSE)
  k <- min(which(cumsum(var_explained(res.PCA))>threshold))
    
  # reconstruction
  Xhat <- res.PCA$x[, 1:k]%*%t(res.PCA$rotation[, 1:k])
  Xhat <- as.matrix(t((t( Xhat) * scale)) + rep(center, each = nrow(Xhat)))

  # computation of reconstruction error
  RE <- matrix(data = NA, nrow= nrow(X), ncol = ncol(X[,const_cols]))
  colnames(RE) <- colnames(X[,const_cols])
  for (i in 1:nrow(X)) {
    comp <- rbind(X[i,const_cols], Xhat[i,])
    res <- t(comp[1,] - comp[2,])^2
    RE[i, ] <- res
  }

  MRE <- apply(RE, 2, mean)
  
  # smoothing step
  MRE_smooth <- ksmooth(1:length(MRE), MRE, "normal", bandwidth = 20)
  MRE <- MRE_smooth$y
  threshold_non_inf <- quantile(MRE, prob = 0.925, na.rm = TRUE)
  index_noise <- which(abs(MRE)>threshold_non_inf)
  
  # diagnostic plot
  wn <- as.numeric(gsub("_", ".", 
                        gsub("x", "", 
                             colnames(X[, const_cols]))))
  plot_df <- data.frame(cbind(wn = wn, mre = MRE))
  diag_plot <- ggplot(data = plot_df, aes(x = wn, y = mre)) +
    theme_minimal() +
    theme(panel.grid.major = element_line(),panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA)) +
    geom_line() +
    geom_hline(yintercept = threshold_non_inf, linetype="dashed", color="red") +
    scale_x_reverse() +
    xlab(expression(paste("Wavenumber in ", bold(cm ^-1)))) +
    ylab("Mean reconstruction error")
  
  show(diag_plot)
  
  return(list(index = colnames(X[, const_cols])[index_noise], X = Xhat))
}
