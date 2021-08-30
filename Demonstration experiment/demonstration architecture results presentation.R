# This script will present the results obtained from the architecture demonstration
# experiment

# Load packages
library("ggplot2")
library("Cairo")

# Add a Cairo window for better rendered images
CairoWin()

# Code snippet to save image in the Cairo window
ggsave(boxplots_sens, filename = 'Figures_demo/somehting.png', type = 'cairo')

# Load in processed results
# "method"      "LFC"         "Fraction.DE" "n"           "sensitivity" "FDR"         "runtime.s"  
arch_demo.df <- readRDS(file ="benchmark_results/architecture_demonstration_results.RDS")

# Convert grouping variables to factors
arch_demo.df$method <- as.factor(arch_demo.df$method)
arch_demo.df$LFC <- as.factor(arch_demo.df$LFC)
arch_demo.df$Fraction.DE <- as.factor(arch_demo.df$Fraction.DE)
arch_demo.df$n <- as.factor(arch_demo.df$n)
levels(arch_demo.df$method)[3] <- "Wilcoxon"

# Univariate analysis
# Boxplots: fdr
boxplots_fdr <- ggplot(data=subset(arch_demo.df, 
                                   subset = method %in% c("DESeq2", "EdgeR", "Wilcoxon")), 
                       aes(y=method, x=FDR, fill = LFC)) +
  geom_boxplot() +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 0.65) +
  facet_grid(n ~ Fraction.DE, labeller = label_both, scales = "fixed") +
  labs(y = "Analysis method") +
  theme(axis.text.x = element_text( size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text( size = 14, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text( size = 14, hjust = .5, vjust = .5, face = "bold"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black", angle = 0))
boxplots_fdr

# Boxplots: sensitivity
boxplots_sens <- ggplot(data=arch_demo.df, aes(y=method, x=sensitivity, fill = LFC)) +
  geom_boxplot() +
  facet_grid(n ~ Fraction.DE, labeller = label_both) +
  labs(y = "Analysis method", x="Sensitivity") +
  theme(axis.text.x = element_text( size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text( size = 14, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text( size = 14, hjust = .5, vjust = .5, face = "bold"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black", angle = 0))
boxplots_sens

# Histograms: fdr
histogram_fdr <- ggplot(data=subset(arch_demo.df, 
                                   subset = Fraction.DE %in% c(0.05, 0.4) & 
                                     n %in% c(5, 100)),
                       aes(x=FDR, fill = LFC)) +
  stat_bin(position = "identity", alpha = 0.75, bins = 30) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 0.65) +
  facet_grid( Fraction.DE + n ~ method, labeller = label_both, scales = "free") +
  labs() +
  theme(axis.text.x = element_text( size = 12, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text( size = 12, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 12, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text( size = 12, hjust = .5, vjust = .5, face = "bold"),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black", angle = -90))
histogram_fdr

# Histograms: sensitivity
histogram_sens <- ggplot(data=subset(arch_demo.df, 
                                    subset = Fraction.DE %in% c(0.05, 0.4) & 
                                      n %in% c(5, 100)),
                        aes(x=sensitivity, fill = LFC)) +
  stat_bin(position = "identity", alpha = 0.75, bins = 30) +
  facet_grid( Fraction.DE + n ~ method, labeller = label_both, scales = "free") +
  labs() +
  theme(axis.text.x = element_text( size = 12, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text( size = 12, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 12, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text( size = 12, hjust = .5, vjust = .5, face = "bold"),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black", angle = -90))
histogram_sens

# Histograms runtime
histogram_runtime <- ggplot(data=arch_demo.df,
                         aes(x=runtime.s, fill = n)) +
  stat_bin(position = "identity", alpha = 0.5, bins = 50) +
  facet_grid( Fraction.DE + LFC ~ method, labeller = label_both, scales = "free") +
  labs( x = "Runtime (s)") +
  theme(axis.text.x = element_text( size = 12, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text( size = 12, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 16, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text( size = 16, hjust = .5, vjust = .5, face = "bold"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 13, colour = "black", angle = -90))
histogram_runtime


# Calculate mean and quantile based intervals for each quantity
mean_and_quantile <- aggregate(.~ method + LFC + Fraction.DE + n,
          arch_demo.df,
          FUN = function(x){
            lower_quantile <- quantile(x, 0.05)
            upper_quantile <- quantile(x, 0.95)
            center_mean <- mean(x)
            
            c(lower_quantile, center_mean, upper_quantile)
          })
mean_and_quantile$sensitivity <- as.data.frame(mean_and_quantile$sensitivity)
mean_and_quantile$FDR <- as.data.frame(mean_and_quantile$FDR)
mean_and_quantile$runtime.s <- as.data.frame(mean_and_quantile$runtime.s)

colnames(mean_and_quantile$sensitivity) <- c("5%", "mean", "95%")
colnames(mean_and_quantile$FDR) <- c("5%", "mean", "95%")
colnames(mean_and_quantile$runtime.s) <- c("5%", "mean", "95%")

mean_and_quantile$x.labels <- paste(paste("LFC", mean_and_quantile$LFC, sep = " = "), 
                                    paste("pDE", mean_and_quantile$Fraction.DE, sep = " = "),
                                    paste("n", mean_and_quantile$n, sep = " = "),
                                          sep = ", ")

# Plot mean and quantile based intervals 
mean_and_quantile_plot <- ggplot(data=subset(mean_and_quantile, 
                                   subset = method %in% c("DESeq2", "EdgeR", "Wilcoxon")), 
                       aes(x=FDR[["mean"]], y = x.labels, color = x.labels)) +
  geom_pointrange(aes(xmin=FDR[["5%"]], xmax=FDR[["95%"]])) +
  #geom_vline(xintercept = 0.05, color = "red", linetype = 'dashed') +
  facet_grid(.~ method, scales = "fixed") +
  labs(y = "", x ="Runtime (s)") +
  theme(axis.text.x = element_text( size = 16, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text( size = 14, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text( size = 14, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text( size = 14, hjust = .5, vjust = .5, face = "bold", angle = 0),
        strip.text.x = element_text(size = 16, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black", angle = 0),
        legend.position = "none")
mean_and_quantile_plot

# Runtime
# Calculate mean and quantile based intervals for runtime
mean_and_quantile_runtime <- aggregate(.~ method + LFC + Fraction.DE + n,
                               arch_demo.df,
                               FUN = function(x){
                                 lower_quantile <- quantile(x[-1], 0.05)
                                 upper_quantile <- quantile(x[-1], 0.95)
                                 center_mean <- mean(x[-1])
                                 
                                 c(lower_quantile, center_mean, upper_quantile)
                               })
mean_and_quantile_runtime$runtime.s <- as.data.frame(mean_and_quantile_runtime$runtime.s)

colnames(mean_and_quantile_runtime$runtime.s) <- c("5%", "mean", "95%")

mean_and_quantile_runtime$x.labels <- paste(paste("LFC", mean_and_quantile_runtime$LFC, sep = " = "), 
                                    paste("pDE", mean_and_quantile_runtime$Fraction.DE, sep = " = "),
                                    paste("n", mean_and_quantile_runtime$n, sep = " = "),
                                    sep = ", ")

# Plot mean and quantile based intervals
mean_and_quantile_runtime_plot <- ggplot(data=subset(mean_and_quantile_runtime, 
                                             subset = method %in% c("DESeq2", "EdgeR", "Wilcoxon")), 
                                 aes(x=runtime.s[["mean"]], y = x.labels, color = x.labels)) +
  geom_pointrange(aes(xmin=runtime.s[["5%"]], xmax=runtime.s[["95%"]])) +
  #geom_vline(xintercept = 0.05, color = "red", linetype = 'dashed') +
  facet_grid(.~ method, scales = "fixed") +
  labs(y = "", x ="Runtime (s)") +
  theme(axis.text.x = element_text( size = 16, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text( size = 14, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text( size = 14, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text( size = 14, hjust = .5, vjust = .5, face = "bold", angle = 0),
        strip.text.x = element_text(size = 16, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black", angle = 0),
        legend.position = "none")
mean_and_quantile_runtime_plot


# Bivariate analysis
# Sensitivity vs. FDR
sens_vs_fdr <- ggplot(data=subset(arch_demo.df, 
                                  subset = Fraction.DE %in% c(0.05,0.15, 0.4) & 
                                    n %in% c(5,25, 100)), aes(y=sensitivity, x=(1-FDR), color = method, size = LFC, shape = LFC)) +
  geom_point(alpha = 0.75) +
  scale_shape_manual(values=c(3, 20)) +
  scale_size_manual(values=c(2,1)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", size = 0.3) +
  facet_grid(n ~ Fraction.DE, labeller = label_both) +
  labs() +
  theme(axis.text.x = element_text( size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text( size = 14, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text( size = 14, hjust = .5, vjust = .5, face = "bold"),
        strip.text.x = element_text(size = 12, colour = "black"),
        strip.text.y = element_text(size = 12, colour = "black", angle = 0))
sens_vs_fdr

# Plot mean and quantile based intervals for both FDR and sensitivity at once
# Re-calculate the mean and its quantiles, this is done because only a subset of
# factors will be used in order not to clutter then final plot
mean_and_quantile <- aggregate(.~ method + LFC + Fraction.DE + n,
                               subset(arch_demo.df, subset = 
                                        Fraction.DE %in% c(0.05,0.15, 0.4) & 
                                        n %in% c(100)),
                               FUN = function(x){
                                 lower_quantile <- quantile(x, 0.05)
                                 upper_quantile <- quantile(x, 0.95)
                                 center_mean <- mean(x)
                                 
                                 c(lower_quantile, center_mean, upper_quantile)
                               })
mean_and_quantile$sensitivity <- as.data.frame(mean_and_quantile$sensitivity)
mean_and_quantile$FDR <- as.data.frame(mean_and_quantile$FDR)
mean_and_quantile$runtime.s <- as.data.frame(mean_and_quantile$runtime.s)

colnames(mean_and_quantile$sensitivity) <- c("5%", "mean", "95%")
colnames(mean_and_quantile$FDR) <- c("5%", "mean", "95%")
colnames(mean_and_quantile$runtime.s) <- c("5%", "mean", "95%")

mean_and_quantile$x.labels <- paste(paste("LFC", mean_and_quantile$LFC, sep = " = "), 
                                    paste("pDE", mean_and_quantile$Fraction.DE, sep = " = "),
                                    paste("n", mean_and_quantile$n, sep = " = "),
                                    sep = ", ")

# Plot the mean and quantile based intervals for both FDR and sensitivity
mean_and_quantile_plot <- ggplot(data=subset(mean_and_quantile, 
                                             subset = method %in% c("DESeq2", "EdgeR", "Wilcoxon")), 
                                 aes(x=FDR[["mean"]], y = sensitivity[["mean"]], color = Fraction.DE, size=LFC)) +
  geom_errorbar(aes(xmin=FDR[["5%"]], xmax=FDR[["95%"]]), size = 0.5) +
  geom_errorbar(aes(ymin=sensitivity[["5%"]], ymax=sensitivity[["95%"]]), size = 0.5) +
  geom_point()+
  #geom_vline(xintercept = 0.05, color = "red", linetype = 'dashed') +
  facet_grid(.~ method, scales = "fixed") +
  labs(y = "Sensitivity", x ="FDR", color="DE fraction") +
  geom_vline(xintercept = 0.05, color = "red", linetype="dashed")+
  theme(axis.text.x = element_text( size = 16, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text( size = 14, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text( size = 14, hjust = .5, vjust = .5, face = "bold", angle = 0),
        strip.text.x = element_text(size = 16, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black", angle = 0))+
  scale_size_manual(values = c(2,4))
  
mean_and_quantile_plot
