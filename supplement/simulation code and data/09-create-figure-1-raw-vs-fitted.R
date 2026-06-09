library("tidyverse")

.plot_theme_mod <- function(plot_base, text_size = 9, 
                            lab_text_size = 8, space_legend = 0.1,
                            neg_legend_margin = -0.35){
  # Function modifies plot aesthetics to make better use of plot real estate
  # neg legend margin pushes the legend further into the direction of the plot
  # to reduce whitespace between legend and plot. 
  out_plot <- plot_base + 
    theme(text = element_text(size = text_size, family = "Ubuntu Light"),
          plot.margin = unit(x = c(0.03, 0.03, 0.03, 0.03), units = "cm"),
          legend.margin=margin(t=neg_legend_margin, r=0.03, b=0.03, l=0.03, 
                               unit="cm"),
          legend.position="bottom",
          legend.title = element_text(size = lab_text_size), 
          legend.text  = element_text(size = lab_text_size),
          legend.key.size = unit(space_legend, "lines"),
          panel.grid.minor = element_line(size = 0.1), 
          panel.grid.major = element_line(size = 0.2)) +
    scale_color_manual(values = c("black","#D55E00", "#CC79A7", 
                                  "#009E73", "#56B4E9")) + 
    scale_shape_manual(values = c(0,1,2,5,6)) +
    scale_linetype_manual(values = c("dotted","dotdash","dashed"))
}

## Load tables with parameters for fitted models -------------------------------
fitted_parameters_anaero <- 
  readRDS("output/fitted-parameters-anaero.RDS")
fitted_parameters_bifido <- 
  readRDS("output/fitted-parameters-bifido.RDS")

## Load raw data ratio tables --------------------------------------------------
ratio_table_anaero <- readRDS("output/ratio_table_anaero.RDS")
ratio_table_bifido <- readRDS("output/ratio_table_bifido.RDS")

## Make Empirical Pairwise Proportion Histograms -------------------------------
## Select Ratio Pairs
tmp1 <- ratio_table_bifido %>%
  filter(., numerator %in% c("Veillonella", "Lactobacillus")) %>%
  mutate(ratio_id = as.integer(ratio_id))
tmp2 <- ratio_table_anaero %>%
  filter(., numerator %in% c("Streptococcus", "Actinomyces")) %>%
  mutate(ratio_id = as.integer(ratio_id))
plotable <- bind_rows(tmp1, tmp2)

plotable <- plotable %>%
  mutate(., taxon_pair = paste(numerator, "to", denominator),
         pairwise_proportion = (num_count / (num_count + denom_count)))

tmp_plot_1 <- ggplot(plotable) + 
  geom_histogram(aes(x = pairwise_proportion), 
                 alpha = 1, size = 0.2, color = "grey45",
                 breaks = seq(0,1, by = 0.05)) +
  facet_wrap(~taxon_pair, ncol = 4) + 
  labs(title = "Empirical Pairwise Proportion histogram by Baseline Pair",
       caption = "64 pairwise zero-count samples removed.") +
  xlab(label = "Empirical Pairwise Proportion") +
  ylab(label = "Density")
tmp_plot_1 <- .plot_theme_mod(tmp_plot_1, text_size = 6)


# Fitted Beta-Binomial Density Plots ----
tmp3 <- filter(fitted_parameters_bifido, ratio_id %in% unique((tmp1$ratio_id))) %>%
  left_join(., unique(select(tmp1, ratio_id, numerator, denominator)), by = "ratio_id")
tmp4 <- filter(fitted_parameters_anaero, ratio_id %in% unique(as.integer(tmp2$ratio_id)))  %>%
  left_join(., unique(select(tmp2, ratio_id, numerator, denominator)), by = "ratio_id")
base_settings <- bind_rows(tmp3, tmp4) %>% 
  mutate(., taxon_pair = paste(numerator, "to", denominator))

x <- seq(0.001, 0.999, length = 1000)
tmp_draws <- vector(mode = "list",  dim(base_settings)[1])
for (i in 1:dim(base_settings)[1]){
  alpha <- base_settings$alpha[i]
  beta <- base_settings$beta[i]
  tmp_draws[[i]] <- 
    tibble(taxon_pair = base_settings$taxon_pair[i],
           x = x, 
           probability_density = dbeta(x, alpha, beta))
}
plotable2 <- bind_rows(tmp_draws, .id = "id")
tmp_plot_2 <- ggplot(plotable2, aes(x, probability_density)) + 
  geom_line(alpha = 1, size = 1, color = "grey45") + 
  labs(title="Fitted Beta Distributions by Baseline Pair",
       caption = "Density plotted for proportion range: 0.001 to 0.999.") + 
  labs(x= "Proportion", 
       y= "Probability Density Function")  +
  xlab(label = "Pairwise Proportion") +
  ylab(label = "Density") +
  facet_wrap(~taxon_pair, ncol = 4)
tmp_plot_2 <- .plot_theme_mod(tmp_plot_2, text_size = 6)

library(gridExtra)
tmp_plot <- arrangeGrob(tmp_plot_1, tmp_plot_2, ncol = 1)
ggsave(filename = "figures/empirical-vs-fitted.pdf", tmp_plot,
       device = cairo_pdf, width = 14.25, height = 7.5, units = "cm", dpi = 500)
