rm(list = ls())
library(tidyverse)
library(ggthemes)

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


## Make Parametric Output plots for Bifidobacterium reference run --------------
## Load Data
output_summary_table <- readRDS("output/bifido-parametric-summary.RDS")

## Create Print friendly factor levels 
tmp_levels <- paste(c("Low Dispersion", "Median Dispersion", "High Dispersion"),
                    " (≈", 
                    sort(unique(round(output_summary_table$disp_control,2))),
                    ")", sep = "")
tmp_disp_labels <- tibble(disp_control = 
                            sort(unique(output_summary_table$disp_control)),
                          disp = factor(tmp_levels, levels = tmp_levels))
tmp_levels <- c("1/999", "1/99", "1/9", "1/1")
tmp_odds_labels <- tibble(odds_control = 
                            sort(unique(output_summary_table$odds_control)),
                          odds = factor(tmp_levels, 
                                        levels = tmp_levels))

tmp_levels <- paste("Baseline Odds = ", c("1/999", "1/99", "1/9", "1/1"), sep = "")
tmp_odds_labels_long <- 
  tibble(odds_control = sort(unique(output_summary_table$odds_control)),
         odds = factor(tmp_levels, levels = tmp_levels))

tmp_effect_labels <- output_summary_table %>%
  mutate(., diff = abs(mu_control - mu_new)) %>%
  mutate(., effect = paste(diff, " (", round(odds_ratio, 2), ")", 
                           sep = "")) %>%
  arrange(., diff) %>%
  mutate(., effect = factor(effect, levels = unique(effect))) %>%
  select(., scenario_id, diff, effect)

## Coverage Plot for Null Scenario 
tmp_plot <- filter(output_summary_table, 
                   method %in% c("bbglm", "lrlm", "lrlm_2")) %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels, by = "odds_control") %>%
  mutate(., diff = round(abs(mu_control - mu_new), 3)) %>%
  filter(., diff == 0) %>%
  arrange(., odds_control, disp_control) %>%
  mutate(., odds_control = round(odds_control, 3), 
         disp_control = round(disp_control, 2)) %>%
  mutate(., odds_control = as.factor(odds_control)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_point(aes(y = coverage_conditional, x = odds, 
                 color = method, pch = method), 
             alpha =0.8, size = 1.25) +
  ylim(c(0,1)) +
  ylab(label = "95% Confidence Interval Coverage") +
  xlab(label = "Baseline Odds") +
  labs(color = "Method:", pch = "Method:") +
  #theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(n_samples~disp) +
  geom_hline(yintercept = 0.95, lty = "dotted", color = "black", size = 0.5)
tmp_plot <- .plot_theme_mod(tmp_plot)  
ggsave(filename = "figures/parametric-coverage-bifido-null-scenario.pdf", tmp_plot, 
       device = cairo_pdf, width = 14.25, height = 12, units = "cm", dpi = 500)
rm(list = c("tmp_plot"))


# Coverage - Various Scenarios 
tmp_plot <- filter(output_summary_table, 
                   method %in% c("bbglm", "lrlm", "lrlm_2")) %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels, by = "odds_control") %>%
  left_join(., tmp_effect_labels, by = "scenario_id") %>%
  filter(., disp == "Median Dispersion (≈0.46)") %>%
  mutate(., diff = round(abs(mu_control - mu_new), 3)) %>%
  arrange(., odds_control, disp_control) %>%
  mutate(., odds_control = round(odds_control, 3), 
         disp_control = round(disp_control, 2)) %>%
  mutate(., odds_control = as.factor(odds_control)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_point(aes(y = coverage_conditional, x = effect, color = method, 
                 pch = method), 
             alpha =0.8, size = 1.25) +
  ylim(c(0,1)) +
  ylab(label = "95% Confidence Interval Coverage") +
  xlab(label = "Absolute Difference in Success Probability Between Groups (Odds-Ratio)") +
  labs(color = "Method:", pch = "Method:") +
  facet_grid(n_samples~odds, scales = "free") +
  geom_hline(yintercept = 0.95, lty = "dotted", color = "black", size = 0.5) +
  theme(axis.text.x = element_text(angle = 90))
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 8)  
ggsave(filename = "figures/parametric-coverage-bifido-median-dispersion-scenarios.pdf", 
       tmp_plot, device = cairo_pdf, width = 14.25, height = 12, units = "cm", 
       dpi = 500)
rm(list = c("tmp_plot"))


## Nominal Error Rates
tmp_plot <- filter(output_summary_table) %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels, by = "odds_control") %>%
  mutate(., diff = round(abs(mu_control - mu_new), 3)) %>%
  filter(., diff == 0) %>%
  arrange(., odds_control, disp_control) %>%
  mutate(., odds_control = round(odds_control, 3), 
         disp_control = round(disp_control, 2)) %>%
  mutate(., odds_control = as.factor(odds_control)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_point(aes(y = power_conditional, x = odds, color = method, pch = method), 
             alpha =0.8, size = 1.25) +
  ylim(c(0,0.1)) +
  ylab(label = "Type-I Error Rate") +
  xlab(label = "Baseline Odds") +
  labs(color = "Method:", pch = "Method:") +
  facet_grid(n_samples~disp) +
  geom_hline(yintercept = 0.05, lty = "dotted", color = "black", size = 0.5)
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 8)  
ggsave(filename = "figures/parametric-type-1-error-bifido.pdf", 
       tmp_plot, device = cairo_pdf, width = 14.25, height = 8, units = "cm", 
       dpi = 500)
rm(list = c("tmp_plot"))


## Power - Median Dispersion across all sample sizes & reference scenarios
tmp_plot <- output_summary_table %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels_long, by = "odds_control") %>%
  left_join(., tmp_effect_labels, by = "scenario_id") %>%
  filter(., disp == "Median Dispersion (≈0.46)") %>%
  mutate(., diff = abs(mu_control - mu_new)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_point(aes(y = power_conditional, x = effect, color = method, pch = method), 
             alpha =0.6, size = 1) +
  geom_line(aes(y = power_conditional, x = effect, color = method, group = method), 
            alpha =0.6, size = 0.4, lty = "solid") +
  geom_hline(yintercept = 0.8, color = "black", lty = "dotted", size = 0.3) +
  scale_y_continuous(breaks = c(0.00, 0.20, 0.40, 0.60, 0.80, 1.00)) +
  ylab(label = "Power (1 - Type-II Error Rate)") +
  xlab(label = "Absolute Difference in Success Probability Between Groups (Odds-Ratio)") +
  labs(color = "Method:", pch = "Method:") +
  facet_grid(n_samples ~ odds, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90))
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 9)  
ggsave(filename = "figures/parametric-power-bifido-median-dispersion.pdf", tmp_plot, 
       device = cairo_pdf, width = 14.25, height = 14.25, 
       units = "cm", dpi = 500)
rm(list = c("tmp_plot"))

## Power - All Dispersion levels across all sample sizes & reference scenarios 
tmp_plot <- output_summary_table %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels_long, by = "odds_control") %>%
  left_join(., tmp_effect_labels, by = "scenario_id") %>%
  mutate(., diff = abs(mu_control - mu_new)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_line(aes(y = power_conditional, x = effect, color = method, 
                group = interaction(method, disp)), 
            alpha =0.5, size = 0.35) +
  geom_point(aes(y = power_conditional, x = effect, color = method, 
                 group = interaction(method, disp), 
                 pch = disp), 
             alpha =0.8, size = 0.8) +
  geom_hline(yintercept = 0.8, color = "black", lty = "dotted", size = 0.3) +
  ylab(label = "Power (1 - Type-II Error Rate)") +
  xlab(label = "Absolute Difference in Success Probability Between Groups (Odds-Ratio)") +
  labs(color = "Method:", pch = "Dispersion:") +
  facet_grid(n_samples ~ odds, scales = "free_x") +
  scale_y_continuous(breaks = c(0.00, 0.20, 0.40, 0.60, 0.80, 1.00)) +
  theme(axis.text.x = element_text(angle = 90))
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 9, lab_text_size = 7)  +
  scale_shape_manual(values = c(1,3,4))
ggsave(filename = "figures/parametric-power-bifido-all.pdf", tmp_plot, 
       device = cairo_pdf, width = 18.00, height = 14.25, 
       units = "cm", dpi = 500)





## Make Parametric Output plots for Anaerostipes reference run -----------------
output_summary_table <- readRDS("output/anaero-parametric-summary.RDS")

## Create Print friendly factor levels 
tmp_levels <- paste(c("Low Dispersion", "Median Dispersion", "High Dispersion"),
                    " (≈", 
                    sort(unique(round(output_summary_table$disp_control,2))),
                    ")", sep = "")
tmp_disp_labels <- tibble(disp_control = 
                            sort(unique(output_summary_table$disp_control)),
                          disp = factor(tmp_levels, levels = tmp_levels))
tmp_levels <- c("1/999", "1/99", "1/9", "1/1")
tmp_odds_labels <- tibble(odds_control = 
                            sort(unique(output_summary_table$odds_control)),
                          odds = factor(tmp_levels, 
                                        levels = tmp_levels))
tmp_levels <- paste("Baseline Odds = ", c("1/999", "1/99", "1/9", "1/1"), 
                    sep = "")
tmp_odds_labels_long <- 
  tibble(odds_control = sort(unique(output_summary_table$odds_control)),
         odds = factor(tmp_levels, levels = tmp_levels))
tmp_effect_labels <- output_summary_table %>%
  mutate(., diff = abs(mu_control - mu_new)) %>%
  mutate(., effect = paste(diff, " (", round(odds_ratio, 2), ")", 
                           sep = "")) %>%
  arrange(., diff) %>%
  mutate(., effect = factor(effect, levels = unique(effect))) %>%
  select(., scenario_id, diff, effect)

## Coverage Plot for Null Scenario 
tmp_plot <- filter(output_summary_table, 
                   method %in% c("bbglm", "lrlm", "lrlm_2")) %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels, by = "odds_control") %>%
  mutate(., diff = round(abs(mu_control - mu_new), 3)) %>%
  filter(., diff == 0) %>%
  arrange(., odds_control, disp_control) %>%
  mutate(., odds_control = round(odds_control, 3), 
         disp_control = round(disp_control, 2)) %>%
  mutate(., odds_control = as.factor(odds_control)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_point(aes(y = coverage_conditional, x = odds, 
                 color = method, pch = method), 
             alpha =0.8, size = 1.25) +
  ylim(c(0,1)) +
  ylab(label = "95% Confidence Interval Coverage") +
  xlab(label = "Baseline Odds") +
  labs(color = "Method:", pch = "Method:") +
  #theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(n_samples~disp) +
  geom_hline(yintercept = 0.95, lty = "dotted", color = "black", size = 0.5)
tmp_plot <- .plot_theme_mod(tmp_plot)  
ggsave(filename = "figures/parametric-coverage-anaero-null-scenario.pdf", tmp_plot, 
       device = cairo_pdf, width = 14.25, height = 12, units = "cm", dpi = 500)
rm(list = c("tmp_plot"))


# Coverage - Various Scenarios 
tmp_plot <- filter(output_summary_table, 
                   method %in% c("bbglm", "lrlm", "lrlm_2")) %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels, by = "odds_control") %>%
  left_join(., tmp_effect_labels, by = "scenario_id") %>%
  filter(., disp == "Median Dispersion (≈5.39)") %>%
  mutate(., diff = round(abs(mu_control - mu_new), 3)) %>%
  arrange(., odds_control, disp_control) %>%
  mutate(., odds_control = round(odds_control, 3), 
         disp_control = round(disp_control, 2)) %>%
  mutate(., odds_control = as.factor(odds_control)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_point(aes(y = coverage_conditional, x = effect, color = method, pch = method), 
             alpha =0.8, size = 1.25) +
  ylim(c(0,1)) +
  ylab(label = "95% Confidence Interval Coverage") +
  xlab(label = "Absolute Difference in Success Probability Between Groups (Odds-Ratio)") +
  labs(color = "Method:", pch = "Method:") +
  facet_grid(n_samples~odds, scales = "free") +
  geom_hline(yintercept = 0.95, lty = "dotted", color = "black", size = 0.5) +
  theme(axis.text.x = element_text(angle = 90))
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 8)  
ggsave(filename = "figures/parametric-coverage-anaero-median-dispersion-scenarios.pdf", 
       tmp_plot, device = cairo_pdf, width = 14.25, height = 12, units = "cm", 
       dpi = 500)
rm(list = c("tmp_plot"))


## Nominal Error Rates
tmp_plot <- filter(output_summary_table) %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels, by = "odds_control") %>%
  mutate(., diff = round(abs(mu_control - mu_new), 3)) %>%
  filter(., diff == 0) %>%
  arrange(., odds_control, disp_control) %>%
  mutate(., odds_control = round(odds_control, 3), 
         disp_control = round(disp_control, 2)) %>%
  mutate(., odds_control = as.factor(odds_control)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_point(aes(y = power_conditional, x = odds, color = method, pch = method), 
             alpha =0.8, size = 1.25) +
  ylim(c(0,0.1)) +
  ylab(label = "Type-I Error Rate") +
  xlab(label = "Baseline Odds") +
  labs(color = "Method:", pch = "Method:") +
  facet_grid(n_samples~disp) +
  geom_hline(yintercept = 0.05, lty = "dotted", color = "black", size = 0.5)
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 8)  
ggsave(filename = "figures/parametric-type-1-error-anaero.pdf", 
       tmp_plot, device = cairo_pdf, width = 14.25, height = 8, units = "cm", 
       dpi = 500)
rm(list = c("tmp_plot"))


## Power - Median Dispersion across all sample sizes & reference scenarios
tmp_plot <- output_summary_table %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels_long, by = "odds_control") %>%
  left_join(., tmp_effect_labels, by = "scenario_id") %>%
  filter(., disp == "Median Dispersion (≈5.39)") %>%
  mutate(., diff = abs(mu_control - mu_new)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_point(aes(y = power_conditional, x = effect, color = method, pch = method), 
             alpha =0.6, size = 1) +
  geom_line(aes(y = power_conditional, x = effect, color = method, group = method), 
            alpha =0.6, size = 0.4, lty = "solid") +
  geom_hline(yintercept = 0.8, color = "black", lty = "dotted", size = 0.3) +
  scale_y_continuous(breaks = c(0.00, 0.20, 0.40, 0.60, 0.80, 1.00)) +
  ylab(label = "Power (1 - Type-II Error Rate)") +
  xlab(label = "Absolute Difference in Success Probability Between Groups (Odds-Ratio)") +
  labs(color = "Method:", pch = "Method:") +
  facet_grid(n_samples ~ odds, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90))
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 9)  
ggsave(filename = "figures/parametric-power-anaero-median-dispersion.pdf", tmp_plot, 
       device = cairo_pdf, width = 14.25, height = 14.25, 
       units = "cm", dpi = 500)
rm(list = c("tmp_plot"))

## Power - All Dispersion levels across all sample sizes & reference scenarios 
tmp_plot <- output_summary_table %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels_long, by = "odds_control") %>%
  left_join(., tmp_effect_labels, by = "scenario_id") %>%
  mutate(., diff = abs(mu_control - mu_new)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_line(aes(y = power_conditional, x = effect, color = method, 
                group = interaction(method, disp)), 
            alpha =0.5, size = 0.35) +
  geom_point(aes(y = power_conditional, x = effect, color = method, 
                 group = interaction(method, disp), 
                 pch = disp), 
             alpha =0.8, size = 0.8) +
  geom_hline(yintercept = 0.8, color = "black", lty = "dotted", size = 0.3) +
  ylab(label = "Power (1 - Type-II Error Rate)") +
  xlab(label = "Absolute Difference in Success Probability Between Groups (Odds-Ratio)") +
  labs(color = "Method:", pch = "Dispersion:") +
  facet_grid(n_samples ~ odds, scales = "free_x") +
  scale_y_continuous(breaks = c(0.00, 0.20, 0.40, 0.60, 0.80, 1.00)) +
  theme(axis.text.x = element_text(angle = 90))
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 9, lab_text_size = 7)  +
  scale_shape_manual(values = c(1,3,4))
ggsave(filename = "figures/parametric-power-anaero-all.pdf", tmp_plot, 
       device = cairo_pdf, width = 18.00, height = 14.25, 
       units = "cm", dpi = 500)
rm(list = c("tmp_plot"))


# Make Sampling Based Output plots for Bifidobacterium reference run -----------
output_summary_table <- readRDS("output/bifido-sampling-summary.RDS")

tmp_plot <- output_summary_table %>%
  arrange(., ratio_id) %>%
  mutate(., ratio_id = as.factor(ratio_id)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("Sample Size =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_histogram(aes(y = ..density.., x = type_I_error_rate_conditional, 
                     color = method, fill = method), 
                 alpha =0.8, size = 0.2, position = "identity", 
                 color = "black", breaks = seq(0,0.1, by = 0.005)) +
  geom_density(aes(x = type_I_error_rate_conditional, color = method, 
                   fill = method), 
               alpha=0.2, bw = 0.005, kernel = "gaussian") +
  geom_vline(xintercept = 0.05, color = "black", lty = "dotted", size = 0.5) + 
  xlab(label = "Type-I Error Rate") +
  ylab(label = "Density") +
  labs(color = "Method:", fill = "Method:") +
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_grid(method ~ n_samples, scales = "free_y")
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 8, lab_text_size = 6, 
                            space_legend = 0.5) +
  scale_color_manual(values = c("black","#D55E00", "#009E73", "#56B4E9")) + 
  scale_fill_manual(values = c("black","#D55E00", "#009E73", "#56B4E9"))

ggsave(filename = "figures/sampling-type-1-error-bifido.pdf", tmp_plot, 
       device = cairo_pdf, width = 14.15, height = 6, 
       units = "cm", dpi = 500)
rm(list = c("tmp_plot"))

# Make Sampling Based Output plots for Anaerostipes reference run --------------
output_summary_table <- readRDS("output/anaero-sampling-summary.RDS")

tmp_plot <- output_summary_table %>%
  arrange(., ratio_id) %>%
  mutate(., ratio_id = as.factor(ratio_id)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("Sample Size =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_histogram(aes(y = ..density.., x = type_I_error_rate_conditional, 
                     color = method, fill = method), 
                 alpha =0.8, size = 0.2, position = "identity", 
                 color = "black", breaks = seq(0,0.15, by = 0.005)) +
  geom_density(aes(x = type_I_error_rate_conditional, color = method, 
                   fill = method), 
               alpha=0.2, bw = 0.005, kernel = "gaussian") +
  geom_vline(xintercept = 0.05, color = "black", lty = "dotted", size = 0.5) + 
  xlab(label = "Type-I Error Rate") +
  ylab(label = "Density") +
  labs(color = "Method:", fill = "Method:") +
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_grid(method ~ n_samples, scales = "free_y")
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 8, lab_text_size = 6, 
                            space_legend = 0.5) +
  scale_color_manual(values = c("black","#D55E00", "#009E73", "#56B4E9")) + 
  scale_fill_manual(values = c("black","#D55E00", "#009E73", "#56B4E9"))

ggsave(filename = "figures/sampling-type-1-error-anaero.pdf", tmp_plot, 
       device = cairo_pdf, width = 14.15, height = 6, 
       units = "cm", dpi = 500)
rm(list = c("tmp_plot"))

## Make Parametric Output plots for Bifidobacterium ext. reference run ---------
## Load Data
output_summary_table <- 
  readRDS("output/bifido-parametric-extended-summary.RDS")

## Create Print friendly factor levels 
tmp_levels <- paste(c("Low Dispersion", "Median Dispersion", "High Dispersion"),
                    " (≈", 
                    sort(unique(round(output_summary_table$disp_control,2))),
                    ")", sep = "")
tmp_disp_labels <- tibble(disp_control = 
                            sort(unique(output_summary_table$disp_control)),
                          disp = factor(tmp_levels, levels = tmp_levels))
tmp_levels <- c("1/999", "1/99", "1/9", "1/1")
tmp_odds_labels <- tibble(odds_control = 
                            sort(unique(output_summary_table$odds_control)),
                          odds = factor(tmp_levels, 
                                        levels = tmp_levels))
tmp_levels <- paste("Baseline Odds = ", c("1/999", "1/99", "1/9", "1/1"), 
                    sep = "")
tmp_odds_labels_long <- 
  tibble(odds_control = sort(unique(output_summary_table$odds_control)),
         odds = factor(tmp_levels, levels = tmp_levels))

tmp_effect_labels <- output_summary_table %>%
  mutate(., diff = abs(mu_control - mu_new)) %>%
  mutate(., effect = paste(diff, " (", round(odds_ratio, 2), ")", 
                           sep = "")) %>%
  arrange(., diff) %>%
  mutate(., effect = factor(effect, levels = unique(effect))) %>%
  select(., scenario_id, diff, effect)

## Power - Median Dispersion across all sample sizes & reference scenarios
tmp_plot <- output_summary_table %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels_long, by = "odds_control") %>%
  left_join(., tmp_effect_labels, by = "scenario_id") %>%
  filter(., disp == "Median Dispersion (≈0.46)") %>%
  mutate(., diff = abs(mu_control - mu_new)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_point(aes(y = power_conditional, x = effect, color = method, pch = method), 
             alpha =0.6, size = 1) +
  geom_line(aes(y = power_conditional, x = effect, color = method, group = method), 
            alpha =0.6, size = 0.4, lty = "solid") +
  geom_hline(yintercept = 0.8, color = "black", lty = "dotted", size = 0.3) +
  scale_y_continuous(breaks = c(0.00, 0.20, 0.40, 0.60, 0.80, 1.00)) +
  ylab(label = "Power (1 - Type-II Error Rate)") +
  xlab(label = "Absolute Difference in Success Probability Between Groups (Odds-Ratio)") +
  labs(color = "Method:", pch = "Method:") +
  facet_grid(n_samples ~ odds, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90))
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 9)  
ggsave(filename = "figures/parametric-power-bifido-ext-median-dispersion.pdf", tmp_plot, 
       device = cairo_pdf, width = 14.25, height = 14.25, 
       units = "cm", dpi = 500)
rm(list = c("tmp_plot"))

## Power - All Dispersion levels across all sample sizes & reference scenarios 
tmp_plot <- output_summary_table %>%
  left_join(., tmp_disp_labels, by = "disp_control") %>%
  left_join(., tmp_odds_labels_long, by = "odds_control") %>%
  left_join(., tmp_effect_labels, by = "scenario_id") %>%
  mutate(., diff = abs(mu_control - mu_new)) %>%
  mutate(., n_samples = 
           factor(n_samples, 
                  labels = paste("n =", sort(unique(n_samples))))) %>%
  ggplot(.) +
  geom_line(aes(y = power_conditional, x = effect, color = method, 
                group = interaction(method, disp)), 
            alpha =0.5, size = 0.35) +
  geom_point(aes(y = power_conditional, x = effect, color = method, 
                 group = interaction(method, disp), 
                 pch = disp), 
             alpha =0.8, size = 0.8) +
  geom_hline(yintercept = 0.8, color = "black", lty = "dotted", size = 0.3) +
  ylab(label = "Power (1 - Type-II Error Rate)") +
  xlab(label = "Absolute Difference in Success Probability Between Groups (Odds-Ratio)") +
  labs(color = "Method:", pch = "Dispersion:") +
  facet_grid(n_samples ~ odds, scales = "free_x") +
  scale_y_continuous(breaks = c(0.00, 0.20, 0.40, 0.60, 0.80, 1.00)) +
  theme(axis.text.x = element_text(angle = 90))
tmp_plot <- .plot_theme_mod(tmp_plot, text_size = 9, lab_text_size = 7)  +
  scale_shape_manual(values = c(1,3,4))
ggsave(filename = "figures/parametric-power-bifido-ext-all.pdf", tmp_plot, 
       device = cairo_pdf, width = 18.00, height = 14.25, 
       units = "cm", dpi = 500)
