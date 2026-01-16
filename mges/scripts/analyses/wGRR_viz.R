# visualization of the wGRR in urban and natural soils

library(tidyverse)
library(car)
library(lme4)
pal1 = c("#FF823B", "#359951")
pal2 = c("#56282D", "#77966D")

grr <- read_csv("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/plas_bact_wgrr.csv")

md <- read_csv("/Volumes/BunnyBike/landuse_functionaltraits/data/LU_metadata.csv")


# Join and filter out the samples with NA wGRR
md_fil <- md %>%
  left_join(.,grr, by = "sample") %>%
  drop_na() %>%
  mutate(wGRR_fixed = case_when(wGRR > 1 ~ 1, 
                          TRUE ~ wGRR))

mean(md_fil$wGRR_fixed)

# Line with error bar
mrelab_lineplot <- md_fil %>%
  dplyr::select(c("urban.natural", "wGRR_fixed")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(wGRR_fixed),
                   wGRR_fixed = mean(wGRR_fixed))

ggplot(md_fil, aes(x = urban.natural, y = wGRR_fixed)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  wGRR_fixed-sd, ymax =  wGRR_fixed+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("wGRR")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plas_bact_wGRR.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(wGRR ~ urban.natural + (1|site), data = md_fil), type = 3)
#                Chisq Df Pr(>Chisq)    
# (Intercept)   82.3469  1     <2e-16 ***
#   urban.natural  0.3032  1     0.5819 

## Model diagnostics
model <- lmer(wGRR ~ urban.natural + (1|site), data = md_fil)
# Residual plots
library(patchwork)

# Extract model data
model_df <- data.frame(
  Fitted = fitted(model),
  Residuals = residuals(model)
)

# Residual plot
p1 <- ggplot(model_df, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  labs(title = "Residuals vs Fitted", 
       x = "Fitted values", 
       y = "Residuals")

# Q-Q plot
p2 <- ggplot(model_df, aes(sample = Residuals)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_bw() +
  labs(title = "Normal Q-Q Plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

# Combine and save
p1 + p2

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/wgrr_model_diagnostics.pdf", 
       width = 10, height = 4, units = "in")

# Shapiro-Wilk test for normality
shapiro.test(residuals(model))
# W = 0.94699, p-value = 0.4436

# Check model (performance)
check_model(model, check=c("linearity", "homogeneity", "normality"))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/wgrr_checkmodel.pdf", device = "pdf", width = 7, height = 5 , units = "in")


