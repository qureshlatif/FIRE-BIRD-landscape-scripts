library(dplyr)
library(ggplot2)
library(cowplot)

#_____________________________ Inputs ________________________________________#
setwd("F:/research stuff/FS_PostDoc/landscapeDS")
dat <- read.csv("T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/HabByYear_DCFonly.csv",
                header = T, stringsAsFactors = F) %>% tbl_df %>%
  slice(5:n()) %>%
  mutate(HAB_MODHI = HAB_MOD + HAB_HIGH)
out_file <- "Malheur_TimeSeries_DCF.tiff"
#_____________________________________________________________________________#

# Potential nest abundance #
p.hab <- ggplot(dat, aes(x = YEAR, y = HAB_MODHI)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = mean(dat$HAB_MODHI), linetype="dotted", size=1) +
  ylab("Habitat acres (Mod-High)") + xlab("Year") +
  scale_x_continuous(breaks = seq(1990, 2015, by = 5)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.title.y=element_text(size=25))

# Potential nest abundance #
p.nest <- ggplot(dat, aes(x = YEAR, y = PNest_mn)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = PNest_95lo, ymax = PNest_95hi), size=1, width=0.2) +
  ylab("Potential number of nests") + xlab("Year") +
  scale_x_continuous(breaks = seq(1990, 2015, by = 5)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.title.y=element_text(size=25))

theme_set(theme_bw())
p <- ggdraw() + 
  draw_plot(p.hab, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(p.nest, x = 0, y = 0, width = 1, height = 0.5)

save_plot(out_file, p, ncol = 1.5, nrow = 3, dpi = 200)
