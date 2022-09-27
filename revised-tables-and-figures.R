library(ggplot2)
library(data.table)
load("cache.rda")

df  = dcast(df_plot, method + year ~variable, value.var = "tot")
df = subset(df, !is.na(agb) & !is.na(awp))
df[, dagb2 := awp-awm]

## New figure 1 ####

# prepare data for the figure
dfig <- subset(df_plot, method == "chave14_h")
# intermediate census years for agb fluxes
dfig[variable != "agb", year := year + 2.5]
dfig[variable == "awp", year := year - 0.05]
dfig[variable == "awm", year := year + 0.05]

dfig$var_group <- as.character(dfig$variable)
dfig$var_group[dfig$variable %in% c("awp", "awm")] <- "flux"

dfig$variable <- factor(dfig$variable, levels = c("agb", "dagb", "awp", "awm"))
levels(dfig$variable) <- toupper(levels(dfig$variable))

dfig$var_group <- factor(dfig$var_group)
levels(dfig$var_group) <- c("'AGB (Mg/ha)'", "Delta*'AGB (Mg/ha/yr)'", "'AWP or AWM (Mg/ha/yr)'")

## panel labels
ann_panels <- data.frame(text = letters[1:3], var_group = factor(levels(dfig$var_group), levels = levels(dfig$var_group)))

## horizontal line
ann_line <- data.frame(line = 0, var_group = factor(levels(dfig$var_group)[2], levels = levels(dfig$var_group)))

ggplot2::ggplot(dfig, ggplot2::aes(x = year, y = tot)) +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = lwr, ymax = upr, color = variable)) +
  ggplot2::labs(x = "", y = "", color = "") +
  ggplot2::geom_text(data = ann_panels, aes(label = text), x = -Inf, y = Inf, hjust = 5, vjust = 0.5, size = 5) +
  ggplot2::geom_hline(data = ann_line, aes(yintercept = line), lty = 2) + ## add horizontal line to delta AGB
  ggplot2::theme_classic() +
  ggplot2::coord_cartesian(clip = "off") +  # allow text to be written outside the plot (panel labels)
  ggplot2::scale_color_discrete(breaks = c("AWP", "AWM"), type = c("black", "black", "cyan4", "coral2")) +
  ggplot2::theme(strip.placement = "outside",
                 strip.background = element_blank(),
                 legend.position = c(0.9, 0.3),
                 plot.title = ggplot2::element_text(hjust = 0.5, vjust = 2, size = 15)) +
  ggplot2::facet_wrap(~var_group, scales = "free_y", strip.position = "left", ncol=1, labeller = label_parsed)

ggsave("figures/new_fig1.png", height = 6, width = 6)


# associated table

# prepare data for the figure
dtab <- subset(df_plot, method == "chave14_h+subs")

dtab[, value := paste0(signif(tot, 3), " (", signif(lwr, 3), "-", signif(upr, 3), ")")]

dtab <- data.table::dcast(dtab, year ~ variable, value.var = "value")

dtab$interval<- c(paste(dtab$year-5, dtab$year, sep = " - ")[-1], NA)

knitr::kable(dtab, col.names = c("Census year", "Census interval", "AGB (Mg/ha)", "AWP (Mg/ha/yr)", "AWM (Mg/ha/yr)"))



# Revised Figure 3 - AGB per group ####

load("cache.rda")

# a - size groups

# change flux year to match middle of census interval
df_size[variable != "agb", year := year + 2.5]

# change variable labels (AGB, AWP) to upper case
levels(df_size$variable) <- c("AGB (Mg/ha)", "AWP (Mg/ha/yr)", "AWM (Mg/ha/yr)")

# revert order of size classes to have bigger trees above
df_size$size <- factor(df_size$size, levels = rev(levels(df_size$size)))

## panel labels
ann_panels <- data.frame(text = letters[1:3], 
                         h = c(4, 3, 4), v = 0.5,
                         variable = factor(levels(df_size$variable), 
                                           levels = levels(df_size$variable)))

fig3a_bis_groups_size <- ggplot(df_size) + 
  geom_col(aes(x = year, y = tot, fill = size)) + 
  ggplot2::geom_text(data = ann_panels, aes(label = text, hjust = h, vjust = v), x = -Inf, y = Inf, size = 5) +
  ggplot2::coord_cartesian(clip = "off") +  # allow text to be written outside the plot (panel labels)
  labs(x = "", y = "", fill = "Diameter\nclass (cm)") +
  scale_y_continuous(expand = c(0, 0)) +  
  scale_fill_brewer(palette = "RdYlBu") +
  ggplot2::theme_classic() +
  ggplot2::theme(strip.placement = "outside", 
                 strip.background = element_blank(), 
                 plot.title = ggplot2::element_text(hjust = 0.5, vjust = 2, size = 15)) +
  ggplot2::facet_wrap(~variable, scales = "free", strip.position = "left")

# b - PFT according to Ruger et al 2020

# change flux year to match middle of census interval
df_pft[variable != "agb", year := year + 2.5]

levels(df_pft$variable) <- c("AGB (Mg/ha)", "AWP (Mg/ha/yr)", "AWM (Mg/ha/yr)")

## panel labels
ann_panels <- data.frame(text = letters[4:6], 
                         h = c(4, 3, 5), v = 0.5,
                         variable = factor(levels(df_pft$variable), 
                                           levels = levels(df_pft$variable)))

# order and rename PFT 
df_pft$PFT[is.na(df_pft$PFT)] <- "NA"
df_pft$PFT <- factor(df_pft$PFT, levels = c("slow", "fast", "LLP", "SLB", "intermediate", "NA"))
levels(df_pft$PFT) <- list("Slow" = "slow", "Fast" = "fast", 
                           "Long-lived pioneers" = "LLP", 
                           "Short-lived breeders" = "SLB", 
                           "Intermediate" = "intermediate", 
                           "No value" = "NA")

fig3b_bis_groups_pft <- ggplot(df_pft) + 
  geom_col(aes(x = year, y = tot, fill = PFT)) + 
  ggplot2::geom_text(data = ann_panels, aes(label = text, hjust = h, vjust = v), x = -Inf, y = Inf, size = 5) +
  ggplot2::coord_cartesian(clip = "off") +  # allow text to be written outside the plot (panel labels)
  labs(x = "", y = "", fill = "PFT \n(Ruger et al., 2020)") +
  scale_y_continuous(expand = c(0, 0)) +  
  scale_fill_manual(values = c("orchid", "gold", "seagreen", "royalblue", "darkorange", "grey")) +
  theme_classic() +
  ggplot2::theme(strip.placement = "outside", 
                 strip.background = element_blank(), 
                 plot.title = ggplot2::element_text(hjust = 0.5, vjust = 2, size = 15)) +
  ggplot2::facet_wrap(~variable, scales = "free", strip.position = "left")

# c - PFT according to Rubio et al 2022

# change flux year to match middle of census interval
df_pft2[variable != "agb", year := year + 2.5]

# change variable labels (AGB, AWP) to upper case
levels(df_pft2$variable) <- c("AGB (Mg/ha)", "AWP (Mg/ha/yr)", "AWM (Mg/ha/yr)")

# order and rename PFT 
df_pft2$FG[is.na(df_pft2$FG)] <- "NA"
df_pft2$FG <- factor(df_pft2$FG, levels = c("1", "2", "3", "NA"))
levels(df_pft2$FG) <- list("Short shade-tolerants" = "1", 
                           "Tall shade-tolerants" = "2", 
                           "Tall pioneers" = "3",
                           "No value" = "NA")
## panel labels
ann_panels <- data.frame(text = letters[7:9], 
                         h = c(4, 3, 7), v = 0.5,
                         variable = factor(levels(df_pft2$variable), 
                                           levels = levels(df_pft2$variable)))

fig3c_bis_groups_pft2 <- ggplot(df_pft2) + 
  geom_col(aes(x = year, y = tot, fill = FG)) + 
  ggplot2::geom_text(data = ann_panels, aes(label = text, hjust = h, vjust = v), x = -Inf, y = Inf, size = 5) +
  ggplot2::coord_cartesian(clip = "off") +  # allow text to be written outside the plot (panel labels)
  labs(x = "", y = "", fill = "FG \n(Rubio & Swenson, 2022)") +
  scale_y_continuous(expand = c(0, 0)) +  
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "grey")) +
  theme_classic() +
  ggplot2::theme(strip.placement = "outside", 
                 strip.background = element_blank(), 
                 plot.title = ggplot2::element_text(hjust = 0.5, vjust = 2, size = 15)) +
  ggplot2::facet_wrap(~variable, scales = "free", strip.position = "left")

fig3a_leg <- ggpubr::as_ggplot(ggpubr::get_legend(fig3a_bis_groups_size))
fig3b_leg <- ggpubr::as_ggplot(ggpubr::get_legend(fig3b_bis_groups_pft))
fig3c_leg <- ggpubr::as_ggplot(ggpubr::get_legend(fig3c_bis_groups_pft2))

ggpubr::ggarrange(fig3a_bis_groups_size + theme(legend.position = "none"), fig3a_leg,
                  fig3b_bis_groups_pft + theme(legend.position = "none"), fig3b_leg,
                  fig3c_bis_groups_pft2 + theme(legend.position = "none"), fig3c_leg,
                  ncol = 2, nrow = 3, widths = c(4,1))

ggsave("figures/fig3_new.png", height = 8, width = 10)


## Revised figure 4 ####

load("cache.rda")
df_crown$agb[df_crown$agb>1000] <- 1000
df_crown$awp[df_crown$awp>30] <- 30
df_crown$awp[df_crown$awp<00] <- 0

## get bci habitats from Harms et al 2001
load("bci_habitat.rda")
# change name of habitat to match publication: first letter in upper case and add spaces
levels(bci_habitat$habitat) <- paste0(toupper(substr(levels(bci_habitat$habitat), 1, 1)), 
                                      substr(levels(bci_habitat$habitat), 2, nchar(levels(bci_habitat$habitat))))
levels(bci_habitat$habitat) <- gsub("Hi_", "High ", levels(bci_habitat$habitat))
levels(bci_habitat$habitat) <- gsub("_", " ", levels(bci_habitat$habitat))
levels(bci_habitat$habitat) <- gsub("Young", "Young forest", levels(bci_habitat$habitat))

## Convert bci_habitat to a raster object
bci_habitat[, habitat_f := as.numeric(as.factor(habitat))]
bci_habitat$secondary = bci_habitat$habitat=="young"
## translate x and y to have center of quadrats
bci_habitat$x = bci_habitat$x+10; bci_habitat$y = bci_habitat$y + 10
r <- raster::rasterFromXYZ(bci_habitat[, c("x", "y", "habitat_f")])

## Extract polygons and make it a simple feature 
habitats <- raster::rasterToPolygons(r, dissolve = TRUE) 
habitats <- sf::st_as_sf(habitats) 
habitats <- merge(habitats, unique(bci_habitat[, c("habitat", "habitat_f")]))
habitats_buff <- sf::st_buffer( habitats, -2)

gagb <- ggplot(df_crown) + 
  geom_raster(aes(x = x, y = y, fill = agb)) +
  geom_sf(data = habitats_buff, aes(color = habitat), size = 1.4, fill = NA) +
  scale_fill_gradient(low = "black", high = "white", 
                      breaks = c(0, 250, 500, 750, 1000),
                      labels = c(0, 250, 500, 750, ">1000")) +
  labs(fill = "AGB\n(Mg/ha)", x = "", y = "") +
  expand_limits(x = c(0, 1000), y = c(0, 500)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_discrete(guide = "none") +
  # scale_color_brewer(palette = "Set2") +
  theme_bw()

gawp <- ggplot(df_crown) + 
  geom_raster(aes(x = x, y = y, fill = awp)) +
  geom_sf(data = habitats_buff, aes(color = habitat), size = 1.4, fill = NA) +
  scale_fill_gradient(low = "black", high = "white",
                      breaks = c(0, 10, 20, 30),
                      labels = c(0, 10, 20, ">30")) +
  labs(fill = "AWP\n(Mg/ha/yr)", x = "", y = "") +
  expand_limits(x = c(0, 1000), y = c(0, 500)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_discrete(guide = "none") +
  # scale_color_brewer(palette = "Set1") +
  theme_bw() 


## stats by habitat
# add habitat information to df_crown
# 1. add quadrat based on x and y
df_crown[, `:=` (q1 = as.character(floor(x/20)), q2 = as.character(floor(y/20)))]
df_crown[nchar(q2) == 1, q2 := paste0("0", q2)]
df_crown[, q20 := as.numeric(paste0(q1, q2))]
df_crown[, `:=`(q1 = NULL, q2 = NULL)]
# merge with bci_habitats
df_crown <- merge(df_crown, bci_habitat[, c("q20", "habitat")], by = "q20")

dfhab2 <- data.table::melt(df_crown, measure.vars = c("agb", "awp", "awm"), id.vars = c("habitat", "q20"))
# group by quadrat
dfhab2 <- dfhab2[, .(value = mean(value)), .(variable, q20, habitat)]

dfhab2 <- dfhab2[, .(mean = mean(value), se = sd(value)/sqrt(length(value))), .(habitat, variable)]
# levels(dfhab2$variable) <- c("AGB (Mg/ha)", "AWP (Mg/ha/yr)")

gagb2 <- ggplot(subset(dfhab2, variable == "agb"), 
                aes(x=habitat, y = mean, ymin = mean-se, ymax = mean+se, col = habitat)) +
  geom_pointrange() +
  labs(x = "", y = "AGB (Mg/ha)") +
  # scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        legend.position = "none")

gawp2 <- ggplot(subset(dfhab2, variable == "awp"), 
                aes(x=habitat, y = mean, ymin = mean-se, ymax = mean+se, col = habitat)) +
  geom_pointrange() +
  labs(x = "", y = "AWP (Mg/ha/yr)") +
  # scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        legend.position = "none")

# habitat map
ghabitat <- ggplot(df_crown) +
  geom_sf(data = habitats, aes(fill = habitat), size = 1.4) + 
  labs(x = "", y = "", fill = "") +
  expand_limits(x = c(0, 1000), y = c(0, 500)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_color_brewer(palette = "Set1") +
  theme_bw() 

gblank <- ggplot() +theme(panel.background = element_blank())

# group graphs
gmap <- ggpubr::ggarrange(gagb, gawp, ghabitat, nrow=3, labels = c("a","b","c"))
gbar <- ggpubr::ggarrange(gagb2, gawp2, gblank, nrow=3, labels = c("d", "e", ""), heights = c(3,3,2))

ggpubr::ggarrange(gmap, gbar, widths = c(2,1), ncol=2)

ggsave("figures/fig4_new.png", height = 10, width = 12)


### New figure 2 - BCI in comparison to other plots ####

if (!dir.exists("MSullivan2020")) {
  # download data from Sullivan et al 2020
  sullivan_url <- "https://forestplots.net/upload/data-packages/sullivan-et-al-2020/MSullivan2020.zip"
  download.file(sullivan_url, destfile = "MSullivan2020.zip", mode = "wb")
  
  # unzip
  utils::unzip("MSullivan2020.zip", exdir = "MSullivan2020")
}

sullivan_mul <- read.csv("MSullivan2020/Data package/MultiCensusPlots.csv")
sullivan_all <- read.csv("MSullivan2020/Data package/AllPlots.csv")

df_mean <- df_plot[method=="chave14_h+subs", .(value = mean(tot)), .(variable)]
df_mean <- data.table::dcast(df_mean, .~variable)

sagb <- ggplot(sullivan_all, aes(x = AGB)) + 
  geom_density(alpha = 0.5, aes(fill = Continent)) +
  geom_vline(col=2, xintercept = df_mean$agb) +
  annotate(geom = "text", col = 2, label = "BCI", x = df_mean$agb*1.25, y = 0, vjust = -1) +
  labs(x = "AGB (Mg/ha)", y = "") +
  expand_limits(x = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.75))

sawp <- ggplot(sullivan_mul, aes(x = AGWP)) + 
  geom_density(alpha = 0.5, aes(fill = Continent)) +
  geom_vline(col=2, xintercept = df_mean$awp) +
  annotate(geom = "text", col = 2, label = "BCI", x = df_mean$awp*0.9, y = 0, vjust = -1) +
  labs(x = "AWP (Mg/ha/yr)", y = "") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "none")

sawm <- ggplot(subset(sullivan_mul, Mort < 22), aes(x = Mort)) + 
  geom_density(alpha = 0.5, aes(fill = Continent)) +
  geom_vline(col=2, xintercept = df_mean$awm) +
  annotate(geom = "text", col = 2, label = "BCI", x = df_mean$awm*0.8, y = 0, vjust = -1) +
  labs(x = "AWM (Mg/ha/yr)", y = "") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "none")

sdagb <- ggplot(subset(sullivan_mul, Mort < 20), aes(x = AGWP-Mort)) + 
  geom_density(alpha = 0.5, aes(fill = Continent)) +
  geom_vline(col=2, xintercept = df_mean$dagb) +
  annotate(geom = "text", col = 2, label = "BCI", x = (df_mean$dagb)+2, y = 0, vjust = -1) +
  labs(x = bquote(Delta*'AGB (Mg/ha/yr)'), y = "") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "none")

ggpubr::ggarrange(sagb, sdagb, sawp,sawm, ncol=2, nrow = 2, labels = "auto")

ggsave("figures/new_fig2.png", height = 6, width = 8)
