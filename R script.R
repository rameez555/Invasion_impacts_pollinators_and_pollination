## R scripts for reproducing the analysis for "Impacts of plant invasion on pollinators and pollination" ###

library(pacman)
pacman::p_load(rgdal,
               ggplot2,
               rnaturalearth,
               dplyr,
               tidyr,
               metafor, # main package for meta-analysis
               broom,
               tidyverse, 
               rmarkdown,
               readxl,
               kableExtra,  # making nice tables
               orchaRd,
               patchwork, # putting ggplots together
               sjPlot,
               gridExtra, 
               jtools,
               ggtext,
               purrr,
               corpcor,
               R.rsp,
               tibble, # nice tables
               ggpubr, # arranging multipe ggplots
               here # making reading files easy
)

source("SAFE_fun.R")  # from Nakagawa et al. 2026 (https://doi.org/10.32942/X24666)

# lnM - delta
##  ln M  (non-SAFE delta-method)  ─────────────────────────────────────
get_lnM_raw <- function(m1, m2, s1, s2, n1, n2) {
  out <- lnM_delta1_indep(m1, m2, s1, s2, n1, n2)   # <- returns c(point, var)
  tibble(
    yi_lnM_raw = out["point"],
    vi_lnM_raw = out["var"]          # ← sampling variance here
  )
}

# lnM - SAFE
##  SAFE ln M  (parametric bootstrap)  ─────────────────────────────────
##        B = 10,000 resamples is a good default.
get_lnM_safe <- function(m1, m2, s1, s2, n1, n2, min_kept   = 1e5, 
                         chunk_init = 5000) {
  out <- safe_lnM_indep(m1, m2, s1, s2, n1, n2)
  tibble(
    yi_lnM_safe = out$point,
    vi_lnM_safe = out$var,           # bootstrap sampling variance
    draws_kept  = out$kept,
    draws_total = out$total
  )
}

# function to modify orchard plots aesthitics

orchard_plot <- function(object, mod = "1", group, data, xlab, N = NULL,
                         alpha = 0.5, angle = 90, cb = TRUE, k = TRUE, g = TRUE,
                         trunk.size = 3, branch.size = 1.2, twig.size = 0.5, errorbar = 0.03, precision.size = 3.5, circle.col = "#999999",
                         transfm = c("none", "tanh", "exp"), condition.lab = "Condition",
                         legend.pos = c("bottom.right", "bottom.left",
                                        "top.right", "top.left",
                                        "top.out", "bottom.out",
                                        "none"), # "none" - no legends
                         k.pos = c("right", "left", "none"),
                         colour = FALSE,
                         fill = TRUE,
                         weights = "prop", by = NULL, at = NULL, upper = TRUE)
{
  ## evaluate choices, if not specified it takes the first choice
  transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  
  if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){
    
    if(mod != "1"){
      results <-  mod_results(object, mod, group, data, N,
                              by = by, at = at, weights = weights, upper = upper)
    } else {
      results <-  mod_results(object, mod = "1", group, data, N,
                              by = by, at = at, weights = weights, upper = upper)
    }
  }
  
  if(any(class(object) %in% c("orchard"))) {
    results <- object
  }
  
  mod_table <- results$mod_table
  
  mod_table$significant <-
    ifelse(mod_table$lowerCL > 0 | mod_table$upperCL < 0,
           "Significant",
           "Not significant")
  
  data_trim <- results$data
  data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)
  
  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
  legend <- "Precision"
  
  if(any(N != "none")){
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
    #latex2exp::TeX()
  }
  
  if(transfm == "tanh"){
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- Zr_to_r(mod_table[,cols])
    data_trim$yi <- Zr_to_r(data_trim$yi)
    label <- xlab
  }else{
    label <- xlab
  }
  
  
  if(transfm == "exp"){
    lnRR_to_perc <- function(df){
      EXP <- function(d) {1-exp(d)}
      return(sapply(df, EXP))
    }
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- lnRR_to_perc(mod_table[,cols])
    data_trim$yi <- lnRR_to_perc(data_trim$yi)
    label <- xlab
  }else{
    label <- xlab
  }
  
  # Add in total effect sizes for each level
  mod_table$K <- as.vector(by(data_trim, data_trim[,"moderator"], function(x) length(x[,"yi"])))
  
  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  mod_table$g <- as.vector(num_studies(data_trim, moderator, stdy)[,2])
  
  # the number of groups in a moderator & data points
  group_no <- length(unique(mod_table[, "name"]))
  
  #data_no <- nrow(data)
  
  # colour blind friendly colours with grey
  cbpl <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  
  # setting fruit colour
  if(colour == TRUE){
    color <- as.factor(data_trim$stdy)
    color2 <- NULL
  }else{
    color <- data_trim$mod
    color2 <- mod_table$name
  }
  
  # whether we fill fruit or not
  if(fill == TRUE){
    fill <- color
  }else{
    fill <- NULL
  }
  
  # whether marginal
  if(names(mod_table)[2] == "condition"){
    
    # the number of levels in the condition
    condition_no <- length(unique(mod_table[, "condition"]))
    
    plot <- ggplot2::ggplot() +
      # pieces of fruit (bee-swarm and bubbles)
      ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21) +
      
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.5) +
      # creating CI
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL), size = branch.size, position = ggplot2::position_dodge2(width = 0.3)) +
      # drowning point estimate and PI
      
      # this will only work for up to 5 different conditions
      # flipping things around (I guess we could do use the same geoms but the below is the original so we should not change)
      ggplot2::scale_shape_manual(values =  20 + (1:condition_no)) + ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1)) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::theme(legend.direction="horizontal") +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
      ggplot2::labs(shape = condition.lab) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
                                                         hjust = 0.5,
                                                         angle = angle))
    
  } else {
    
    plot <- ggplot2::ggplot() +
      # pieces of fruit (bee-swarm and bubbles)
      ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21, col="#999999") + 
      
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.3) +
      # drowning point estimate and PI
      geom_linerange(
        data=mod_table,
        inherit.aes=FALSE,
        aes(
          x=name,
          ymin=lowerPR,
          ymax=upperPR
        ),
        colour="grey60",
        linewidth=twig.size
      ) +
      # creating CI # geom_linerange
      ggplot2::geom_errorbar(
        data=mod_table,
        inherit.aes=FALSE,
        aes(
          x=name,
          ymin=lowerCL,
          ymax=upperCL
        ),
        colour="black",
        linewidth=branch.size,
        width=errorbar
      ) +
      
      # Non-significant
      geom_point(
        data=subset(mod_table, significant=="Not significant"),
        aes(x=name, y=estimate),
        shape=23,
        fill="white",
        colour="black",
        size=trunk.size+0.5
      ) +
      # change the point estimate to white color
      geom_point(
        data=subset(mod_table, significant=="Significant"),
        aes(x=name, y=estimate),
        shape=23,
        fill="black",
        colour="black",
        size=trunk.size+0.5
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::theme(legend.direction="horizontal") +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
                                                         hjust = 0.5,
                                                         angle = angle)) + #+
      #ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
      ggplot2::theme(panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line.x = element_line(colour = "black"),
                     axis.ticks.y = element_blank())
    
  }
  
  # adding legend
  if(legend.pos == "bottom.right"){
    plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
  } else if ( legend.pos == "bottom.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
  } else if ( legend.pos == "top.right") {
    plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
  } else if (legend.pos == "top.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
  } else if (legend.pos == "top.out") {
    plot <- plot + ggplot2::theme(legend.position="top")
  } else if (legend.pos == "bottom.out") {
    plot <- plot + ggplot2::theme(legend.position="bottom")
  } else if (legend.pos == "none") {
    plot <- plot + ggplot2::theme(legend.position="none")
  }
  
  # putting colors in
  if(cb == TRUE){
    plot <- plot +
      ggplot2::scale_fill_manual(values = cbpl) +
      ggplot2::scale_colour_manual(values = cbpl)
  }
  
  # putting k and g in
  if(k == TRUE && g == FALSE && k.pos == "right"){
    plot <- plot +
      ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                        label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = precision.size)
  } else if(k == TRUE && g == FALSE && k.pos == "left") {
    plot <- plot +  ggplot2::annotate('text', y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                      label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "left", size = precision.size)
  } else if (k == TRUE && g == TRUE && k.pos == "right"){
    # get group numbers for moderator
    plot <- plot + #ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
      #label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
      #parse = TRUE, hjust = "right", size = 3.5)  
      ggplot2::annotate('text', y = 3, x = (seq(1, group_no, 1)+0.3),
                        label= paste("italic(N)[obs]==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = precision.size) +
      ggplot2::annotate('text', y = 3, x = (seq(1, group_no, 1)+0.2),
                        label= paste("italic(N)[Studies]==", mod_table$g[1:group_no]), parse = TRUE, hjust = "right", size = precision.size)
    
  } else if (k == TRUE && g == TRUE && k.pos == "left"){
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate('text',  y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                     label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, hjust = "left", size = 3.5)
  }
  
  
  return(plot)
}

# loading data file for analysis
dat <- read_csv("data_all_formatted.csv", show_col_types = FALSE)


## Number of observations and studies

dat %>% group_by(Continent) %>% summarise(`Number of observations` = n(),
                                          `Number of studies` = n_distinct(Paper_ID))

dat %>% group_by(Climate_region) %>% summarise(`Number of observations` = n(),
                                               `Number of studies` = n_distinct(Paper_ID))

dat %>% group_by(Ecosystem) %>% summarise(`Number of observations` = n(),
                                          `Number of studies` = n_distinct(Paper_ID))

dat %>% group_by(Study_design) %>% summarise(`Number of observations` = n(),
                                             `Number of studies` = n_distinct(Paper_ID))

dat %>% group_by(Comparison) %>% summarise(`Number of observations` = n(),
                                           `Number of studies` = n_distinct(Paper_ID))

# Meta-analysis
## Calculating effect sizes 
dat <-  escalc(measure = "SMD", 
               m1i = m1i, 
               m2i = m2i, 
               sd1i = sd1i, 
               sd2i = sd2i, 
               n1i = n1i, 
               n2i = n2i,
               data = dat, append = TRUE,
               var.names = c("yi_d",   "vi_d"))


# lnM
dat <- dat %>% mutate(
  lnM_raw = pmap_dfr(
    list(m1i, m2i, sd1i, sd2i, n1i, n2i),
    get_lnM_raw)
) %>% unnest(lnM_raw)

# SAFE is stochastic so set a seed
set.seed(123)

dat <- dat %>% mutate(
  lnM_safe = pmap_dfr(
    list(m1i, m2i, sd1i, sd2i, n1i, n2i),
    get_lnM_safe)
) %>% unnest(lnM_safe)
## we need to report this online 
summary_stats <- list(
  n_total               = nrow(dat),
  
  ## delta-method lnM
  n_lnM_raw_point_NA    = sum(is.na(dat$yi_lnM_raw)),
  n_lnM_raw_var_NA      = sum(is.na(dat$vi_lnM_raw)),
  
  ## SAFE lnM
  n_lnM_safe_point_NA   = sum(is.na(dat$yi_lnM_safe)),
  n_lnM_safe_var_NA     = sum(is.na(dat$vi_lnM_safe)),
  
  ## rescued comparisons (point estimate only)
  n_rescued_by_SAFE     = sum(is.na(dat$yi_lnM_raw) & !is.na(dat$yi_lnM_safe)),
  
  ## rescued variances
  n_var_rescued_by_SAFE = sum(is.na(dat$vi_lnM_raw) & !is.na(dat$vi_lnM_safe))
)

stopifnot(exists("dat"), is.data.frame(dat))

n_total <- nrow(dat)

summary_tbl <- tibble::tibble(
  Metric = c(
    "Total effect sizes",
    "Delta plug-in lnM: point estimate unavailable (undefined)",
    "Delta plug-in lnM: sampling variance unavailable",
    "SAFE lnM: point estimate unavailable",
    "SAFE lnM: sampling variance unavailable",
    "SAFE rescue: point estimate available when Delta plug-in is undefined",
    "SAFE rescue: sampling variance available when Delta plug-in is undefined"
  ),
  Count = c(
    n_total,
    sum(is.na(dat$yi_lnM_raw)),
    sum(is.na(dat$vi_lnM_raw)),
    sum(is.na(dat$yi_lnM_safe)),
    sum(is.na(dat$vi_lnM_safe)),
    sum(is.na(dat$yi_lnM_raw) & !is.na(dat$yi_lnM_safe)),
    sum(is.na(dat$vi_lnM_raw) & !is.na(dat$vi_lnM_safe))
  )
) %>%
  dplyr::mutate(
    Percent = dplyr::if_else(
      Metric == "Total effect sizes",
      NA_real_,
      100 * Count / n_total
    )
  )
# Headline numbers people care about
n_saved <- sum(!is.na(dat$yi_lnM_safe))
p_saved <- 100 * n_saved / n_total


summary_tbl %>%
  kableExtra::kable(
    digits = 1,
    caption = "Data availability for lnM estimates: Delta plug-in vs SAFE (counts and % of all effect sizes).",
    align = c("l", "r", "r")
  ) %>%
  kableExtra::kable_styling(full_width = FALSE)%>%
  kableExtra::scroll_box(width = "100%", height = "400px")


# INTERCEPT MODELS
full_SMD <- rma.mv(yi_d, vi_d, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat, test   = "t")
full_lnM_safe <- rma.mv(yi_lnM_safe, vi_lnM_safe, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat, test   = "t")

full_SMD_results <- mod_results(full_SMD, mod = "1",  group = "Paper_ID", data = dat)
full_lnM_SAFE_results <- mod_results(full_lnM_safe, mod = "1",  group = "Paper_ID", data = dat)

p1<-orchard_plot(full_SMD_results, 
                 #mod = "1", 
                 alpha = 0.4, angle = 90,
                 xlab = "SMD", 
                 group = "Paper_ID", k = TRUE, g = TRUE,
                 trunk.size = 1.5, branch.size = 0.7, twig.size = 0.5, errorbar = 0.05, precision.size = 3.5, 
                 legend.pos = "bottom.right", fill = FALSE,
                 data = dat) +
  scale_y_continuous(expand = c(0,0)) +
  #theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)) + 
  scale_x_discrete(labels = c("Overall")) +
  #paletteer::scale_fill_paletteer_d("nord::aurora") 
  scale_colour_manual(values = "#67A3C2") +
  scale_y_continuous(limits = c(-3,3))

p2<-orchard_plot(full_lnM_SAFE_results, 
                 #mod = "1", 
                 alpha = 0.4, angle = 90,
                 xlab = "lnM", 
                 group = "Paper_ID", k = TRUE, g = TRUE,
                 trunk.size = 1.5, branch.size = 0.7, twig.size = 0.5, errorbar = 0.05, precision.size = 3.5, 
                 legend.pos = "bottom.right", fill = FALSE,
                 data = dat) +
  scale_y_continuous(expand = c(0,0)) +
  #theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)) + 
  scale_x_discrete(labels = c("")) +
  #paletteer::scale_fill_paletteer_d("nord::aurora") 
  #scale_colour_manual(values = "#D8867D") +
  scale_y_continuous(limits = c(-3,3))

png(filename = "Fig_S3.png", width = 5, height = 4, units = "in", type = "windows", res = 400)
(p1 | p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()


# estimating I2 as measure of heterogeneity in percentage
SMD<-round(i2_ml(full_SMD), 2)
lnM_SAFE<-round(i2_ml(full_lnM_safe), 2)

# Meta-regression models
## Moderator: Response variable classification
class_SMD <- rma.mv(yi_d, vi_d, mods = ~ Response_variable_classification -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat, test = "t")
class_lnM_SAFE <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Response_variable_classification -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat, test = "t")

class_SMD_results <- mod_results(class_SMD, mod = "Response_variable_classification",  group = "Paper_ID", data = dat)
class_lnM_SAFE_results <- mod_results(class_lnM_SAFE, mod = "Response_variable_classification",  group = "Paper_ID", data = dat)


p1<- orchard_plot(class_SMD_results, 
                  #mod = "1", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 0.7, twig.size = 0.5, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "bottom.right",
                  data = dat, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(class_lnM_SAFE_results, 
                  #mod = "1", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 0.7, twig.size = 0.5, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "bottom.right",
                  data = dat, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) + 
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure 2.png", width = 7, height = 5, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## subsetting data: pollinator and pollination
dat_pollinator <- dat %>% filter(Response_variable_classification == "Pollinator")

## Moderator: Pollinator taxonomic order (pollinator dataset)
ord_SMD_pollinator <- rma.mv(yi_d, vi_d, mods = ~ Order -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")
ord_lnM_SAFE_pollinator <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Order -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")

ord_SMD_pollinator_results <- mod_results(ord_SMD_pollinator, mod = "Order",  group = "Paper_ID", data = dat_pollinator)
ord_lnM_SAFE_pollinator_results <- mod_results(ord_lnM_SAFE_pollinator, mod = "Order",  group = "Paper_ID", data = dat_pollinator)

p1<- orchard_plot(ord_SMD_pollinator_results, 
                  mod = "Order", 
                  alpha = 0.9, angle = 45,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE, 
                  trunk.size = 2, branch.size = 1, twig.size = 0.8, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))


p2<- orchard_plot(ord_lnM_SAFE_pollinator_results, 
                  mod = "Order", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 2, branch.size = 1, twig.size = 0.8, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_y_continuous(limits = c(-3,3))


png(filename = "Figure S4.png", width = 8, height = 10, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Invader growth form (pollinator dataset)
inv_SMD_pollinator <- rma.mv(yi_d, vi_d, mods = ~ Invader_Growth_form -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")
inv_lnM_SAFE_pollinator <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Invader_Growth_form -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")

inv_SMD_pollinator_results <- mod_results(inv_SMD_pollinator, mod = "Invader_Growth_form",  group = "Paper_ID", data = dat_pollinator)
inv_lnM_SAFE_pollinator_results <- mod_results(inv_lnM_SAFE_pollinator, mod = "Invader_Growth_form",  group = "Paper_ID", data = dat_pollinator)

p1<- orchard_plot(inv_SMD_pollinator_results, 
                  mod = "Invader_Growth_form", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE, 
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.8, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(inv_lnM_SAFE_pollinator_results, 
                  mod = "Invader_Growth_form", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.8, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) +  
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S5.png", width = 6, height = 5, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: response metric (pollinator dataset)
metric_SMD_pollinator <- rma.mv(yi_d, vi_d, mods = ~ Response_variable_broad_cat -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")
metric_lnM_SAFE_pollinator <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Response_variable_broad_cat -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")

metric_SMD_pollinator_results <- mod_results(metric_SMD_pollinator, mod = "Response_variable_broad_cat",  group = "Paper_ID", data = dat_pollinator)
metric_lnM_SAFE_pollinator_results <- mod_results(metric_lnM_SAFE_pollinator, mod = "Response_variable_broad_cat",  group = "Paper_ID", data = dat_pollinator)

p1<- orchard_plot(metric_SMD_pollinator_results, 
                  mod = "Response_variable_broad_cat", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "bottom.right",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(metric_lnM_SAFE_pollinator_results, 
                  mod = "Response_variable_broad_cat", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "bottom.right",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) + 
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure 3.png", width = 7, height = 6, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Study_design (pollinator dataset)
des_SMD_pollinator <- rma.mv(yi_d, vi_d, mods = ~ Study_design -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")
des_lnM_SAFE_pollinator <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Study_design -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")

des_SMD_pollinator_results <- mod_results(des_SMD_pollinator, mod = "Study_design",  group = "Paper_ID", data = dat_pollinator)
des_lnM_SAFE_pollinator_results <- mod_results(des_lnM_SAFE_pollinator, mod = "Study_design",  group = "Paper_ID", data = dat_pollinator)

p1<- orchard_plot(des_SMD_pollinator_results, 
                  mod = "Study_design", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE, 
                  trunk.size = 1.5, branch.size = 0.9, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(des_lnM_SAFE_pollinator_results, 
                  mod = "Study_design", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 0.9, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S6.png", width = 5.5, height = 4, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Paired comparison (pollinator dataset)
com_SMD_pollinator <- rma.mv(yi_d, vi_d, mods = ~ Comparison -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")
com_lnM_SAFE_pollinator <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Comparison -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")

com_SMD_pollinator_results <- mod_results(com_SMD_pollinator, mod = "Comparison",  group = "Paper_ID", data = dat_pollinator)
com_lnM_SAFE_pollinator_results <- mod_results(com_lnM_SAFE_pollinator, mod = "Comparison",  group = "Paper_ID", data = dat_pollinator)

p1<- orchard_plot(com_SMD_pollinator_results, 
                  mod = "Comparison", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE, 
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.8, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(com_lnM_SAFE_pollinator_results, 
                  mod = "Comparison", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.8, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S7.png", width = 7, height = 8, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Climatic region (pollinator dataset)
clim_SMD_pollinator <- rma.mv(yi_d, vi_d, mods = ~ Climate_region -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")
clim_lnM_SAFE_pollinator <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Climate_region -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollinator, test = "t")

clim_SMD_pollinator_results <- mod_results(clim_SMD_pollinator, mod = "Climate_region",  group = "Paper_ID", data = dat_pollinator)
clim_lnM_SAFE_pollinator_results <- mod_results(clim_lnM_SAFE_pollinator, mod = "Climate_region",  group = "Paper_ID", data = dat_pollinator)

p1<- orchard_plot(clim_SMD_pollinator_results, 
                  mod = "Climate_region", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 0.9, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(clim_lnM_SAFE_pollinator_results, 
                  mod = "Climate_region", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 0.9, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollinator, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S8.png", width = 6, height = 4, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## subsetting data: pollination
dat_pollination <- dat %>% filter(Response_variable_classification == "Pollination")

## Moderator: response variable (pollination dataset)
var_SMD_pollination <- rma.mv(yi_d, vi_d, mods = ~ Response_variable_broad_cat -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")
var_lnM_SAFE_pollination <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Response_variable_broad_cat -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")

var_SMD_pollination_results <- mod_results(var_SMD_pollination, mod = "Response_variable_broad_cat",  group = "Paper_ID", data = dat_pollination)
var_lnM_SAFE_pollination_results <- mod_results(var_lnM_SAFE_pollination, mod = "Response_variable_broad_cat",  group = "Paper_ID", data = dat_pollination)

p1<- orchard_plot(var_SMD_pollination_results, 
                  mod = "Response_variable_broad_cat", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(var_lnM_SAFE_pollination_results, 
                  mod = "Response_variable_broad_cat", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure 3.png", width = 7, height = 7, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Climatic region (pollination dataset)
clim_SMD_pollination <- rma.mv(yi_d, vi_d, mods = ~ Climate_region -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")
clim_lnM_SAFE_pollination <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Climate_region -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")

clim_SMD_pollination_results <- mod_results(clim_SMD_pollination, mod = "Climate_region",  group = "Paper_ID", data = dat_pollination)
clim_lnM_SAFE_pollination_results <- mod_results(clim_lnM_SAFE_pollination, mod = "Climate_region",  group = "Paper_ID", data = dat_pollination)

p1<- orchard_plot(clim_SMD_pollination_results, 
                  mod = "Climate_region", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 0.9, twig.size = 0.6, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(clim_lnM_SAFE_pollination_results, 
                  mod = "Climate_region", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 0.9, twig.size = 0.6, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S13.png", width = 7, height = 4, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Invader growth form (pollination dataset)
inv_SMD_pollination <- rma.mv(yi_d, vi_d, mods = ~ Invader_Growth_form -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")
inv_lnM_SAFE_pollination <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Invader_Growth_form -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")

inv_SMD_pollination_results <- mod_results(inv_SMD_pollination, mod = "Invader_Growth_form",  group = "Paper_ID", data = dat_pollination)
inv_lnM_SAFE_pollination_results <- mod_results(inv_lnM_SAFE_pollination, mod = "Invader_Growth_form",  group = "Paper_ID", data = dat_pollination)

p1<- orchard_plot(inv_SMD_pollination_results, 
                  mod = "Invader_Growth_form", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(inv_lnM_SAFE_pollination_results, 
                  mod = "Invader_Growth_form", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S9.png", width = 7, height = 5, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Response plant growth form (pollinator dataset)
np_SMD_pollination <- rma.mv(yi_d, vi_d, mods = ~ Response_Plant_Growth_form -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")
np_lnM_SAFE_pollination <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Response_Plant_Growth_form -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")

np_SMD_pollination_results <- mod_results(np_SMD_pollination, mod = "Response_Plant_Growth_form",  group = "Paper_ID", data = dat_pollination)
np_lnM_SAFE_pollination_results <- mod_results(np_lnM_SAFE_pollination, mod = "Response_Plant_Growth_form",  group = "Paper_ID", data = dat_pollination)

p1<- orchard_plot(np_SMD_pollination_results, 
                  mod = "Response_Plant_Growth_form", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(np_lnM_SAFE_pollination_results, 
                  mod = "Response_Plant_Growth_form", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S10.png", width = 6, height = 4, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Study design (pollination dataset)
des_SMD_pollination <- rma.mv(yi_d, vi_d, mods = ~ Study_design -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")
des_lnM_SAFE_pollination <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Study_design -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")

des_SMD_pollination_results <- mod_results(des_SMD_pollination, mod = "Study_design",  group = "Paper_ID", data = dat_pollination)
des_lnM_SAFE_pollination_results <- mod_results(des_lnM_SAFE_pollination, mod = "Study_design",  group = "Paper_ID", data = dat_pollination)

p1<- orchard_plot(des_SMD_pollination_results, 
                  mod = "Study_design", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(des_lnM_SAFE_pollination_results, 
                  mod = "Study_design", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S14.png", width = 6, height = 4, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Paired comparison (pollination dataset)
com_SMD_pollination <- rma.mv(yi_d, vi_d, mods = ~ Comparison -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")
com_lnM_SAFE_pollination <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Comparison -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")

com_SMD_pollination_results <- mod_results(com_SMD_pollination, mod = "Comparison",  group = "Paper_ID", data = dat_pollination)
com_lnM_SAFE_pollination_results <- mod_results(com_lnM_SAFE_pollination, mod = "Comparison",  group = "Paper_ID", data = dat_pollination)

p1<- orchard_plot(com_SMD_pollination_results, 
                  mod = "Comparison", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(com_lnM_SAFE_pollination_results, 
                  mod = "Comparison", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S15.png", width = 8, height = 8, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Flower symmetry similarity (pollination dataset)
fss_SMD_pollination <- rma.mv(yi_d, vi_d, mods = ~ Flower_symmetry_similarity -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")
fss_lnM_SAFE_pollination <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Flower_symmetry_similarity -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")

fss_SMD_pollination_results <- mod_results(fss_SMD_pollination, mod = "Flower_symmetry_similarity",  group = "Paper_ID", data = dat_pollination)
fss_lnM_SAFE_pollination_results <- mod_results(fss_lnM_SAFE_pollination, mod = "Flower_symmetry_similarity",  group = "Paper_ID", data = dat_pollination)

p1<- orchard_plot(fss_SMD_pollination_results, 
                  mod = "Flower_symmetry_similarity", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(fss_lnM_SAFE_pollination_results, 
                  mod = "Flower_symmetry_similarity", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S12.png", width = 5, height = 4, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()

## Moderator: Flower colour similarity (pollination dataset)
fcs_SMD_pollination <- rma.mv(yi_d, vi_d, mods = ~ Flower_colour_similarity -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")
fcs_lnM_SAFE_pollination <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Flower_colour_similarity -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_pollination, test = "t")

fcs_SMD_pollination_results <- mod_results(fcs_SMD_pollination, mod = "Flower_colour_similarity",  group = "Paper_ID", data = dat_pollination)
fcs_lnM_SAFE_pollination_results <- mod_results(fcs_lnM_SAFE_pollination, mod = "Flower_colour_similarity",  group = "Paper_ID", data = dat_pollination)

p1<- orchard_plot(fcs_SMD_pollination_results, 
                  mod = "Flower_colour_similarity", 
                  alpha = 0.9, angle = 90,
                  xlab = "SMD", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) + 
  scale_y_continuous(limits = c(-3,3))

p2<- orchard_plot(fcs_lnM_SAFE_pollination_results, 
                  mod = "Flower_colour_similarity", 
                  alpha = 0.9, angle = 90,
                  xlab = "lnM", #transfm = "exp",
                  group = "Paper_ID", k = TRUE, g = TRUE,
                  trunk.size = 1.5, branch.size = 1, twig.size = 0.7, errorbar = 0.05, precision.size = 3.5,
                  legend.pos = "none",
                  data = dat_pollination, fill = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_y_continuous(limits = c(-3,3))

png(filename = "Figure S11.png", width = 6, height = 4, units = "in", type = "windows", res = 400)
(p1| p2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()


# Publication bias
## effective sample size
dat$n0 <- (dat$n1i * dat$n2i) / (dat$n1i + dat$n2i)
dat$nSE <- 1 / sqrt(dat$n0)

## effers regression
eggers_SMD <- rma.mv(yi_d, vi_d, mods = ~ nSE, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat, test = "t")
eggers_lnM_SAFE <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ nSE, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat, test = "t")

eggers_SMD_res <- coef(summary(eggers_SMD))
est  <- eggers_SMD_res["nSE", "estimate"]
lcl  <- eggers_SMD_res["nSE", "ci.lb"]
ucl  <- eggers_SMD_res["nSE", "ci.ub"]
pval <- eggers_SMD_res["nSE", "pval"]

stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*", "")))

lab <- sprintf(
  "Slope = %.2f%s; 95%% CI = [%.2f, %.2f]; p = %.3f",
  est, stars, lcl, ucl, pval
)

p1 <- bubble_plot(
  eggers_SMD, mod = "nSE",
  ylab = "SMD",
  xlab = "sqrt(inverse sample size)",
  group = "Paper_ID", legend.pos = "none", k.pos = "bottom.left"
) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = lab,
    hjust = -0.05,
    vjust = 1.3,
    size = 4,
    #fontface = "bold"
  ) +
  scale_colour_manual(values = rep("grey20", 8)) +
  scale_fill_manual(values = "#CC6677") +
  scale_y_continuous(limits = c(-3,3)) +
  theme(axis.line = element_line(color = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

eggers_lnM_SAFE <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ nSE, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat, test = "t")

eggers_lnM_SAFE_res <- coef(summary(eggers_lnM_SAFE))
est  <- eggers_lnM_SAFE_res["nSE", "estimate"]
lcl  <- eggers_lnM_SAFE_res["nSE", "ci.lb"]
ucl  <- eggers_lnM_SAFE_res["nSE", "ci.ub"]
pval <- eggers_lnM_SAFE_res["nSE", "pval"]

stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*", "")))

lab <- sprintf(
  "Slope = %.2f%s; 95%% CI = [%.2f, %.2f]; p = %.3f",
  est, stars, lcl, ucl, pval
)

p2 <- bubble_plot(
  eggers_lnM_SAFE, mod = "nSE",
  ylab = "lnM",
  xlab = "sqrt(inverse sample size)",
  group = "Paper_ID", legend.pos = "none", k.pos = "bottom.right"
) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = lab,
    hjust = -0.05,
    vjust = 1.3,
    size = 4,
    #fontface = "bold"
  ) +
  scale_colour_manual(values = rep("grey20", 8)) +
  scale_fill_manual(values = "#CC6677") +
  scale_y_continuous(limits = c(-3,3)) +
  theme(axis.line = element_line(color = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())




## time lag bias
time_SMD <- rma.mv(yi_d, vi_d, mods = ~ Publication_Year, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat, test = "t")

time_SMD_res <- coef(summary(time_SMD))
est  <- time_SMD_res["Publication_Year", "estimate"]
lcl  <- time_SMD_res["Publication_Year", "ci.lb"]
ucl  <- time_SMD_res["Publication_Year", "ci.ub"]
pval <- time_SMD_res["Publication_Year", "pval"]

stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*", "")))

lab <- sprintf(
  "Slope = %.3f%s; 95%% CI = [%.3f, %.3f]; p = %.3f",
  est, stars, lcl, ucl, pval
)
p3 <- bubble_plot(
  time_SMD, mod = "Publication_Year",
  ylab = "SMD",
  xlab = "Publication year",
  group = "Paper_ID", legend.pos = "none", k.pos = "bottom.left"
) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = lab,
    hjust = -0.01,
    vjust = 1.3,
    size = 4,
    #fontface = "bold"
  ) +
  scale_colour_manual(values = rep("grey20", 8)) +
  scale_fill_manual(values = "#CC6677") +
  scale_y_continuous(limits = c(-3,3)) +
  theme(axis.line = element_line(color = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

time_lnM_SAFE <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods = ~ Publication_Year, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat, test = "t")

time_lnM_SAFE_res <- coef(summary(time_lnM_SAFE))
est  <- time_lnM_SAFE_res["Publication_Year", "estimate"]
lcl  <- time_lnM_SAFE_res["Publication_Year", "ci.lb"]
ucl  <- time_lnM_SAFE_res["Publication_Year", "ci.ub"]
pval <- time_lnM_SAFE_res["Publication_Year", "pval"]

stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*", "")))

lab <- sprintf(
  "Slope = %.3f%s; 95%% CI = [%.3f, %.3f]; p = %.3f",
  est, stars, lcl, ucl, pval
)

p4 <- bubble_plot(
  time_lnM_SAFE, mod = "Publication_Year",
  ylab = "lnM",
  xlab = "Publication year",
  group = "Paper_ID", legend.pos = "none", k.pos = "bottom.left"
) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = lab,
    hjust = -0.01,
    vjust = 1.3,
    size = 4,
    #fontface = "bold"
  ) +
  scale_colour_manual(values = rep("grey20", 8)) +
  scale_fill_manual(values = "#CC6677") +
  scale_y_continuous(limits = c(-3,3)) +
  theme(axis.line = element_line(color = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

png(filename = "Figure S2.png", width = 9, height = 8, units = "in", type = "windows", res = 400)
(p1| p2| p3| p4) + plot_layout(ncol = 2, nrow = 2) +plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
dev.off()


# Sensitivity analysis
 
## variance analogue
dat$n0 <- (dat$n1i * dat$n2i) / (dat$n1i + dat$n2i)
dat$nV <- 1 / dat$n0
full_lnM_safe_sen <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods   = ~ nV, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat, test   = "t")
summary(full_lnM_safe_sen)

## studies with total sample size n1 + n2 ≥ 40
dat_sens <- dat %>% dplyr::filter((n1i + n2i) >= 40)
dim(dat_sens)

### intercept model
full_lnM_safe_sen <- rma.mv(yi_lnM_safe, vi_lnM_safe, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_sens, test   = "t")
summary(full_lnM_safe_sen)
### for pollinators and pollinators
full_lnM_safe_res_sen <- rma.mv(yi_lnM_safe, vi_lnM_safe, mods   = ~ Response_variable_classification -1, random = list(~1 | Paper_ID, ~1 | Observation_ID, ~1 | Invsive_plant_ID, ~1 | Cohort_ID), method = "REML", data = dat_sens, test   = "t")
full_lnM_safe_res_sen_res <- mod_results(full_lnM_safe_res_sen, mod = "Response_variable_classification",  group = "Paper_ID", data = dat_sens)

# session info
library(pander)
sessionInfo() %>% pander()