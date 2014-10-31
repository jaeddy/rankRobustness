library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)

tidy_data <- function(df) {
    df %>%
        gather(key, value, -subtype) %>%
        separate(key, into = c("method", "type", "bc"), sep = "_") %>%
        dcast(subtype + method + bc ~ type) %>%
        group_by(method, bc) %>%
        mutate(acc = (sens + spec)/2) %>%
        melt(id.vars = c("subtype", "method", "bc"), variable.name = "type")
}

plot_data <- function(df) {
    df %>%
        mutate(process = relevel(factor(bc), "nbc"),
               method = relevel(factor(method), "rank")) %>%
        filter(type == "acc") %>%
        ggplot(aes(x = subtype, y = value)) +
        geom_bar(aes(fill = method, alpha = bc), 
                 stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"),
                         name = "Method",
                         breaks = c("gep", "dirac", "rank"),
                         labels = c("GEV", "DIRAC", "Rank")) +
        scale_alpha_discrete(range = c(0.4, 1), guide = "none") +
        ylab("Balanced Accuracy") +
        xlab("Subtype") +
        coord_flip() +
        theme_bw(base_size = 8) +
        theme(legend.position = "top")
}

plot_data_wide <- function(df) {
    df %>%
        mutate(process = relevel(factor(bc), "nbc"),
               method = relevel(factor(method), "gep")) %>%
        filter(type == "acc") %>%
        ggplot(aes(x = subtype, y = value)) +
        geom_bar(aes(fill = method, alpha = process), 
                 stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"),
                          name = "Method",
                          breaks = c("gep", "dirac", "rank"),
                          labels = c("GEV", "DIRAC", "Rank")) +
        scale_alpha_discrete(range = c(0.4, 1), guide = "none") +
        ylab("Balanced Accuracy") +
        xlab("Subtype") +
        theme_bw(base_size = 8) +
        theme(legend.position = "top")
}

get_overall <- function(df) {
    df %>%
        group_by(method, bc, type) %>%
        mutate(overall = mean(value)) %>% 
        dcast(method + bc + type + overall ~ subtype) %>%
        melt(id.vars = c("method", "bc", "type"),
             variable.name = "subtype")
}

format_data <- function(df) {
    df %>%
        dcast(subtype + method + bc ~ type) %>%
        mutate(bc = ifelse(bc == "bc", "yes", "no")) %>%
        select(-acc)
}

# Validation of signatures with pediatric 1 AML dataset
subtypes <- c("MLL", "t(8;21)", "inv(16)", "t(15;17)", "NUP98.NSD1",
              "CEPBA.dm", "NPM1")

dirac_sens_disc <- c(83, 86, 75, 89, 13, 52, 30)
dirac_sens_val <- c(91, 80, 75, 100, 20, 20, 86)
dirac_spec_disc <- c(97, 100, 97, 100, 99, 100, 99)
dirac_spec_val <- c(100, 100, 100, 100, 100, 100, 99)

gep_sens_disc <- c(93, 98, 97, 92, 23, 58, 53)
gep_sens_val <- c(95, 100, 100, 100, 20, 60, 85)
gep_spec_disc <- c(98, 99, 100, 100, 99, 99, 99)
gep_spec_val <- c(97, 100, 98, 100, 100, 100, 98)

rank_sens_disc <- c(94, 98, 95, 95, 40, 91, 54)
rank_sens_val <- c(100, 60, 100, 100, 40, 40, 29)
rank_spec_disc <- c(98, 100, 99, 100, 99, 99, 98)
rank_spec_val <- c(100, 100, 98, 100, 98, 100, 98)


# Validation of signatures with pediatric 2 AML dataset
subtypes_ped <- c("MLL", "t(8;21)", "inv(16)", "t(15;17)")

dirac_sens_nbc <- c(86, 25, 25, 100)
dirac_sens_bc <- c(87, 75, 81, 86)
dirac_spec_nbc <- c(86, 100, 98, 82)
dirac_spec_bc <- c(77, 100, 96, 98)

gep_sens_nbc <- c(93, 100, 100, 100)
gep_sens_bc <- c(93, 100, 94, 100)
gep_spec_nbc <- c(6, 13, 0, 23)
gep_spec_bc <- c(62, 100, 92, 100)

rank_sens_nbc <- c(80, 100, 44, 86)
rank_sens_bc <- c(87, 100, 94, 71)
rank_spec_nbc <- c(95, 100, 100, 80)
rank_spec_bc <- c(67, 100, 100, 100)

ped2AMLstats <- data.frame(subtype = subtypes_ped,
                           dirac_sens_nbc = dirac_sens_nbc,
                           dirac_sens_bc = dirac_sens_bc,
                           dirac_spec_nbc = dirac_spec_nbc,
                           dirac_spec_bc = dirac_spec_bc,
                           gep_sens_nbc = gep_sens_nbc,
                           gep_sens_bc = gep_sens_bc,
                           gep_spec_nbc = gep_spec_nbc,
                           gep_spec_bc = gep_spec_bc,
                           rank_sens_nbc = rank_sens_nbc,
                           rank_sens_bc = rank_sens_bc,
                           rank_spec_nbc = rank_spec_nbc,
                           rank_spec_bc = rank_spec_bc)

ped2AMLplot <- tidy_data(ped2AMLstats)

cairo_ps(filename = "figures/acc_fig_ped.eps", 
         width = 3.35, height = 7, pointsize = 8)

plot_data(get_overall(ped2AMLplot))
ggsave(filename = "figures/acc_fig_ped.pdf",
       width = 8.5, height = 10, units = "cm", dpi = 600)
dev.off()

ped2AMLprint <- format_data(ped2AMLplot)

# Validation of signatures with adult AML dataset
subtypes_ad <- c("MLL", "t(8;21)", "inv(16)", "t(15;17)", "CEPBA.dm", "NPM1")

dirac_sens_nbc <- c(88, 61, 79, 76, 96, 11)
dirac_sens_bc <- c(79, 90, 88, 80, 73, 22)
dirac_spec_nbc <- c(68, 100, 100, 100, 100, 100)
dirac_spec_bc <- c(73, 98, 98, 100, 100, 100)

gep_sens_nbc <- c(100, 100, 100, 100, 100, 100)
gep_sens_bc <- c(97, 100, 95, 100, 69, 34)
gep_spec_nbc <- c(26, 53, 0, 9, 0, 0)
gep_spec_bc <- c(64, 100, 99, 99, 100, 100)

rank_sens_nbc <- c(85, 100, 98, 96, 96, 62)
rank_sens_bc <- c(97, 100, 98, 96, 80, 35)
rank_spec_nbc <- c(93, 100, 98, 100, 96, 100)
rank_spec_bc <- c(66, 100, 99, 99, 100, 100)

adAMLstats <- data.frame(subtype = subtypes_ad,
                           dirac_sens_nbc = dirac_sens_nbc,
                           dirac_sens_bc = dirac_sens_bc,
                           dirac_spec_nbc = dirac_spec_nbc,
                           dirac_spec_bc = dirac_spec_bc,
                           gep_sens_nbc = gep_sens_nbc,
                           gep_sens_bc = gep_sens_bc,
                           gep_spec_nbc = gep_spec_nbc,
                           gep_spec_bc = gep_spec_bc,
                           rank_sens_nbc = rank_sens_nbc,
                           rank_sens_bc = rank_sens_bc,
                           rank_spec_nbc = rank_spec_nbc,
                           rank_spec_bc = rank_spec_bc)

adAMLplot <- tidy_data(adAMLstats)

plot_data(get_overall(adAMLplot))
ggsave(filename = "figures/acc_fig_ad.pdf",
       width = 8.5, height = 13.2, units = "cm", dpi = 600)


