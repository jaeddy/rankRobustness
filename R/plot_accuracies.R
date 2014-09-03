library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

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

ped1AMLstats <- data.frame(subtype = subtypes,
                           dirac_sens_disc = dirac_sens_disc,
                           dirac_sens_val = dirac_sens_val,
                           dirac_spec_disc = dirac_spec_disc,
                           dirac_spec_val = dirac_spec_val,
                           gep_sens_disc = gep_sens_disc,
                           gep_sens_val = gep_sens_val,
                           gep_spec_disc = gep_spec_disc,
                           gep_spec_val = gep_spec_val,
                           rank_sens_disc = rank_sens_disc,
                           rank_sens_val = rank_sens_val,
                           rank_spec_disc = rank_spec_disc,
                           rank_spec_val = rank_spec_val)

ped1AMLplot <- ped1AMLstats %>%
    gather(key, value, -subtype) %>%
    separate(key, into = c("method", "type", "test"), sep = "_") %>%
    dcast(subtype + method + test ~ type) %>%
    group_by(method, test) %>%
    mutate(acc = (sens + spec)/2) %>%
    melt(id.vars = c("subtype", "method", "test"), variable.name = "type")

ggplot(ped1AMLplot, aes(x = subtype, y = value)) +
    geom_bar(aes(fill = method), stat = "identity", position = "dodge") +
    facet_grid( ~ type + test)




