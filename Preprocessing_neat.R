#Cleaning Data in order to select samples for sequencing 
#Tuesday 19th May
library(readr)
library(dplyr)
#importing in the metadata file and telling R to ignore the first row so that i can make the first row into column names 
WickedSchisto <- read.csv("Desktop/Project/WickedSchisto.csv", header = FALSE)
View(WickedSchisto)
#removing the first row with the leftover entry
clean_WickedSchisto <- WickedSchisto[-1, ]
#checking
head(clean_WickedSchisto)
#writing out this cleaned file to reclean in excel (avoiding having to transpose the data)
write.csv(clean_WickedSchisto, "clean_WickedSchisto.csv", row.names = FALSE, col.names = FALSE)
#reloading the clean file without the default column and row headers given by R
final_WS <- read.csv("Desktop/Project/clean_WickedSchisto.csv", header = TRUE)
final_WS <- final_WS[, -ncol(final_WS)]
head(final_WS)
All_CAL <- read.csv("Desktop/Project/All_CAL.csv", header = TRUE)
final_All <- left_join(All_CAL, final_WS, by = "ID")
head(final_All)
clean_WS <- final_All %>%
  filter(
    if_all(everything() & -ID & -CAL..ug.g.feces., ~ .x != 999 & .x != "999" & .x != 99 & .x != "99"), 
    if_all(everything() & -ID & -Age & -Height & -Weight & -CAL..ug.g.feces., ~ .x != 99 & .x != "99")
    )
head(clean_WS)


#creating a dataframe with all the malaria positive samples
Malaria <- read.csv("Desktop/Project/Malaria_positive.csv", header=TRUE)
Malaria$ID <- as.numeric(Malaria$ID)
head(Malaria)
Uninfected <- read.csv("Desktop/Project/Malaria_negative.csv", header=TRUE)
Uninfected$ID <- as.numeric(Uninfected$ID)
head(Uninfected)


Malaria_pos_data <- left_join(Malaria, final_WS, by = "ID")
head(Malaria_pos_data)
clean_Malaria_data <- Malaria_pos_data %>%
  filter(
    if_all(everything() & -ID & -CAL..ug.g.feces., ~ .x != 999 & .x != "999" & .x != 99 & .x != "99"), 
    if_all(everything() & -ID & -Age & -Height & -Weight & -CAL..ug.g.feces., ~ .x != 99 & .x != "99")
  )
head(clean_Malaria_data)

Uninfected_data1 <- left_join(Uninfected, final_WS, by = "ID")
head(Uninfected_data1)
Uninfected_data1$CAL..ug.g.feces. <- as.numeric(Uninfected_data1$CAL..ug.g.feces.)
clean_Uninfected_data <- Uninfected_data1 %>%
  filter(
    if_all(everything() & -ID & -CAL..ug.g.feces., ~ .x != 999 & .x != "999" & .x != 99 & .x != "99"), 
    if_all(everything() & -ID & -Age & -Height & -Weight & -CAL..ug.g.feces., ~ .x != 99 & .x != "99")
  )
head(clean_Uninfected_data)

#finding the median and the extremes of uninfected samples based on Calprotectin levels 
median(clean_Uninfected_data$CAL..ug.g.feces., na.rm = TRUE)
min(clean_Uninfected_data$CAL..ug.g.feces., na.rm = TRUE)
max(clean_Uninfected_data$CAL..ug.g.feces., na.rm = TRUE)
mean(clean_Uninfected_data$CAL..ug.g.feces., na.rm = TRUE)

#finding the median and the extremes of infected samples based on Calprotectin levels 
clean_Malaria_data$CAL..ug.g.feces. <- as.numeric(clean_Malaria_data$CAL..ug.g.feces.)
median(clean_Malaria_data$CAL..ug.g.feces., na.rm = TRUE)
min(clean_Malaria_data$CAL..ug.g.feces., na.rm = TRUE)
max(clean_Malaria_data$CAL..ug.g.feces., na.rm = TRUE)
mean(clean_Malaria_data$CAL..ug.g.feces., na.rm = TRUE)

#finding the samples wiht the highest and lowest calprotectin levels in the infected samples 
Malaria_top_10 <- clean_Malaria_data %>%
  filter(CAL..ug.g.feces. >= quantile(CAL..ug.g.feces., 0.90, na.rm =TRUE))
Malaria_top_10
Malaria_low_10 <- clean_Malaria_data %>%
  filter(CAL..ug.g.feces. <= quantile(CAL..ug.g.feces., 0.10, na.rm =TRUE))   
Malaria_low_10

#finding the samples with the highest and lowest calprotectin levels in the group of uninfected
Uninfected_top_10 <- clean_Uninfected_data %>%
  filter(CAL..ug.g.feces. >= quantile(CAL..ug.g.feces., 0.90, na.rm =TRUE))
Uninfected_top_10  
Uninfected_low_10 <- clean_Uninfected_data %>%
  filter(CAL..ug.g.feces. <= quantile(CAL..ug.g.feces., 0.10, na.rm =TRUE))
Uninfected_low_10


clean_Malaria_data$CAL..ug.g.feces. <- as.numeric(clean_Malaria_data$CAL..ug.g.feces.)
cor.test(clean_Malaria_data$BRISTOL, clean_Malaria_data$CAL..ug.g.feces., method = "spearman", use = "complete.obs")
library(ggpubr)
ggboxplot(clean_Malaria_data, 
          x = "BRISTOL", 
          y = "CAL..ug.g.feces.",
          color = "BRISTOL",
          add = "jitter",
          legend = "none",
          outlier.shape = NA) +
  stat_cor(aes(x = as.numeric(BRISTOL), y = CAL..ug.g.feces.),
           method = "spearman",
           cor.coef.name = "rho") +
  scale_y_log10() +
  labs(
    title = "Fecal Calprotectin concentration by Bristol Stool Score for Malaria infected individuals",
    subtitle = "Spearman rank correlation",
    x = "Bristol Stool Scale",
    y = "Calprotectin Level (µg/g of feces, log10 scale)"
  )

ggboxplot(clean_WS, 
          x = "BRISTOL", 
          y = "CAL..ug.g.feces.",
          color = "BRISTOL",
          add = "jitter",
          legend = "none",
          outlier.shape = NA) +
  stat_cor(aes(x = as.numeric(BRISTOL), y = CAL..ug.g.feces.),
           method = "spearman",
           cor.coef.name = "rho") +
  scale_y_log10() +
  labs(
    title = "Fecal Calprotectin concentration by Bristol Stool Score for all samples",
    subtitle = "Spearman rank correlation",
    x = "Bristol Stool Scale",
    y = "Calprotectin Level (µg/g of feces, log10 scale)"
  )

ggboxplot(clean_Uninfected_data, 
          x = "BRISTOL", 
          y = "CAL..ug.g.feces.",
          color = "BRISTOL",
          legend = "none",
          outlier.shape = NA) +
  stat_cor(aes(x = as.numeric(BRISTOL), y = CAL..ug.g.feces.),
           method = "spearman",
           cor.coef.name = "rho",
           label.x = 4.5,         # Moves it horizontally across the screen
           label.y = log10(600)) +
  geom_beeswarm() +
  scale_y_log10() +
  labs(
    title = "Fecal Calprotectin concentration by Bristol Stool Score for uninfected samples",
    subtitle = "Spearman rank correlation",
    x = "Bristol Stool Scale",
    y = "Calprotectin Level (µg/g of feces, log10 scale)"
  )

selected_infected_samples <- c("31303", "62002", "43802", "60503", "34002", "61702", "62301", "50101", "32602", "53601", "40901", "50901", "30801", "51301", "63801", "61003", "60102", "52501", "40303", "51102", "51003", "53301", "51502", "31702", "52702", "61503", "53103", "53401", "51702", "62102", "33301", "63403", "31203", "51803", "52801", "51701")
selected_uninfected <- c("33101", "63601", "32901", "53701", "61201", "40202", "30301", "53001", "61101", "61501", "50301", "50801", "62101", "63201", "63003", "41202", "63301", "43401", "32704", "42501", "51902", "41201", "40201", "44102", "30602", "30401", "43101", "33102", "33603", "62202")

library(dplyr)


malaria_pos_labels <- clean_Malaria_data  %>%
  select(ID) %>% 
  mutate(MALARIA = "Positive")


malaria_neg_labels <- clean_Uninfected_data %>%
  select(ID) %>% 
  mutate(MALARIA = "Negative")  

malaria_master_key <- bind_rows(malaria_pos_labels, malaria_neg_labels)

#MALARIA ONLY CALPROTECTIN ############################
file_with_malaria <- clean_WS %>%
  left_join(malaria_master_key, by = "ID")
file_with_malaria %>%
  filter(!is.na(MALARIA), MALARIA != "") %>%
  ggboxplot( 
          x = "MALARIA", 
          y = "CAL..ug.g.feces.",
          color = "MALARIA",
          legend = "none",
          fill = "white") +
  geom_beeswarm(aes(color = MALARIA)) +
  scale_color_manual(values = c("Positive" = "salmon", "Negative" = "pink")) +
  scale_y_log10() +
  stat_compare_means(
    method = "wilcox.test"
  ) +
  theme_bw() +
  labs(
    title = "Calprotectin Levels (µg per g feaces) for infected and uninfected samples",
    x = "Malaria Status",
    y = "Calprotectin Level (µg/g of feces, log10 scale)",
    caption = "Wilcoxon, p = 0.79"
  )

library(dplyr)
head(selected_infected_samples)

plot_ready_data <- file_with_malaria %>%
  mutate(Infection_Group = case_when(
    MALARIA == "Positive" & (SMA1A > 0 | SMA1B > 0)  ~ "Malaria+ / Schisto+",
    MALARIA == "Positive" & (SMA1A == 0 & SMA1B == 0) ~ "Malaria+ / Schisto-",
    MALARIA == "Negative" & (SMA1A > 0 | SMA1B > 0)  ~ "Malaria- / Schisto+",
    MALARIA == "Negative" & (SMA1A == 0 & SMA1B == 0) ~ "Malaria- / Schisto-",
    TRUE ~ NA_character_ 
  )) %>%
  filter(!is.na(Infection_Group))
mutate(Point_Style = case_when(
  ID %in% selected_infected_samples ~ "Selected1",
  ID %in% selected_uninfected        ~ "Selected2",
  TRUE                              ~ Infection_Group
))

head(plot_ready_data)
ANOVA_model <- aov(log10(CAL..ug.g.feces.) ~ Infection_Group, data = plot_ready_data)
summary(ANOVA_model)
TukeyHSD(ANOVA_model)

library(car)
qqPlot(ANOVA_model$residuals,
       id = FALSE # id = FALSE to remove point identification
)

head(plot_ready_data)

ggplot(plot_ready_data, aes(x = Infection_Group, y = CAL..ug.g.feces.)) +
  geom_boxplot(aes(color = Infection_Group), fill = "white", outlier.shape = NA) +
  geom_beeswarm(aes(color = Infection_Group)) +
  scale_color_manual(values = c(
    "Selected Malaria Samples" = "blue",
    "Selected Uninfected Samples" = "red",
    "Malaria+ / Schisto+" = "lightgreen",
    "Malaria+ / Schisto-" = "lightblue",
    "Malaria- / Schisto+" = "salmon",
    "Malaria- / Schisto-" = "pink"
  )) +
  theme_bw() +
  scale_y_log10() +
  labs(
    title = "Fecal Calprotectin Levels Across Infection Cohorts",
    x = "Infection Status Group",
    y = "Faecal Calprotectin (µg/g of feces, log10 scale)",
    caption = "One-way ANOVA, p = 0.009")

#plotting only uninfected malaria to look at whether bristol stool is correlated with schisto diagnosis
plot_ready3 <- plot_ready_data %>%
  filter(!if_all(
    c(SMA1A, SMA1B), 
    ~ .x == 0 | is.na(.x)
  ))




head(plot_ready3)
ggboxplot(plot_ready3, 
          x = "BRISTOL", 
          y = "CAL..ug.g.feces.",
          color = "BRISTOL",
          legend = "none",
          outlier.shape = NA) +
  stat_cor(aes(x = as.numeric(BRISTOL), y = CAL..ug.g.feces.),
           method = "spearman",
           cor.coef.name = "rho",
           label.x = 4.5,         # Moves it horizontally across the screen
           label.y = log10(600)) +
  geom_beeswarm() +
  scale_y_log10() +
  labs(
    title = "Fecal Calprotectin concentration by Bristol Stool Score for Schisto samples",
    subtitle = "Spearman rank correlation",
    x = "Bristol Stool Scale",
    y = "Calprotectin Level (µg/g of feces, log10 scale)"
  )

library(ggbeeswarm)
plot_uninfected <- plot_ready_data %>%
  filter(MALARIA == "Negative")

ggplot(plot_uninfected, aes(x = Infection_Group, y = CAL..ug.g.feces.)) +
  geom_boxplot(aes(color = Infection_Group), fill = "white", outlier.shape = NA) +
  geom_beeswarm(aes(color = Infection_Group)) +
  scale_color_manual(name = "Schistomaniasis status", values = c(
    "Malaria- / Schisto+" = "salmon",
    "Malaria- / Schisto-" = "pink"
  ),
  labels = c(
    "Malaria- / Schisto+" = "Schisto Infected",
    "Malaria- / Schisto-" = "Uninfected Control"
  )) +
  theme_bw() +
  scale_y_log10() +
  labs(
    x = "Schistomaniasis Status",
    y = "Calprotectin Level (µg/g of feces, log10 scale)",
    caption = "Wilcoxon rank-sum test, p = 0.0049") +
  scale_x_discrete(
    labels = c(
      "Malaria- / Schisto+" = "Schisto Infected",
      "Malaria- / Schisto-" = "Uninfected Control"
    )
  )


       
qPCR <- read.csv("Desktop/Research Project/qPCR.csv", header = TRUE)
View(qPCR)
qPCR$ID <- as.numeric(qPCR$ID)
qPCR$CAL..ug.g.feces. <- as.numeric(qPCR$CAL..ug.g.feces.)
final_qPCR <- left_join(qPCR, Malaria_pos_data, by = "ID")
View(final_qPCR)
final_clean_qPCR <- final_qPCR %>%
  filter(
    if_all(everything() & -ID & -CAL..ug.g.feces., ~ .x != 999 & .x != "999" & .x != 99 & .x != "99"), 
    if_all(everything() & -ID & -Age & -Height & -Weight & -CAL..ug.g.feces., ~ .x != 99 & .x != "99"),
    if_all(everything(), ~ .x != "#N/A")
  )
final_clean_qPCR$Average.CT.values <- as.numeric(final_clean_qPCR$Average.CT.values)
final_clean_qPCR$CAL..ug.g.feces. <- as.numeric(final_clean_qPCR$CAL..ug.g.feces.)
View(final_clean_qPCR) 
ggplot(final_clean_qPCR, aes(x=Average.CT.values, y=CAL..ug.g.feces.)) + 
  geom_point() +
  scale_y_log10()

#Monday 25th May 2026
#filtering out all of the samples that have enteropathogens from the dataset
plot_ready2 <- plot_ready_data %>%
  filter(
  SMA1A == 0 & SMA1B == 0 & ASC1A == 0 & ASC1B == 0 & HW1A == 0 & HW1B == 0 & TT1A == 0 & TT1B == 0 & OTHERS == 0 & SM2A == 0 & SM2B == 0 & ASC2A == 0 & ASC2B == 0 & HW2A == 0 & HW2B == 0 & TT2A == 0 & TT2B == 0 & OTHER == 0)
View(plot_ready2)

plot_ready3 <- plot_ready_data %>%
  filter_out(
    SMA1A == 0 & SMA1B == 0 & ASC1A == 0 & ASC1B == 0 & HW1A == 0 & HW1B == 0 & TT1A == 0 & TT1B == 0 & OTHERS == 0 & SM2A == 0 & SM2B == 0 & ASC2A == 0 & ASC2B == 0 & HW2A == 0 & HW2B == 0 & TT2A == 0 & TT2B == 0 & OTHER == 0)

#separating this new dataset into malaria infected and uninfected 
Malaria_filtered <- plot_ready2 %>%
  filter(
    MALARIA == "Positive")
View(Malaria_filtered)

Uninfected_filtered <- plot_ready2 %>%
  filter(MALARIA == "Negative")
View(Uninfected_filtered)

#distribution of malaria by village - including uninfected individuals - all enteropathogen negative
ggplot(plot_ready2, aes(fill = MALARIA, x= COMMUNITY)) +
  geom_bar(position = "dodge") +
  scale_y_continuous(breaks = seq(0, 40, by = 5)) +
  xlab("Village") +
  ylab("No. of Samples") +
  theme_bw() +
  labs(fill = "Malaria State",
       title = "Malaria status of Enteropathogen-negative individuals across Ugandan Villages")

#boxplot showing the distribution of calprotectin levels of enteropathogen negative but malaria positive individuals across the villages sampled
library(hrbrthemes)
library(viridis)
Malaria_filtered %>%
  ggplot(aes(x = COMMUNITY, y = CAL..ug.g.feces., fill = COMMUNITY)) +
    geom_boxplot(aes(color = COMMUNITY), outlier.shape = NA) +
    geom_jitter(aes(color = COMMUNITY), size=0.8, alpha=0.9) +
    scale_color_manual(values = c(
      "NTOKOL" = "springgreen4",
      "WALUJJO" = "steelblue",
      "WAMBETTE" = "hotpink3"
    )) +
  scale_fill_manual(values = c(
    "NTOKOL" = "lightgreen",
    "WALUJJO" = "lightblue",
    "WAMBETTE" = "pink"
  )) +
  theme_bw() +
    scale_y_log10() +
    xlab("Village Community") +
    ylab("Calprotectin Levels (log10)") +
    labs(fill = "Village Community", color = "Village Community") +
    ggtitle("Fecal Calprotectin levels in Enteropathogen-negative individuals with Malaria across Ugandan Villages")

plot_ready2 %>%
  ggplot(aes(x = COMMUNITY, y = CAL..ug.g.feces., fill = MALARIA)) +
  geom_boxplot(aes(color = MALARIA), outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(), aes(color = MALARIA), size=0.8) +
  scale_color_manual(values = c(
    "Positive" = "steelblue",
    "Negative" = "hotpink3"
  )) +
  scale_fill_manual(values = c(
    "Positive" = "lightblue",
    "Negative" = "pink"
  )) +
  theme_bw() +
  scale_y_log10() +
  xlab("Village Community") +
  ylab("Calprotectin Levels (log10)") +
  labs(fill = "Malaria Status", color = "Malaria Status") +
  ggtitle("Fecal Calprotectin levels in Enteropathogen-negative individuals across Ugandan Villages")


#separating the data based on community:

#getting ntokol data (malaria positive and uninfected)
Ntokol_all_data <- plot_ready2 %>%
  filter(
    COMMUNITY == "NTOKOL")
View(Ntokol_all_data)


ggboxplot(Ntokol_all_data, 
          x = "BRISTOL", 
          y = "CAL..ug.g.feces.",
          fill = "MALARIA",
          add = "jitter",
          outlier.shape = NA) +
  stat_cor(aes(x = as.numeric(BRISTOL), y = CAL..ug.g.feces.),
           method = "spearman",
           cor.coef.name = "rho",
           label.x = 4.5, 
           label.y = log10(600)) +
  scale_y_log10() +
  labs(fill = "Malaria Status",
       title = "Fecal Calprotectin concentration by Bristol Stool Score for Enteropathogen-negative samples within Ntokol Village",
       x = "Bristol Stool Scale",
       y = "Calprotectin Level (µg/g of feces, log10 scale)",
       subtitle = "Ntokol Community") +
  theme(panel.grid.major = element_line(linetype = "dashed", colour = "gray90"), panel.grid.minor = element_line(linetype = "dashed", colour = "gray90"))
head(Ntokol_all_data)

ggboxplot(Ntokol_all_data, 
  x = "MALARIA", 
  y = "CAL..ug.g.feces.",
  color = "MALARIA",
  legend = "none",
  fill = "white") +
  geom_beeswarm(aes(color = MALARIA)) +
  scale_color_manual(values = c("Positive" = "salmon", "Negative" = "pink")) +
  scale_y_log10() +
  theme_bw() +
  labs(
    x = "Malaria Status",
    y = "Calprotectin Level (µg/g of feces, log10 scale)",
    caption = "Wilcoxon, p = 0.34"
  ) +
  theme(panel.grid.major = element_line(linetype = "dashed", colour = "gray90"), panel.grid.minor = element_line(linetype = "dashed", colour = "gray90"))




Walujjo_all_data <- plot_ready2 %>%
  filter(
    COMMUNITY == "WALUJJO")

ggboxplot(Walujjo_all_data, 
          x = "BRISTOL", 
          y = "CAL..ug.g.feces.",
          fill = "MALARIA",
          add = "jitter",
          outlier.shape = NA) +
  stat_cor(aes(x = as.numeric(BRISTOL), y = CAL..ug.g.feces.),
           method = "spearman",
           cor.coef.name = "rho",
           label.x = 4.5, 
           label.y = log10(600)) +
  scale_y_log10() +
  labs(fill = "Malaria Status",
       x = "Bristol Stool Scale",
       y = "Calprotectin Level (µg/g of feces, log10 scale)", 
       title = "Fecal Calprotectin concentration by Bristol Stool Score for Enteropathogen-negative samples within Walujjo Village",
       subtitle = "Walujjo Community")

Wambette_all_data <- plot_ready2 %>%
  filter(
    COMMUNITY == "WAMBETTE")

ggboxplot(Wambette_all_data, 
          x = "BRISTOL", 
          y = "CAL..ug.g.feces.",
          fill = "MALARIA",
          add = "jitter",
          outlier.shape = NA) +
  stat_cor(aes(x = as.numeric(BRISTOL), y = CAL..ug.g.feces.),
           method = "spearman",
           cor.coef.name = "rho",
           label.x = 4.5, 
           label.y = log10(600)) +
  scale_y_log10() +
  labs(fill = "Malaria Status",
       x = "Bristol Stool Scale",
       y = "Calprotectin Level (µg/g of feces, log10 scale)",
       title = "Fecal Calprotectin concentration by Bristol Stool Score for Enteropathogen-negative samples within Wambette Village",
       subtitle = "Wambette Community")

#one big view of the three boxplot graphs produced above  
ggboxplot(plot_ready2, 
          x = "BRISTOL", 
          y = "CAL..ug.g.feces.",
          fill = "MALARIA",
          add = "jitter",
          outlier.shape = NA) +
  stat_cor(aes(x = as.numeric(BRISTOL), y = CAL..ug.g.feces.),
           method = "spearman",
           cor.coef.name = "rho") +
  scale_y_log10() +
  labs(fill = "Malaria Status",
       x = "Bristol Stool Scale",
       y = "Calprotectin Level (µg/g of feces, log10 scale)",
       title = "Fecal Calprotectin concentration by Bristol Stool Score for Enteropathogen-negative samples across villages") +
facet_wrap(~COMMUNITY, ncol = 3)

#the same data but just in bargraph form (this one is plus the uninfected):
ggplot(plot_ready2, aes(fill = MALARIA, x= BRISTOL)) +
  geom_bar(position=position_dodge2(preserve='single')) +
  scale_x_discrete(limits = as.character(1:7), drop = FALSE) +
  scale_y_discrete(limits = as.character(1:12), drop = FALSE) +
  xlab("Bristol Stool Scale") +
  ylab("No. of Samples") +
  labs(fill = "Malaria State") +
  ggtitle("Distribution of diarrhea (as measured with the Bristol Stool Chart) for Enteropathogen-negative individuals across Ugandan Villages") +
  facet_wrap(~COMMUNITY, ncol = 3)

write.csv(plot_ready2, "Desktop/Research Project/filtered_master_data.csv", row.names = FALSE, col.names = FALSE)


#matching controls:
library(dplyr)
library(MatchIt)

Ntokol_all_data <- Ntokol_all_data %>%
  mutate(
    Age = as.numeric(Age),
    Sex = as.factor(SEX)
  )

Ntokol_match_data <- Ntokol_all_data %>%
  mutate(MALARIA = ifelse(MALARIA == "Positive", 1, 0))

match_model <- matchit(MALARIA ~ Age + CAL..ug.g.feces.,
                       data = Ntokol_match_data,
                       method = "nearest",
                       distance = "mahalanobis",
                       ratio = 1)
ntokol_matched <- match.data(match_model)
View(ntokol_matched)
summary(match_model, un = TRUE)
