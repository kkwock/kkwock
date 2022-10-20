library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

csv_data <- read.csv("/Users/kkwock/OneDrive - Guardant Health/Data/2021Nov_lgc_tsvspark.csv", header = FALSE)

names(csv_data) <- csv_data[1,]
csv_data <- csv_data[-1,]

#convert datatype
csv_data <- type.convert(csv_data, as.is=TRUE)

#change levels of groups
csv_data$Group <- factor(csv_data$Group, levels = c("Raw", "DIL1", "DIL2", "DIL3", "DIL4", "DIL5"))


# Scatter plot
plot <- csv_data %>% group_by(Group) %>% ggplot(aes(x=Spark, y=Tapestation, color=Group)) + 
  geom_point() + coord_equal() + theme(axis.ticks = element_blank()) + 
  stat_smooth(method = "lm", se=FALSE, color="orange", formula = y ~ x) + 
  labs(x = "Spark (nM)", y="Tapestation (nM)", title = "Tapestation Vs. Spark Quantitation") + 
  ylim(0, 300) + geom_hline(yintercept=100, linetype="dashed", color = "red")
  
  
plot + geom_text(x = 1050, y =250, label = c(as.expression((R^2~"= 0.4909"))), show.legend = FALSE)
