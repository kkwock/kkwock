library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

csv_data <- read.csv("/Users/kkwock/OneDrive - Guardant Health/Data/7435_Dilutions/7435-Dilutions_sparkvts.csv", header = FALSE)


names(csv_data) <- csv_data[1,]
csv_data <- csv_data[-1,]

#convert datatype
csv_data <- type.convert(csv_data, as.is=TRUE)

#change levels of groups
csv_data$Group <- factor(csv_data$Group, levels = c("Raw", "DIL1", "DIL2", "DIL3"))


# Scatter plot
plot <- csv_data %>% group_by(Group) %>% ggplot(aes(x=Spark, y=Tapestation, color=Group)) + 
  geom_point() + coord_equal() + theme(axis.ticks = element_blank())  + 
  labs(x = "Spark (nM)", y="Tapestation (nM)", title = "Tapestation Vs. Spark Quantitation") + 
  ylim(0, 300) + geom_hline(yintercept=100, linetype="dashed", color = "red")

plot

plot + geom_text(x = 1050, y =250, label = c(as.expression((R^2~"= 0.4909"))), show.legend = FALSE)



###

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

csv_data <- read.csv("/Users/kkwock/Desktop/LGC_convertedsparkdata.csv", header = FALSE)

names(csv_data) <- csv_data[1,]
csv_data <- csv_data[-1,]

#convert datatype
csv_data <- type.convert(csv_data, as.is=TRUE)

#change levels of groups
csv_data$Group <- factor(csv_data$Group, levels = c("Raw (1.5)", "DIL1 (1.3)", "DIL2 (1.1)", "DIL3 (1.0)", "DIL4 (0.9)", "DIL5 (0.8)"))


# Scatter plot
plot <- csv_data %>% group_by(Group) %>% ggplot(aes(x=Spark, y=Tapestation, color=Group)) + 
  geom_point() + coord_equal() + theme(axis.ticks = element_blank()) + 
  labs(x = "Spark (nM)", y="Tapestation (nM)", title = "Concentration (Spark) vs. PCR Productivity (TS): LGC Spark Actual (nM) x Purity") + 
  geom_hline(yintercept=100, linetype="dashed", color = "red")  + labs(colour="Group (µM)") 

#plot based on type
typeplot <- csv_data %>% group_by(Group) %>% ggplot(aes(x=Spark, y=Tapestation, color=Type)) + 
  geom_point(aes(shape=Group)) + coord_equal() + theme(axis.ticks = element_blank(), legend.position = "bottom") +
  labs(x = "Spark (nM)", y="Tapestation (nM)", title = "Concentration (Spark) vs. PCR Productivity (TS)", subtitle= as.expression(R^2~"=0.6349")) + geom_hline(yintercept=100, linetype="dashed", color = "red")  + labs(colour="Type") 



plot+ stat_smooth(geom="smooth", se=True, method="lm", function = y~m)

plot

histogram <- plot <- csv_data %>% group_by(Group) %>% ggplot(aes(x=Spark, y=Tapestation, color=Group)) + 
  geom_point(aes(shape=Type)) + coord_equal() + theme(axis.ticks = element_blank(), legend.position = "bottom") +
  labs(x = "Spark (nM)", y="Tapestation (nM)", title = "Concentration (Spark) vs. PCR Productivity (TS)", subtitle= "IDT (0018) vs. LGC") + geom_hline(yintercept=100, linetype="dashed", color = "red")  + labs(colour="Group (µM)") 
