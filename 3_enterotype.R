#
library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(readxl)
library(Cairo)
library(tidyverse)

##### k=5  MDD
data <-  read.delim ("DMM.txt")
data <- as.data.frame(data)
rownames(data) = data[,1]
data <- data[,-1]
data1=data[c(grep("D",colnames(data)))]

count1 <- as.matrix(t(data1))

fit1 <- lapply(1:5, dmn, count = count1, verbose=TRUE)

lplc <- sapply(fit1, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit1, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit1, BIC) # AIC / BIC / Laplace
best <- fit1[[which.min(unlist(lplc))]]
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")

mixturewt(best)

ass <- apply(mixture(best), 1, which.max)

d1 <- melt(fitted(best))
colnames(d1) <- c("OTU", "cluster", "value")
d1 <- subset(d1, cluster == 1) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

d2 <- melt(fitted(best))
colnames(d2) <- c("OTU", "cluster", "value")
d2 <- subset(d2, cluster == 2) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

d3 <- melt(fitted(best))
colnames(d3) <- c("OTU", "cluster", "value")
d3 <- subset(d3, cluster == 3) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

d4 <- melt(fitted(best))
colnames(d4) <- c("OTU", "cluster", "value")
d4 <- subset(d4, cluster == 4) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

d5 <- melt(fitted(best))
colnames(d5) <- c("OTU", "cluster", "value")
d5 <- subset(d5, cluster == 5) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

write_csv(d1,"D1.csv")
write_csv(d2,"D2.csv")
write_csv(d3,"D3.csv")
write_csv(d4,"D4.csv")
write_csv(d5,"D5.csv")

gp <- best@group %>% as.data.frame()
write.csv(gp,"D_best_group.csv",row.names=TRUE)

#HC+MDD
##### k=5 C D
data <-  read.delim("DMM.txt")
data <- as.data.frame(data)
rownames(data) = data[,1]
data <- data[,-1]

count <- as.matrix(t(data))

fit <- lapply(1:5, dmn, count = count, verbose=TRUE)
# AIC / BIC / Laplace
lplc <- sapply(fit, laplace)
aic  <- sapply(fit, AIC) 
bic  <- sapply(fit, BIC) 
best <- fit[[which.min(unlist(lplc))]]
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")

mixturewt(best)

ass <- apply(mixture(best), 1, which.max)

d1 <- melt(fitted(best))
colnames(d1) <- c("OTU", "cluster", "value")
d1 <- subset(d1, cluster == 1) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

d2 <- melt(fitted(best))
colnames(d2) <- c("OTU", "cluster", "value")
d2 <- subset(d2, cluster == 2) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

d3 <- melt(fitted(best))
colnames(d3) <- c("OTU", "cluster", "value")
d3 <- subset(d3, cluster == 3) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

d4 <- melt(fitted(best))
colnames(d4) <- c("OTU", "cluster", "value")
d4 <- subset(d4, cluster == 4) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

d5 <- melt(fitted(best))
colnames(d5) <- c("OTU", "cluster", "value")
d5 <- subset(d5, cluster == 5) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.8)) 

write_csv(d1,"CD1.csv")
write_csv(d2,"CD2.csv")
write_csv(d3,"CD3.csv")
write_csv(d4,"CD4.csv")
write_csv(d5,"CD5.csv")

gp2 <- best@group %>% as.data.frame()
write.csv(gp2,"CD_best_group.csv",row.names=TRUE)
