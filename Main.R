# Loading of input data

b <- read.table(
  "~/PCAwFCM/Pollen2014.txt",
  sep = ',',
  header = T,
  row.names = 1
)
lb <- read.table("~/PCAwFCM/SupplementaryLabels.txt", sep = ',', header = T)
D <- log2(as.matrix(b) + 1) # log transformation of count data
Input <- t(D) # data matrix, cells in rows, genes in columns
true_tissue_cls <- lb[, 6] # true data partition K=4
true_cell_cls <- lb[, 4]   # true data partition K=11

# Required library

library(ggplot2)
library(ppclust)
library(factoextra)
library(cluster)
library(fclust)
library(mnormt)
require(mgcv)
library("pcaMethods")
library("pcaReduce")
source("~/PCAwFCM/SourceFun.R", echo=TRUE)
# Final output with all cluster labels with different number of sampling run
Output_SF <- Fpca_reduce(Input,
                          nbt = 1,
                          q = 30)


# Evaluation of Rand Index of unequal dimension of Output_SF matrix

N <- length(Output_SF)
M <- dim(Output_SF[[1]])[2]

K11 <- c() # number of cells
K4 <- c()  # number of tissue

for (n in 1:N){
  cls_cell <- c()
  cls_tissue <- c()
  labels <- c()

  for (m in 1:M){
    cls_cell <- c(cls_cell, adjustedRandIndex(Output_SF[[n]][,m], true_cell_cls))
    cls_tissue <- c(cls_tissue, adjustedRandIndex(Output_SF[[n]][,m], true_tissue_cls))
    labels <- c(labels, length(unique(Output_SF[[n]][,m])))
  }

  K11 <- cbind(K11, cls_cell)
  K4 <- cbind(K4, cls_tissue)
}

# Cluster plot on the basis of Rand Index

plot(
  labels,
  K11[, 1],
  col = "cornflowerblue",
  type = "l",
  lty = 3,
  lwd = 0.5,
  main = "",
  xlab = "Number of clusters",
  ylab = "ARANDI score",
  ylim = c(0, 1),
  cex.lab = 1.7,
  cex.axis = 1.5,
  font.lab = 2,
  bty = 'n'
)
for (i in 1:N) {
  lines(
    labels,
    K11[, i],
    col = "cornflowerblue",
    type = "l",
    lty = 3,
    lwd = 1.3
  )
}
for (i in 1:N) {
  lines(
    labels,
    K4[, i],
    col = "yellowgreen",
    type = "l",
    lty = 3,
    lwd = 1.3
  )
}
for (i in 1:N) {
  points(
    labels[which(labels == 11)],
    K11[which(labels == 11), i],
    cex = 1.4,
    pch = 21,
    bg = "cornflowerblue",
    col = "royalblue3"
  )
  points(
    labels[which(labels == 4)],
    K4[which(labels == 4), i],
    cex = 1.4,
    pch = 21,
    bg = "yellowgreen",
    col = "olivedrab4"
  )
}

temp <-
  legend(
    "topright",
    border = NULL,
    legend = c("K = 4", "K = 11"),
    col = c("olivedrab4", "royalblue3"),
    text.width = strwidth("00,000"),
    lty = c(0, 0),
    pch = c(21, 21),
    pt.bg = c("yellowgreen", "cornflowerblue"),
    xjust = 0,
    yjust = 0,
    bty = "n"
  )
dev.off()


