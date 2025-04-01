args <- commandArgs(trailingOnly = TRUE)
cat("Received arguments:\n")
print(args)

if (length(args) < 3) stop("Expected 3 arguments: expr_file, label_file, output_plot")

expr_file <- args[1]
label_file <- args[2]
output_plot <- args[3]
output_file <- args[4]  # <- RDS output path
cat("Checking file existence:\n")
print(file.exists(expr_file))
print(file.exists(label_file))

cat("Loading required packages...\n")
# ... your package install_if_missing() logic here ...

# Loading of input data

b <- read.table(expr_file, sep = ',', header = TRUE, row.names = 1)
lb <- read.table(label_file, sep = ',', header = TRUE)

# b <- read.table(expr_file, sep = '\t', header = TRUE, row.names = 1)
# lb <- read.table(label_file, sep = '\t', header = TRUE)


D <- log2(as.matrix(b) + 1) # log transformation of count data
Input <- t(D) # data matrix, cells in rows, genes in columns
true_tissue_cls <- lb[, 6] # true data partition K=4
true_cell_cls <- lb[, 4]   # true data partition K=11

required_packages <- c(
  "ggplot2", "ppclust", "factoextra", "cluster", "fclust", 
  "mnormt", "mgcv", "pcaMethods", "pcaReduce", "mclust"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))

    if (pkg == "pcaReduce") {
      # Install from GitHub if pcaReduce
      if (!requireNamespace("devtools", quietly = TRUE))
        install.packages("devtools", repos = "https://cloud.r-project.org")
      devtools::install_github("JustinaZ/pcaReduce")
    } else {
      tryCatch({
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }, error = function(e) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        BiocManager::install(pkg, ask = FALSE)
      })
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


# Install and load all
invisible(lapply(required_packages, install_if_missing))

# Load external functions
source("~/PCAwFCM/scripts/SourceFun.R", echo=TRUE) # adjust path if needed

# Run clustering to get Final output with all cluster labels with various sampling run (nbt)
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
# Plot results
png(output_plot, width=1000, height=800)

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

saveRDS(Output_SF, file = output_file)
