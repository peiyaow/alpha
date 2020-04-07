library(huge)
library(igraph)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(clime)
library(readr)
library(xtable)

setwd("~/Documents/GitHub/alpha/")
load("./data/ADNI2_clean3.RData")
source("./function/main_function.R")
ROI_names <- read_csv("./data/ROI_names.csv", col_names = FALSE)$X1

colnames(X) = ROI_names

n = dim(X)[1]
p = dim(X)[2]

label.level = levels(label)
n_label = length(label.level)

X.list = lapply(label.level, function(l) X[label == l,])
Y.list = lapply(label.level, function(l) Y[label == l])
label.list = lapply(label.level, function(l) label[label == l])
n.vec = as.vector(table(label))

X.mean.list = lapply(X.list, colMeans)
X.list = lapply(1:n_label, function(ix) sweep(X.list[[ix]], 2, X.mean.list[[ix]]))

X2U.list = lapply(1:n_label, function(ix) X2U1(X.list[[ix]], plot = F))

H.list = lapply(X2U.list, function(list) list$H)
P.list = lapply(X2U.list, function(list) list$P)
K.list = lapply(X2U.list, function(list) list$K)

F.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
U.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$U)
L.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$L)
FL.list = lapply(1:n_label, function(ix) F.list[[ix]]%*%L.list[[ix]])

PY.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.list[[ix]])
HY.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.list[[ix]])

# # WLS
# OLS.WLS.U = lm(Y~., data = data.U, weight = w)
# ridge.WLS.U = cv.glmnet(x = U, y = HY, weights = w, alpha = 0)
# beta_u = OLS.WLS.U$coefficients[-1]
# beta_u = coef(ridge.WLS.U, s=ridge.WLS.U$lambda.min)[-1]
# 
# lambda = 0.1
# beta.NC = solve(t(L.list[[1]])%*%L.list[[1]] + lambda*diag(p))%*%t(L.list[[1]])*ml.lm.F.list[[1]]$coefficients[-1]
# 
# beta.AD = solve(t(L.list[[5]])%*%L.list[[5]] + lambda*diag(p))%*%t(L.list[[5]])*ml.lm.F.list[[5]]$coefficients[-1]

# beta.data = data.frame(beta.NC, beta.AD, beta_u)
# pairs(beta.data)
# plot(beta.NC)
# plot(beta.AD)


# plot(beta_u, beta.NC)
# plot(beta_u, beta.AD)
# plot(beta.NC, beta.AD)
# 
# par(mfrow = c(3,1))
# hist(beta.NC)
# hist(beta.AD)
# hist(beta_u)


# ------------------------------- plot PY and HY ----------------------------------
dat2 = data.frame(Y = do.call(c, Y.list), PY = do.call(c, PY.list), HY = do.call(c, HY.list), group = c(rep("NC", n.vec[1]), rep("SMC", n.vec[2]), rep("eMCI", n.vec[3]), rep("lMCI", n.vec[4]), rep("AD", n.vec[5])))
index = 1:n

p1 = ggplot(dat2, aes(x= index, y=Y, color = group)) + geom_point()+ theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p2 = ggplot(dat2, aes(x= index, y=PY, color = group)) + geom_point()+ theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p3 = ggplot(dat2, aes(x= index, y=HY, color = group)) + geom_point()+ theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
grid.arrange(p1, p2, p3, nrow = 3)

# go to tomat.R to plot X

# ----------------------------------- graph ---------------------------------------
plot.adjacency = function(graph.FL.NC, ix, method = "clime"){
  if (method == "mb"){
    igraph.FL.NC = graph_from_adjacency_matrix(as.matrix(graph.FL.NC$path[[ix]]), mode = "undirected")
    node_list <- get.data.frame(igraph.FL.NC, what = "vertices")
    edge_list <- get.data.frame(igraph.FL.NC, what = "edges")
    edge_list1 = edge_list[,2:1]
    colnames(edge_list1) = c("from", "to")
    edge_list = rbind(edge_list, edge_list1)
  }else if (method == "clime"){
    igraph.FL.NC = graph_from_adjacency_matrix(as.matrix(graph.FL.NC$Omegalist[[ix]]!=0), mode = "undirected")
    node_list <- get.data.frame(igraph.FL.NC, what = "vertices")
    edge_list <- get.data.frame(igraph.FL.NC, what = "edges")
    edge_list1 = edge_list[,2:1]
    colnames(edge_list1) = c("from", "to")
    edge_list = rbind(edge_list, edge_list1)
  }
  
  
  g = ggplot(edge_list, aes(x = from, y = to)) +
    geom_raster() +
    theme_bw() +
    # Because we need the x and y axis to display every node,
    # not just the nodes that have connections to each other,
    # make sure that ggplot does not drop unused factor levels
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    theme(
      # Rotate the x-axis lables so they are legible
      axis.text.x = element_text(angle = 270, hjust = 0),
      # Force the plot into a square aspect ratio
      aspect.ratio = 1,
      # Hide the legend (optional)
      legend.position = "none", 
      axis.title.x=element_blank(),
      axis.title.y=element_blank())
  return(g)
}

cov.FL.NC = t(FL.list[[1]])%*%FL.list[[1]]/n.vec[1]
cov.X.NC = t(X.list[[1]])%*%X.list[[1]]/n.vec[1]
cov.U.NC = t(U.list[[1]])%*%U.list[[1]]/n.vec[1]

cov.FL.SMC = t(FL.list[[2]])%*%FL.list[[2]]/n.vec[2]
cov.X.SMC = t(X.list[[2]])%*%X.list[[2]]/n.vec[2]
cov.U.SMC = t(U.list[[2]])%*%U.list[[2]]/n.vec[2]

cov.FL.eMCI = t(FL.list[[3]])%*%FL.list[[3]]/n.vec[3]
cov.X.eMCI = t(X.list[[3]])%*%X.list[[3]]/n.vec[3]
cov.U.eMCI = t(U.list[[3]])%*%U.list[[3]]/n.vec[3]

cov.FL.lMCI = t(FL.list[[4]])%*%FL.list[[4]]/n.vec[4]
cov.X.lMCI = t(X.list[[4]])%*%X.list[[4]]/n.vec[4]
cov.U.lMCI = t(U.list[[4]])%*%U.list[[4]]/n.vec[4]

cov.FL.AD = t(FL.list[[5]])%*%FL.list[[5]]/n.vec[5]
cov.X.AD = t(X.list[[5]])%*%X.list[[5]]/n.vec[5]
cov.U.AD = t(U.list[[5]])%*%U.list[[5]]/n.vec[5]

# cov.U = t(U)%*%U

# # huge
# graph.FL.NC = huge(cov.FL.NC, sym = "and", method = "mb")
# graph.X.NC = huge(cov.X.NC, sym = "and", method = "mb")
# graph.U.NC = huge(cov.U.NC, sym = "and", method = "mb")
# 
# graph.X.SMC = huge(cov.X.SMC, sym = "and", method = "mb")
# graph.U.SMC = huge(cov.U.SMC, sym = "and", method = "mb")
# 
# graph.X.eMCI = huge(cov.X.eMCI, sym = "and", method = "mb")
# graph.U.eMCI = huge(cov.U.eMCI, sym = "and", method = "mb")
# 
# graph.X.lMCI = huge(cov.X.lMCI, sym = "and", method = "mb")
# graph.U.lMCI = huge(cov.U.lMCI, sym = "and", method = "mb")
# 
# graph.FL.AD = huge(cov.FL.AD, sym = "and", method = "mb")
# graph.X.AD = huge(cov.X.AD, sym = "and", method = "mb")
# graph.U.AD = huge(cov.U.AD, sym = "and", method = "mb")

# clime
standardize = F
linsolver = "simplex"
graph.FL.NC = clime(x = FL.list[[1]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = .1)
plot.adjacency(graph.FL.NC, ix = 1)



graph.X.NC = clime(x = X.list[[1]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)
graph.U.NC = clime(x = U.list[[1]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)
plot.adjacency(graph.FL.NC, ix = 3)
graph.FL.SMC = clime(x = FL.list[[2]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)
plot.adjacency(graph.FL.SMC, ix = 6)

graph.X.SMC = clime(x = X.list[[2]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)
graph.U.SMC = clime(x = U.list[[2]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)

graph.FL.eMCI = clime(x = FL.list[[3]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)
plot.adjacency(graph.FL.eMCI, ix = 6)

graph.X.eMCI = clime(x = X.list[[3]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)
graph.U.eMCI = clime(x = U.list[[3]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)

graph.FL.lMCI = clime(x = FL.list[[4]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)
graph.X.lMCI = clime(x = X.list[[4]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)
graph.U.lMCI = clime(x = U.list[[4]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)

graph.FL.AD = clime(x = FL.list[[5]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)
graph.X.AD = clime(x = X.list[[5]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)
graph.U.AD = clime(x = U.list[[5]], nlambda = 10, standardize = standardize, linsolver = linsolver, lambda.min = 0.1)

plot.adjacency(graph.FL.NC, ix = 5)
plot.adjacency(graph.FL.AD, ix = 5)
4

lambda_ix = 6
grid.arrange(plot.adjacency(graph.X.NC, ix = lambda_ix),
             plot.adjacency(graph.X.SMC, ix = lambda_ix),
             plot.adjacency(graph.X.eMCI, ix = lambda_ix),
             plot.adjacency(graph.X.lMCI, ix = lambda_ix),
             plot.adjacency(graph.X.AD, ix = lambda_ix),
             plot.adjacency(graph.FL.NC, ix = lambda_ix),
             plot.adjacency(graph.FL.SMC, ix = lambda_ix),
             plot.adjacency(graph.FL.eMCI, ix = lambda_ix),
             plot.adjacency(graph.FL.lMCI, ix = lambda_ix),
             plot.adjacency(graph.FL.AD, ix = lambda_ix),
             plot.adjacency(graph.U.NC, ix = lambda_ix),
             plot.adjacency(graph.U.SMC, ix = lambda_ix),
             plot.adjacency(graph.U.eMCI, ix = lambda_ix),
             plot.adjacency(graph.U.lMCI, ix = lambda_ix),
             plot.adjacency(graph.U.AD, ix = lambda_ix),
             nrow = 3)
              
grid.arrange(plot.adjacency(graph.X.NC, ix = 5),
             plot.adjacency(graph.X.SMC, ix = 5),
             plot.adjacency(graph.X.eMCI, ix = 5),
             plot.adjacency(graph.X.lMCI, ix = 6),
             plot.adjacency(graph.X.AD, ix = 6),
             plot.adjacency(graph.U.NC, ix = 6),
             plot.adjacency(graph.U.SMC, ix = 6),
             plot.adjacency(graph.U.eMCI, ix = 6),
             plot.adjacency(graph.U.lMCI, ix = 6),
             plot.adjacency(graph.U.AD, ix = 6),
             plot.adjacency(graph.FL.NC, ix = 6),
             plot.adjacency(graph.FL.SMC, ix = 6),
             plot.adjacency(graph.FL.eMCI, ix = 6),
             plot.adjacency(graph.FL.lMCI, ix = 6),
             plot.adjacency(graph.FL.AD, ix = 6),
             nrow = 3)

lambda_ix = 6
# upper.tri(matrix(1, nrow = p, ncol = p))
ROI_matrix = cbind(ROI_names[which((graph.FL.AD$Omegalist[[lambda_ix]]==0)*(graph.FL.NC$Omegalist[[lambda_ix]]!=0)*upper.tri(matrix(1, nrow = p, ncol = p))!=0, arr.ind = T)[,1]], 
                   ROI_names[which((graph.FL.AD$Omegalist[[lambda_ix]]==0)*(graph.FL.NC$Omegalist[[lambda_ix]]!=0)*upper.tri(matrix(1, nrow = p, ncol = p))!=0, arr.ind = T)[,2]])

xtable(ROI_matrix)

lambda_ix = 6
grid.arrange(plot.adjacency(graph.X.NC, ix = lambda_ix),
             plot.adjacency(graph.X.SMC, ix = lambda_ix), 
             plot.adjacency(graph.X.eMCI, ix = lambda_ix),
             plot.adjacency(graph.X.lMCI, ix = lambda_ix),
             plot.adjacency(graph.X.AD, ix = lambda_ix),
             plot.adjacency(graph.U.NC, ix = lambda_ix),
             plot.adjacency(graph.U.SMC, ix = lambda_ix), 
             plot.adjacency(graph.U.eMCI, ix = lambda_ix),
             plot.adjacency(graph.U.lMCI, ix = lambda_ix),
             plot.adjacency(graph.U.AD, ix = lambda_ix),
             #plot.adjacency(graph.U.NC, ix = lambda_ix),
             #             plot.adjacency(graph.U.SMC, ix = lambda_ix), 
             #             plot.adjacency(graph.U.eMCI, ix = lambda_ix),
             #             plot.adjacency(graph.U.lMCI, ix = lambda_ix),
             #plot.adjacency(graph.U.AD, ix = lambda_ix),
             #plot.adjacency(graph.FL.NC, ix = lambda_ix),
             #             plot.adjacency(graph.U.SMC, ix = lambda_ix), 
             #             plot.adjacency(graph.U.eMCI, ix = lambda_ix),
             #             plot.adjacency(graph.U.lMCI, ix = lambda_ix),
             #plot.adjacency(graph.FL.AD, ix = lambda_ix),
             nrow = 2)



# g.k = graph.adjacency(as.matrix(graph.FL.NC$path[[2]]), mode = "undirected", 
#                       diag = FALSE)
# plot.igraph(g.k, edge.arrow.size=0.05, vertex.size = 5, vertex.label.cex = 0.5)
# 
# g.k = graph.adjacency(as.matrix(graph.FL.AD$path[[2]]), mode = "undirected", 
#                       diag = FALSE)
# plot.igraph(g.k, edge.arrow.size=0.05, vertex.size = 5, vertex.label.cex = 0.5)
# 
# g.k = graph.adjacency(as.matrix(graph.U$path[[2]]), mode = "undirected", 
#                       diag = FALSE)
# plot.igraph(g.k, edge.arrow.size=0.05, vertex.size = 5, vertex.label.cex = 0.5)



