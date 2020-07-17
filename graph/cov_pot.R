library(gridExtra)
library(ggplot2)
library(readr)
library(reshape2)
library(scico)

setwd("~/Documents/GitHub/alpha/")
load("./data/ADNI2_clean3.RData")
source("./function/main_function.R")
source("./function/sim_function.R")

ROI_names <- read_csv("./data/ROI_names.csv", col_names = FALSE)$X1
#colnames(X) = ROI_names
colnames(X) = NULL

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



cov.FL.NC = t(FL.list[[1]])%*%FL.list[[1]]/n.vec[1]
cov.U.NC = POET(t(U.list[[1]]), K = 0, C = 0.5, thres = "soft", matrix = "vad")$SigmaY
cov.X.NC = cov.U.NC + cov.FL.NC


cov.FL.SMC = t(FL.list[[2]])%*%FL.list[[2]]/n.vec[2]
cov.U.SMC = POET(t(U.list[[2]]), K = 0, C = 0.5, thres = "soft", matrix = "vad")$SigmaY
cov.X.SMC = cov.U.SMC + cov.FL.SMC


cov.FL.eMCI = t(FL.list[[3]])%*%FL.list[[3]]/n.vec[3]
cov.U.eMCI = POET(t(U.list[[3]]), K = 0, C = 0.5, thres = "soft", matrix = "vad")$SigmaY
cov.X.eMCI = cov.U.eMCI + cov.FL.eMCI


cov.FL.lMCI = t(FL.list[[4]])%*%FL.list[[4]]/n.vec[4]
cov.U.lMCI = POET(t(U.list[[4]]), K = 0, C = 0.5, thres = "soft", matrix = "vad")$SigmaY
cov.X.lMCI = cov.U.lMCI + cov.FL.lMCI


cov.FL.AD = t(FL.list[[5]])%*%FL.list[[5]]/n.vec[5]
cov.U.AD = POET(t(U.list[[5]]), K = 0, C = 0.5, thres = "soft", matrix = "vad")$SigmaY
cov.X.AD = cov.U.AD + cov.FL.AD

index = which(as.logical(upper.tri(cov.FL.NC) + diag(p)))

summary(cov.U.NC[index])
summary(cov.U.SMC[index])
summary(cov.U.eMCI[index])
summary(cov.U.lMCI[index])
summary(cov.U.AD[index])
# (-0.5 , 2.5)

hist(cov.X.NC[index])
hist(cov.X.SMC[index])
hist(cov.X.eMCI[index])
hist(cov.X.lMCI[index])
hist(cov.X.AD[index])

summary(cov.X.NC[index])
summary(cov.X.SMC[index])
summary(cov.X.eMCI[index])
summary(cov.X.lMCI[index])
summary(cov.X.AD[index])
# (-0.5, 2.5)

summary(cov.FL.NC[index])
summary(cov.FL.SMC[index])
summary(cov.FL.eMCI[index])
summary(cov.FL.lMCI[index])
summary(cov.FL.AD[index])
# (-0.5, 0.55)


# U
# changing scale
cov.U.NC[cov.U.NC > 0] = (cov.U.NC[cov.U.NC > 0])/2.5*0.5
cov.U.SMC[cov.U.SMC > 0] = (cov.U.SMC[cov.U.SMC > 0])/2.5*0.5
cov.U.eMCI[cov.U.eMCI > 0] = (cov.U.eMCI[cov.U.eMCI > 0])/2.5*0.5
cov.U.lMCI[cov.U.lMCI > 0] = (cov.U.lMCI[cov.U.lMCI > 0])/2.5*0.5
cov.U.AD[cov.U.AD > 0] = (cov.U.AD[cov.U.AD > 0])/2.5*0.5

melted_cov.U.NC <- melt(cov.U.NC)
melted_cov.U.SMC <- melt(cov.U.SMC)
melted_cov.U.eMCI <- melt(cov.U.eMCI)
melted_cov.U.lMCI <- melt(cov.U.lMCI)
melted_cov.U.AD <- melt(cov.U.AD)

# U
# range for cov.U across groups (-0.5 , 2.5)
pu.NC = ggplot(data = melted_cov.U.NC, aes(x=Var1, y = Var2, fill=value)) + geom_tile(show.legend = F) + scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", lim = c(-0.5, 0.5)) + xlab("ROI Index") + ylab("ROI Index") 

pu.SMC = ggplot(data = melted_cov.U.SMC, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = F)+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", limits = c(-0.5, 0.5))

pu.eMCI = ggplot(data = melted_cov.U.eMCI, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = F)+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", limits = c(-0.5, 0.5))+ xlab("ROI Index") + ylab("ROI Index") 

pu.lMCI = ggplot(data = melted_cov.U.lMCI, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = F)+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", limits = c(-0.5, 0.5))

pu.AD = ggplot(data = melted_cov.U.AD, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = F)+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", limits = c(-0.5, 0.5))+ xlab("ROI Index") + ylab("ROI Index") 


# changing scale
# range for cov.FL across groups (-0.5 , 0.55)
cov.FL.NC[cov.FL.NC > 0] = (cov.FL.NC[cov.FL.NC > 0])/.55*0.5
cov.FL.SMC[cov.FL.SMC > 0] = (cov.FL.SMC[cov.FL.SMC > 0])/.55*0.5
cov.FL.eMCI[cov.FL.eMCI > 0] = (cov.FL.eMCI[cov.FL.eMCI > 0])/.55*0.5
cov.FL.lMCI[cov.FL.lMCI > 0] = (cov.FL.lMCI[cov.FL.lMCI > 0])/.55*0.5
cov.FL.AD[cov.FL.AD > 0] = (cov.FL.AD[cov.FL.AD > 0])/.55*0.5

melted_cov.FL.NC <- melt(cov.FL.NC)
melted_cov.FL.SMC <- melt(cov.FL.SMC)
melted_cov.FL.eMCI <- melt(cov.FL.eMCI)
melted_cov.FL.lMCI <- melt(cov.FL.lMCI)
melted_cov.FL.AD <- melt(cov.FL.AD)

ggplot(data = melted_cov.FL.NC, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 0.55), palette = "roma", lim = c(-0.5, 0.5)) 

ggplot(data = melted_cov.FL.SMC, aes(x=Var1, y=Var2, fill=value)) + geom_tile()+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 0.55), palette = "roma", limits = c(-0.5, 0.5))

ggplot(data = melted_cov.FL.eMCI, aes(x=Var1, y=Var2, fill=value)) + geom_tile()+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 0.55), palette = "roma", limits = c(-0.5, 0.5))

ggplot(data = melted_cov.FL.lMCI, aes(x=Var1, y=Var2, fill=value)) + geom_tile()+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 0.55), palette = "roma", limits = c(-0.5, 0.5))

ggplot(data = melted_cov.FL.AD, aes(x=Var1, y=Var2, fill=value)) + geom_tile()+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 0.55), palette = "roma", limits = c(-0.5, 0.5))

# X
# changing scale
cov.X.NC[cov.X.NC > 0] = (cov.X.NC[cov.X.NC > 0])/2.5*0.5
cov.X.SMC[cov.X.SMC > 0] = (cov.X.SMC[cov.X.SMC > 0])/2.5*0.5
cov.X.eMCI[cov.X.eMCI > 0] = (cov.X.eMCI[cov.X.eMCI > 0])/2.5*0.5
cov.X.lMCI[cov.X.lMCI > 0] = (cov.X.lMCI[cov.X.lMCI > 0])/2.5*0.5
cov.X.AD[cov.X.AD > 0] = (cov.X.AD[cov.X.AD > 0])/2.5*0.5


melted_cov.X.NC <- melt(cov.X.NC)
melted_cov.X.SMC <- melt(cov.X.SMC)
melted_cov.X.eMCI <- melt(cov.X.eMCI)
melted_cov.X.lMCI <- melt(cov.X.lMCI)
melted_cov.X.AD <- melt(cov.X.AD)

# X
# range for cov.X across groups (-0.5 , 2.5)
px.NC = ggplot(data = melted_cov.X.NC, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = F) + scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", lim = c(-0.5, 0.5)) + xlab("ROI Index") + ylab("ROI Index") 

px.SMC = ggplot(data = melted_cov.X.SMC, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = F)+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", limits = c(-0.5, 0.5))

px.eMCI = ggplot(data = melted_cov.X.eMCI, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = F)+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", limits = c(-0.5, 0.5))+ xlab("ROI Index") + ylab("ROI Index") 

px.lMCI = ggplot(data = melted_cov.X.lMCI, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = F)+ scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", limits = c(-0.5, 0.5))

px.AD = ggplot(data = melted_cov.X.AD, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = T) + 
  scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", limits = c(-0.5, 0.5),
                   guide = guide_colorbar(title = NULL,
                                          label = T,
                                          draw.ulim = TRUE, 
                                          draw.llim = TRUE,
                                          # here comes the code change:
                                          frame.colour = NULL,
                                          ticks = TRUE, 
                                          nbin = 10,
                                          label.position = "bottom",
                                          barwidth =20,
                                          barheight = 1, 
                                          direction = 'horizontal'))+ xlab("ROI Index") + ylab("ROI Index") +
  theme(legend.text = element_text(size=15))


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(px.AD)

px.AD = ggplot(data = melted_cov.X.AD, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = F) + 
  scale_fill_scico(breaks = c(-0.5, 0, 0.5), labels=c(-0.5, 0, 2.5), palette = "roma", limits = c(-0.5, 0.5),
                   guide = guide_colorbar(title = NULL,
                                          label = T,
                                          draw.ulim = TRUE, 
                                          draw.llim = TRUE,
                                          # here comes the code change:
                                          frame.colour = NULL,
                                          ticks = TRUE, 
                                          nbin = 10,
                                          label.position = "right",
                                          barwidth = 1.3,
                                          barheight = 13, 
                                          direction = 'vertical'))+ xlab("ROI Index") + ylab("ROI Index") 


# cor.FL.NC = cov2cor(cov.FL.NC)
# cor.X.NC = cov2cor(cov.X.NC)
# cor.U.NC = cov2cor(cov.U.NC)
# 
# cor.FL.SMC = cov2cor(cov.FL.SMC)
# cor.X.SMC = cov2cor(cov.X.SMC)
# cor.U.SMC = cov2cor(cov.U.SMC)
# 
# cor.FL.eMCI = cov2cor(cov.FL.eMCI)
# cor.X.eMCI = cov2cor(cov.X.eMCI)
# cor.U.eMCI = cov2cor(cov.U.eMCI)
# 
# cor.FL.lMCI = cov2cor(cov.FL.lMCI)
# cor.X.lMCI = cov2cor(cov.X.lMCI)
# cor.U.lMCI = cov2cor(cov.U.lMCI)
# 
# cor.FL.AD = cov2cor(cov.FL.AD)
# cor.X.AD = cov2cor(cov.X.AD)
# cor.U.AD = cov2cor(cov.U.AD)



# grid.arrange(px.NC, px.SMC, px.eMCI, px.lMCI, px.AD, 
#              pu.NC, pu.SMC, pu.eMCI, pu.lMCI, pu.AD, 
#              nrow = 2)

grid.arrange(px.NC, px.eMCI, px.AD, 
             pu.NC, pu.eMCI, pu.AD, 
             mylegend,
             nrow = 3, heights = c(10,10,1), layout_matrix = rbind(c(1,2,3), c(4,5,6), c(7,7,7)))






