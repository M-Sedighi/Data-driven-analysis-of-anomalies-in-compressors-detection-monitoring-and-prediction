

##############################################
# R script for asset anomaly detection method
##############################################

# defining the input training dataframe, input timeframe, and known failures dataframe

X <- data_B
gtime <- time_B
fX <- subset(data_B,class_B==1)

# preprocessing

Xs <- scale(X)
d <- ncol(Xs)
names(X)[d]<-paste("12_MDI")
f <- nrow (fX)
Xm <- matrix(ncol=d, nrow=1)
Xdev <- matrix(ncol=d, nrow=1)
for (i in 1:d){
  Xm[,i] <- mean(X[,i])
  Xdev[,i] <- sd(X[,i])
}

# process variables definition 

a <- 0.95   #control charts threshold
SF <- 0.91   #GSOM spread factor
AT <- 3   #BKDE anomaly score threshold 


# Color palette definition

pretty_palette <- c('lightskyblue', 'red3', 'slateblue3', 'darkorange1',
                    'mediumseagreen', 'gold1', 'chocolate4','salmon','navy')

heatc <- c('green','red')

#################################################
# Bivariate Kernel Density Estimation
#################################################

library(KernSmooth)
library(lattice)
library(fields)
library(rgl)
library(scatterD3)
library(dygraphs)
library(R.basic)

# BKDE distributions for every possible combination of two variables in dataframe
Y <- X
n <- nrow(Y)
LI=list()
Lii <- data.frame(NULL)
k=1
for (i in 2:d){
  for (j in 1:(i-1)){
    h1=sd(Y[,j])/(n^(1/6))
    h2=sd(Y[,i])/(n^(1/6))
    M=cbind(Y[,j],Y[,i])
    L=bkde2D(M,bandwidth=c(h1,h2),gridsize=c(256,256))
    LI[[k]]=L
    t=interp.surface(list(x=L$x1,y=L$x2,z=L$fhat),loc=M)
    tstar=max(L$fhat)
    AST=-log(t/tstar)
    Lii <- rbind(Lii,new=c(i,j,k))
    k=k+1
  } 
}
k=1
Lii <- as.matrix(Lii)
colnames(Lii) <- c('i','j','k')

#################################################
# Growing Self-organizing map
#################################################

library(GrowingSOM)
library(KernSmooth)
library(lattice)
library(fields)
library(rgl)
library(scatterD3)
library(dygraphs)

# training the GSOM and plotting the related figures

Gtrain_map <- train.gsom(Xs, spreadFactor=SF)

jpeg(filename = paste0("GSOM_2grid.jpg"),
     width = 650, height = 650, units = "px", pointsize = 18,
     quality = 100)
plot(Gtrain_map, type="count")
dev.off()

# Some Plots
jpeg(filename = paste0("GSOM_training2.jpg"),
     width = 800, height = 600, units = "px", pointsize = 18,
     quality = 100)
plot(Gtrain_map, type = "training",lwd=2.5)
dev.off()

# Estimate the number of clusters that suits the grid

wss <- (nrow(Gtrain_map$nodes$codes)-1)*sum(apply(Gtrain_map$nodes$codes,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(Gtrain_map$nodes$codes,
                                     centers=i)$withinss)
par(mar=c(5.1,4.1,4.1,2.1))
jpeg(filename = paste0("WCSS.jpg"),
     width = 600, height = 600, units = "px", pointsize = 16,
     quality = 100)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares", main="Within cluster sum of squares (WCSS)")
dev.off()

gcodesmat <- Gtrain_map$nodes$codes
Gtest_map <- map.gsom(Gtrain_map, Xs)
Gtest_node <- Gtest_map$mapped$bmn

fXs <- matrix(nrow=f,ncol=d)
for (j in 1:f){
  for (i in 1:d){
    fXs[j,i] <- (fX[j,i]-Xm[,i])/Xdev[,i]
  }
}

Gfail_map <- map.gsom(Gtrain_map, fXs)
Gfail_node <- Gfail_map$mapped$bmn

# comparison of the data time-series with different clustering variables

CLUSTER_data <- list()
for (i in 1:9){
  gsom_cluster <- cutree(hclust(dist(gcodesmat)), i)
  Gtest_cluster <- gsom_cluster[Gtest_node]
  Gtest_col <- pretty_palette[Gtest_cluster]
  CLUSTER_data[[i]] <- Gtest_cluster
  jpeg(filename = paste0("MDI_",i,"_clusters.jpg"),
       width = 800, height = 600, units = "px", pointsize = 16,
       quality = 100)
  plot(gtime, X[,d], pch=20 ,col=Gtest_col)
  dev.off()
}

NC <- 8

# clustering the grid

gsom_cluster <- cutree(hclust(dist(gcodesmat)), NC)
gsom_col <- pretty_palette[gsom_cluster]
Gtest_cluster <- gsom_cluster[Gtest_node]
Gtest_col <- pretty_palette[Gtest_cluster]
Gfail_cluster <- gsom_cluster[Gfail_node]
gcodesreal <- gcodesmat
for (j in 1:d){
  gcodesreal[,j] <- (gcodesmat[,j]*Xdev[,j])+Xm[,j]
}

Nd <-nrow(Gtrain_map$nodes$position) 


hc <- hclust(dist(gcodesmat))
dend <- as.dendrogram(hc,hang=0)


jpeg(filename = paste0("dendro15.jpg"),
     width = 1100, height = 1250, units = "px", pointsize = 20,
     quality = 100)
plot(as.phylo(hc), type = "fan", tip.color = gsom_col,edge.width=3, label.offset = 0.3, 
     cex = 1.3, col = "red",no.margin=TRUE,font=2)
dev.off()
# Finding the center of each cluster and checking for their anomaly according to BKDE

cl_mean <- matrix(data=0, nrow=NC, ncol=d)
Cl_mat <- as.matrix(gsom_cluster)
for (i in 1:NC){
  cl <- subset(gcodesreal, Cl_mat==i)
  for (j in 1:d){
    cl_mean[i,j] <- mean(cl[,j])
  }
}

Cl_check <- matrix(data=0, nrow=NC,ncol=nrow(Lii))
Cl_val <- matrix(data=0, nrow=NC,ncol=nrow(Lii))

k=1
for (i in 2:d){
  for (j in 1:(i-1)){
    M=cbind(cl_mean[,j],cl_mean[,i])
    L=LI[[k]]
    t=interp.surface(list(x=L$x1,y=L$x2,z=L$fhat),loc=M)
    tstar=max(L$fhat)
    AST=-log(t/tstar)
    Cl_check[,k] <- AST > AT 
    Cl_val[,k] <- AST
    k=k+1
  } 
}
k=1

Cl_check[is.na(Cl_check)] <- 1
Cl_con <- list()
for (c in 1:NC){
  cll_ch <- matrix(data = NA,ncol=d,nrow=d)
  for (i in 2:d){
    for (j in 1:(i-1)){
      cll_ch[j,i] <- Cl_check[c,k]
      k <- k+1
    }
  }
  k=1
  Cl_con[[c]] <- cll_ch
  colnames(Cl_con[[c]]) <- colnames(X)
  rownames(Cl_con[[c]]) <- colnames(X)
}

# plotting the cluster anomaly control boards based on BKDE

for (i in 1:NC){
  
  test <- Cl_con[[i]]
  jpeg(filename = paste0("Cluster2_",i,"_control.jpg"),
       width = 450, height = 450, units = "px", pointsize = 16,
       quality = 100)
  heatmap(t(test),col=heatc,Rowv = NA,Colv = NA,symm = TRUE,
          main=paste0("Cluster ",i," Control"),margins = c(6, 6))
  dev.off()
}

# plotting GSOM nodes on top of variable biplot with color coding for clusters

for (i in 2:d){
  for (j in 1:(i-1)){
    jpeg(filename = paste0("GSOM_variables_",j,"_",i,".jpg"),
         width = 800, height = 600, units = "px", pointsize = 16,
         quality = 100)
    plot(X[,c(j,i)],pch='.',col='gray40')
    points(gcodesreal[,c(j,i)],pch=paste0(gsom_cluster),cex=1.5,col=gsom_col)
    dev.off()
  }
}

####################################
# Principal component analysis
####################################

library(ICSNP)
library(devtools)
library(ggplot2)
library(ggbiplot)

# performing PCA on the training set a whole

m <- nrow(X)
T2 <- matrix(ncol=1, nrow=n)
Q <- matrix(ncol=1, nrow=n)
comp_PCA <- princomp (X, scores=TRUE, cor=TRUE)


jpeg(filename = "PCs.jpg",
     width = 500, height = 500, units = "px", pointsize = 16,
     quality = 100)
plot(comp_PCA,main="PCA components")
dev.off()

k <- 0
c <- comp_PCA$sdev[1]

while (c >= 1){
  k <- k+1
  c <- comp_PCA$sdev[k+1]
}
if (k==1){k <- 2}

jpeg(filename = "PCs_sdesv.jpg",
     width = 600, height = 500, units = "px", pointsize = 15,
     quality = 100)
plot(comp_PCA$sdev, pch=16,cex=1.6,col="blue", ylab="Standard deviation", xlab="PC count", cex.lab=1.4,
     main="Retained Principal components")
abline(h=1, col="red", lwd=2)
text(comp_PCA$sdev[1:k], labels=round(comp_PCA$sdev[1:k], digits = 3), cex= 1, pos=4)
dev.off()

# Calculating the parameters for control charts

V <- (comp_PCA$sdev^2)
P <- comp_PCA$loadings[,1:k]
Pt <- t(P)
L <- diag(V[1:k])
L1 <- solve(L)
TR3 <- diag(d)
F <- qf(a, df1=k, df2=m-k) 
Ta2 <- (k*((m^2)-1)/(m*(m-k)))*F


teta <- matrix(data=0,ncol=1, nrow=3)
for (i in 1:3){
  for (j in (k+1):d){
    teta[i] <- teta[i]+V[j]^i
  }
}

h0 <- 1-((2*(teta[1]*teta[3]))/(3*(teta[2]^2)))
za <- qnorm (a, mean=0, sd=1)
Qa <- teta[1]*(((za*((2*teta[2]*(h0^2))^0.5)/teta[1])+1+((teta[2]*h0*(1-h0))/(teta[2]^2)))^(1/h0))


for (j in 1:n){
  Bt <- Xs[j,]
  B <- t(Bt)
  
  #T2 Hotelling statistic
  Y <- Pt %*% Bt
  Yt <- t(Y)
  TR <- L1 %*% Y
  T2[j] <- Yt %*% TR
  #Q-statistic
  TR2 <- P %*% Pt
  TR4 <- TR3 - TR2
  TR5 <- TR4 %*% Bt
  Q[j,] <- B %*% TR5
}

# plotting the control charts for all of the data

jpeg(filename = "Hotelling.jpg",
     width = 1200, height = 600, units = "px", pointsize = 18,
     quality = 100)
plot.ts(T2, main="Hotelling test",pch=20,cex=0.8,xlab="sample index",ylab="T^2 score")
abline(h=Ta2, col="red")
dev.off()

jpeg(filename = "Q-statistics.jpg",
     width = 1200, height = 600, units = "px", pointsize = 18,
     quality = 100)
plot.ts(Q, main="Q-statistics",pch=20,cex=0.8, xlab="sample index", ylab="Q-score")
abline(h=Qa, col="red")
dev.off()

nQ <- log(Q/Qa)
nT2<- log(T2/Ta2)

Ylm <- max(abs(nQ),abs(nT2))

jpeg(filename = "anomaly comparison.jpg",
     width = 800, height = 600, units = "px", pointsize = 18,
     quality = 100)
plot(nQ,type="o", main="normalized anomaly score",pch='.',col='chartreuse',
     xlab="sample index", ylab="Normalized anomaly score",ylim=c(-Ylm,Ylm))
points(nT2,type="o", main="Hotelling test",col='coral',pch='.',xlab="sample index",
       ylab="T^2 score")
legend((n-12500),Ylm, c("Q-statistics","T squared"),
       lty=c(1,1),lwd=c(2.5,2.5),col=c('chartreuse','coral')) 
abline(h=0,col='blue')
dev.off()

# control charts plot with color coding clusters

jpeg(filename = paste0("Qstat_",NC,"_clusters.jpg"),
     width = 1200, height = 600, units = "px", pointsize = 18,
     quality = 100)
plot(Q, main="Q-statistic test",pch=19, col=Gtest_col,
     cex=0.8,xlab="sample index",ylab="Q-statistic score")
abline(h=Qa, col="black", lwd=2)
dev.off()

jpeg(filename = paste0("Hotelling_",NC,"_clusters.jpg"),
     width = 1200, height = 600, units = "px", pointsize = 18,
     quality = 100)
plot(T2, main="Hotelling test",pch=19, col=Gtest_col,
     cex=0.8,xlab="sample index",ylab="T^2 score")
abline(h=Ta2, col="black",lwd=2)
dev.off()

# plotting PC biplot with color coding for clusters

for (i in 2:k){
  for (j in 1:(i-1)){
    qplot(comp_PCA$scores[,c(i)],comp_PCA$scores[,c(j)],xlab=paste("Component",i),ylab=paste("Component",j),col = I(Gtest_col) ,asp=1)
    ggsave(paste0("_components_",i,"_",j,".jpg"), device = "jpeg")
  }
}

ggbiplot(comp_PCA, obs.scale = 1, var.scale = 1, ellipse = TRUE,var.axes = FALSE)
ggsave("biplot.jpg", device = "jpeg")

Xg <- matrix(data=1,nrow=m,ncol=1)
ggbiplot(comp_PCA, obs.scale = 1, var.scale = 1, ellipse = TRUE,
         groups = Xg,ellipse.prob = 0.955,
         var.axes=FALSE)
ggsave("biplcot.jpg", device = "jpeg")



# Distinguishing the clusters as normal operation and abnormal                       

PCAs <- list()
Anom_cluster <- list()
Anom <- list()
Anomt <- list()
QAnom <- list()
QAnomt <- list()
Tas2 <- matrix(data=0, nrow=NC, ncol=1)
Qas2 <- matrix(data=0, nrow=NC, ncol=1)
NNC <- c(1,2)         # normal cluster numbers
ANC <- c(3,4,5,6,7,8) # anomaly cluster numbers

# Developing control charts for normal cluster

for (i in 1:length(NNC)){
  
  clusp <- subset(X, Gtest_cluster == NNC[i])
  clust <- subset(gtime, Gtest_cluster == NNC[i])
  clusps <- scale(clusp)
  PCAs[[i]] <- princomp (clusp, scores=TRUE, cor=FALSE)
    
  m <- nrow(clusp)
  T2 <- matrix(ncol=1, nrow=m)
  Q <- matrix(ncol=1, nrow=m)
    
  jpeg(filename = paste0("PCs_cluster_",i,".jpg"),
       width = 500, height = 500, units = "px", pointsize = 16,
       quality = 100)
  plot(PCAs[[i]])
  dev.off()
    
  k <- 0
  c <- PCAs[[i]]$sdev[1]
    
  while (c >= 1){
      k <- k+1
      c <- PCAs[[i]]$sdev[k+1]
  }
  if (k==1){k <- 2}
    
  jpeg(filename = paste0("PCs_sdev_cluste_",i,".jpg"),
       width = 700, height = 600, units = "px", pointsize = 12,
       quality = 100)
  plot(PCAs[[i]]$sdev, pch=16,cex=1.6,col="blue", ylab="Standard deviation", xlab="PC count", cex.lab=1.4)
  abline(h=1, col="red", lwd=2)
  text(PCAs[[i]]$sdev[1:k], labels=round(PCAs[[i]]$sdev[1:k], digits = 3), cex= 1, pos=4)
  dev.off()
    
  V <- (PCAs[[i]]$sdev^2)
  P <- PCAs[[i]]$loadings[,1:k]
  Pt <- t(P)
  L <- diag(V[1:k])
  L1 <- solve(L)
  TR3 <- diag(d)
  F <- qf(a, df1=k, df2=m-k) 
  Tas2[i,] <- (k*((m^2)-1)/(m*(m-k)))*F
    
  teta <- matrix(data=0,ncol=1, nrow=3)
  for (t in 1:3){
    for (j in (k+1):d){
        teta[t] <- teta[t]+V[j]^i
    }
  }
    
  h0 <- 1-((2*(teta[1]*teta[3]))/(3*(teta[2]^2)))
  za <- qnorm (a, mean=0, sd=1)
  Qas2[i,] <- teta[1]*(((za*((2*teta[2]*(h0^2))^0.5)/teta[1])+1+((teta[2]*h0*(1-h0))/(teta[2]^2)))^(1/h0))
    
  for (j in 1:m){
    Bt <- clusps[j,]
    B <- t(Bt)
      
    #T2
    Y <- Pt %*% Bt
    Yt <- t(Y)
    TR <- L1 %*% Y
    T2[j] <- Yt %*% TR
    #Q
    TR2 <- P %*% Pt
    TR4 <- TR3 - TR2
    TR5 <- TR4 %*% Bt
    Q[j,] <- B %*% TR5
    }
    
  jpeg(filename = paste0("Hotelling_cluster_",i,".jpg"),
       width = 1200, height = 600, units = "px", pointsize = 18,
       quality = 100)
  plot.ts(T2, main="Hotelling test",pch=20,cex=0.8,xlab="sample index",ylab="T^2 score")
  abline(h=Tas2[i], col="red")
  dev.off()
    
  jpeg(filename = paste0("Q-statistics_cluster_",i,".jpg"),
       width = 1200, height = 600, units = "px", pointsize = 18,
       quality = 100)
  plot.ts(Q, main="Q-statistics",pch=20,cex=0.8, xlab="sample index", ylab="Q-score")
  abline(h=Qas2[i,], col="red")
  dev.off()
    
  Anom[[NNC[i]]] <- subset(clusp, T2 > Tas2[i])
  Anomt[[NNC[i]]] <- subset(clust, T2 > Tas2[i])
    
  QAnom[[NNC[i]]] <- subset(clusp, Q > Qas2[i])
  QAnomt[[NNC[i]]] <- subset(clust, Q > Qas2[i])
} 

# Adding the abnormal clusters to the anomaly values from normal clusters

for (i in 1:(length(ANC))){
  clusp <- subset(X, Gtest_cluster == ANC[i])
  clust <- subset(gtime, Gtest_cluster == ANC[i])
  clusps <- scale(clusp)
  
  Anom[[ANC[i]]] <- clusp
  Anomt[[ANC[i]]] <- clust
  
  QAnom[[ANC[i]]] <- clusp
  QAnomt[[ANC[i]]] <- clust
  
}

# generating biplots with the anomaly points color coded according to cluster

for (i in 2:d){
  for (j in 1:(i-1)){
    jpeg(filename = paste0("cluster_T2_anomalies_",j,"_",i,".jpg"),
         width = 800, height = 600, units = "px", pointsize = 16,
         quality = 100)
    plot(X[,c(j,i)],pch='.',col='gray60')
    for (z in 1:NC){
      points(Anom[[z]][,c(j,i)],pch=19,cex=1,col=pretty_palette[z])
    }
    points(fX[,c(j,i)],pch='+',cex=1.2,col='black')
    dev.off()
  }
}

# Heat-map biplots for BKDE distributions

NX <- subset(X,Gtest_cluster %in% NNC)



# BKDE distributions for every possible combination of two variables in dataframe
Y <- NX
n <- nrow(Y)
LI=list()
Lii <- data.frame(NULL)
k=1
for (i in 2:d){
  for (j in 1:(i-1)){
    h1=sd(Y[,j])/(n^(1/6))
    h2=sd(Y[,i])/(n^(1/6))
    M=cbind(Y[,j],Y[,i])
    L=bkde2D(M,bandwidth=c(h1,h2),gridsize=c(256,256))
    LI[[k]]=L
    t=interp.surface(list(x=L$x1,y=L$x2,z=L$fhat),loc=M)
    tstar=max(L$fhat)
    AST=-log(t/tstar)
    Lii <- rbind(Lii,new=c(i,j,k))
    k=k+1
  } 
}
k=1
Lii <- as.matrix(Lii)
colnames(Lii) <- c('i','j','k')


k=1
for (i in 2:d){
  for (j in 1:(i-1)){
    L=LI[[k]]  
    tstar=max(L$fhat)
    jpeg(filename = paste0("Kernel_map_test_",i,"_",j,".jpg"),
         width = 800, height = 600, units = "px", pointsize = 16,
         quality = 100)
    contour(x=L$x1,y=L$x2,z=-log((L$fhat)/tstar),
            levels=seq(0,6,0.1),labcex = 1.8,lwd = 8,
            drawlabels = FALSE,
            vfont = c("sans serif", "bold"), 
            col=rgb(1,1-(seq(0,6,0.1)/6),0,alpha=0.2),
            xlab=colnames(X[j]),ylab=colnames(X[i])) 
    #points(Y[,j],Y[,i],pch='.',cex=0.5,col='black')
    dev.off()
    k=k+1
  }
}


# Heat-maps with anomaly plot combination based on T2
k=1

for (i in 2:d){
  for (j in 1:(i-1)){
    
    L=LI[[k]]  
    tstar=max(L$fhat)
    jpeg(filename = paste0("Kernel_map_all_",i,"_",j,".jpg"),
         width = 1200, height = 800, units = "px", pointsize = 16,
         quality = 100)
    contour(x=L$x1,y=L$x2,z=-log((L$fhat)/tstar),
            levels=seq(0,6,0.1),labcex = 1.8,lwd = 8,
            drawlabels = FALSE,
            vfont = c("sans serif", "bold"), 
            col=rgb(1,1-(seq(0,6,0.1)/6),0,alpha=0.1),
            xlab=colnames(X[j]),ylab=colnames(X[i])) 
    points(X[,j],X[,i],pch='.',cex=0.5,col='gray60')
    for (z in 1:NC){
      points(Anom[[z]][,c(j,i)],pch=20,cex=1.6,col=pretty_palette[z])
    }
    points(fX[,c(j,i)],pch='+',cex=1,col='black')
    dev.off()
    k=k+1
  }
}


# Heat-maps with anomaly plot combination based on Q statistics
k=1

for (i in 2:d){
  for (j in 1:(i-1)){
    
    L=LI[[k]]  
    tstar=max(L$fhat)
    jpeg(filename = paste0("Kernel_map_all_",i,"_",j,".jpg"),
         width = 1200, height = 800, units = "px", pointsize = 16,
         quality = 100)
    contour(x=L$x1,y=L$x2,z=-log((L$fhat)/tstar),
            levels=seq(0,6,0.1),labcex = 1.8,lwd = 8,
            drawlabels = FALSE,
            vfont = c("sans serif", "bold"), 
            col=rgb(1,1-(seq(0,6,0.1)/6),0,alpha=0.1),
            xlab=colnames(X[j]),ylab=colnames(X[i])) 
    points(X[,j],X[,i],pch='.',cex=0.5,col='gray60')
    for (z in 1:NC){
      points(Anom[[z]][,c(j,i)],pch=20,cex=1.6,col=pretty_palette[z])
    }
    points(fX[,c(j,i)],pch='+',cex=1,col='black')
    dev.off()
    k=k+1
  }
}
k=1

