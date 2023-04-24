# David C. Green
# 2023-03-21
# fungal microbiome analysis
# ---------------------------------------------------------------------
# Libraries
library(tidyverse)
library(ggthemes)
library(patchwork)
library(igraph)
library(plotly)
library(rethinking)
library(beepr)

# ---------------------------------------------------------------------
# Sources
apr2018OTU <- read.csv('Apr_2018_OTU_table_filtered_genus_minus_unidentified.csv')
apr2018tax <- read.csv('Apr_2018_taxonomy_sorted_genus_minus_unidentified.csv')
sep2018OTU <- read.csv('Sep_2018_OTU_table_genus_minus_unidentified.csv')
sep2018tax <- read.csv('Sep_2018_taxonomy_genus_minus_unidentified.csv')
source('fmb_functions.R')
# ---------------------------------------------------------------------



# inducing a graph from data for april data
#A1
A1_apr_18 <- data.frame(apr2018OTU$X.NAME, apr2018OTU$A1)
g <- graph_from_data_frame(A1_apr_18)
lab <- seq(1, 481)

plot(g, 
     mark.shape = 0.01,
     xlim = c(-0.6,0.6),
     ylim = c(-0.6, 0.6),
     vertex.size = 0.01,
     edge.width = 0.2,
     edge.arrow.size = 0.2,
     vertex.label.cex = 0.4,
     vertex.label.color = 'black',
     vertex.color = 'grey')


r_list <- NULL
for(i in 1:nrow(apr2018OTU)) {
  a <- i
  b <- i + 1
  r_list[i] <- predTest(a,b)
}
tfR <- r_list > 0.6
dfR <- data.frame(r_list, tfR)


fit1 <- QUADfit(1, 2, apr2018OTU)
df1 <- predTest(1, 2)
xseq <- seq( from=min(apr2018OTU[1,])-0.15 , to=max(apr2018OTU[2,])+0.15 , length.out=30 )
mu <- link(fit1 , data=list(cor=xseq) )
mu_mean <- apply(mu,2,mean)
plot( apr2018OTU[1,] ~ apr2018OTU[2,] , data=apr2018OTU)
lines( xseq , mu_mean , lwd=2 )


#create correlation matrix
corMat <- correlationMatrix(apr2018OTU)
#set blank matrix
corMatTF <- matrix(nrow = nrow(apr2018OTU), ncol = nrow(apr))
#set r value threshold
r <- corMat > 0.8 & corMat < 1
#create TF 0/1 matrix
for(i in 1:481) {
  for(j in 1:481) {
    if(r[i,j] == T) {corMatTF[i,j] <- 1} else {corMatTF[i,j] <- 0}
  }
}
write(corMatTF, file = 'corrMatrix')

# make genus list
gen <- apr2018tax$Genus
# list of all genera
genL <- distinct(apr2018tax, across(contains('Genus')))
# genus list as df
genA <- data.frame(gen)

# vertex colors function
#create gradient
vgradient <- c('white', 'purple')
#generate color list
VgenCol <- vcolors(genL, gen, vgradient)
    
means <- avgAbundance(apr2018OTU)

meanAbCols <- vcolors2(means, vgradient)



g <- graph_from_adjacency_matrix(corMatTF)  
lay <- layout_(g, on_grid(width = 20, height = 0))
# lay <- layout_(g, with_graphopt())
deg <- degree(g, v = V(g), mode = 'out')
vertex_attr(g) <- list(color = meanAbCols,
                       names = gen,
                       degree = deg)
Vattr <- vertex.attributes(g)
VattrDF <- data.frame(name = Vattr$names, color = Vattr$color, degree = Vattr$degree)

plot(g,
     layout = lay,
     mark.shape = 0.01,
     edge.width = 0.2,
     edge.arrow.size = 0.01,
     vertex.label.cex = 0.4,
     vertex.label.color = 'black',
     vertex.size = Vattr$degree*3)

gheatmap <- as_adjacency_matrix(g, sparse=F)
colnames(gheatmap) <- V(g)$names
rownames(gheatmap) <- V(g)$names

palf <- colorRampPalette(c('white', 'purple')) 
heatmap(gheatmap[,17:1], Rowv = NA, Colv = NA, col = palf(100), 
        scale="none", margins=c(10,10) )



iso <- which(degree(g) == 0)
g2 <- delete_vertices(g, iso)
lay2 <- lay[-iso,]
ns <- c(10,10,10)
deg2 <- deg
for(i in 1:length(deg2)) {
  ifelse(deg2[i] == 0, deg2[i] <- NA, next)
}
degdf <- data.frame(deg2)
degdf <- na.omit(degdf)
deg2 <- degdf$deg2

V2attrDF <- filter(VattrDF, VattrDF$degree > 0)

plot(g2,
     layout = lay2,
     mark.shape = 1,
     edge.width = 0.2,
     edge.arrow.size = 0.1,
     vertex.label.cex = deg2*0.02,
     vertex.label.color = 'black',
     vertex.size = V2attrDF[,3]*2,
     vertex.label =  V2attrDF$name)

degg2 <- degree(g2)
tmax <- centr_degree_tmax(g2, loop = FALSE)
centralize(degg2, tmax)


legend_image <- as.raster(matrix(collist(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)



df <- data.frame(1:148, 1:148)
plot(df, col = VgenCol)

ego1 <- make_ego_graph(g2, order = 1, nodes = V(g2), mode = 'out')
g227 <- ego1[[227]]
V227attr <- vertex_attr(g227)
lay3 <- layout_(g227, as_tree())
plot(g227,
     layout = lay3,
     mark.shape = 1,
     edge.width = 0.2,
     edge.arrow.size = 0.1,
     vertex.label.color = 'black',
     vertex.label =  V227attr$names,
     vertex.label.cex = 0.5,
     main = 'Asscoiated Genera')


g235 <- ego1[[235]]
V235attr <- vertex_attr(g235)
lay3 <- layout_(g235, as_tree())
plot(g235,
     layout = lay3,
     mark.shape = 1,
     edge.width = 0.2,
     edge.arrow.size = 0.1,
     vertex.label.color = 'black',
     vertex.label =  V235attr$names,
     vertex.label.cex = 0.5,
     main = 'Asscoiated Genera')

ego2 <- make_ego_graph(g2, order = 2, nodes = V(g2), mode = 'out')
g235_2 <- ego2[[235]]
V235attr <- vertex_attr(g235_2)
lay3 <- layout_(g235_2, as_tree())
plot(g235_2,
     layout = lay3,
     mark.shape = 1,
     edge.width = 0.2,
     edge.arrow.size = 0.1,
     vertex.label.color = 'black',
     vertex.label =  V235attr$names,
     vertex.label.cex = 0.5,
     main = 'Asscoiated Genera')









# September data

corMat2 <- correlationMatrix(sep2018OTU)
#set blank matrix
corMatTF2 <- matrix(nrow = nrow(sep2018OTU), ncol = nrow(sep2018OTU))
#set r value threshold
r <- corMat2 > 0.8 & corMat2 < 1
#create TF 0/1 matrix
for(i in 1:481) {
  for(j in 1:481) {
    if(r[i,j] == T) {corMatTF2[i,j] <- 1} else {corMatTF2[i,j] <- 0}
  }
}

Fcount <- 0
Tcount <- 0
for(i in 1:nrow(sep2018OTU)) {
  for(j in 1:length(sep2018OTU)) {
    ifelse(is.na(sep2018OTU[i,j]) == TRUE, Tcount <- Tcount + 1, Fcount <- Fcount + 1)
  }
}

for(i in 1:length(sep2018OTU)) {
  ifelse(is.na(sep2018OTU[50,i]) == TRUE, Tcount <- Tcount + 1, Fcount <- Fcount + 1)
}
