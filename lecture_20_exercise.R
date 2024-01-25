
################################################################################

library(NbClust)
library(factoextra)
library(cluster)

### Acute lymphoblastic leukemia (ALL) data
data <- read.table( "ALL_data.data" )

dat <- data[,1:(ncol(data)-1)]
# the gene expression

cl <- data$cl
# the true cell type

dim( dat )
# more features than samples

# create new data set with only the 100 genes with highest variation
genes.var <- apply( dat , 2 , var )
# calcola la var per ogni colonna (dim2)
# variation of genes, see ?apply
# how much variation do we have for each gene
genes.var.select <- order(-genes.var)[1:100]
# see ?order
# ordina in ordine crescente ed estrae gli indici
dat.s <- dat[,genes.var.select]
dim(dat.s)

# Do cluster analysis on the complete and reduced data sets. 
# Can you recover the true cell types? 
# Try both hierarchical clustering and k-means. Which k is optimal?
# Try also k-medoids. Which k is optimal?


################################################################################

fviz_nbclust( dat , kmeans , method = "silhouette" , k.max = 15 )
# k_opt = 3

fviz_nbclust( dat , kmeans , method = "wss" , k.max = 15 )

k_dat <- kmeans( dat , centers = 3 , nstart = 1 )

x.dat <- table( k_dat$cluster , cl )
x.dat

k_dat$betweenss / k_dat$totss

sil.dat <- silhouette( k_dat$cluster , dist( dat ) )  # 0.1
fviz_silhouette( sil.dat )


fviz_cluster( k_dat , dat , ellipse.type = "norm" )



################################################################################


k_dat.s <- kmeans( dat.s , centers = 2 , nstart = 1 )

x.dat.s <- table( k_dat.s$cluster , cl )
x.dat.s

k_dat.s$betweenss / k_dat.s$totss


sil.dat.s <- silhouette( k_dat.s$cluster , dist( dat.s ) )  # 0.32
fviz_silhouette( sil.dat.s )


fviz_nbclust( dat.s , kmeans , method = "silhouette" , k.max = 15 )
# k_opt = 2

fviz_nbclust( dat.s , kmeans , method = "wss" , k.max = 15 )

fviz_cluster( k_dat.s , dat.s , ellipse.type = "norm" )




################################################################################


dist_mat <- dist( dat )
# image( as.matrix( dist_mat ) )

hclust_avg <- hclust( dist_mat , method = 'complete' )
plot( hclust_avg , labels = cl )

rect.hclust( hclust_avg , k = 2 , border = 2:6 )


cut_avg <- cutree( hclust_avg , k = 2 )
#cut_avg
table( cut_avg , cl )



################################################################################


dist_mat.s <- dist( dat.s )
image( as.matrix( dist_mat.s ) )

hclust_avg.s <- hclust( dist_mat.s , method = 'complete' )
plot( hclust_avg.s )

rect.hclust( hclust_avg.s , k = 2 , border = 2:6 )


cut_avg.s <- cutree( hclust_avg.s , k = 2 )
#cut_avg
table( cut_avg.s , cl )

# why looking genes with higher variation we are able to predict different cell
# types? because there are genes expressed in all cells and some others
# are specialized only in some cell types



################################################################################


kmed <- pam( dat , k = 2 )

table( kmed$cluster , cl )

# molto meglio di k-means !!



################################################################################


kmed.s <- pam( dat.s , k = 2 )

table( kmed.s$cluster , cl )

fviz_nbclust( dat, cluster::pam , method = "silhouette" )

# il migioramento che ho nella predizione usando data.s non è dovuto alla riduzione
# della dimensionalità dei dati

