library(igraph)
library(ggraph)
library(tidygraph)
library(qgraph)
library(cowplot)

layout=matrix(c(1,3,
                0,3,
                1,4,
                1,1,
                0,1,
                1,0,
                3,3,
                3,4,
                4,3,
                3,1,
                4,1,
                3,0), byrow = T, ncol=2)


memberships=c(rep("A",3), rep("B",3), rep("C",3), rep("D",3))
memberships


expand.grid(memberships, memberships)

memb.mat=outer(memberships, memberships, "==")
memb.mat

mat1=apply(memb.mat, c(1,2), function(x) if(x==TRUE) x=0.3 else x=0)

g1=graph_from_adjacency_matrix(mat1, mode="undirected", weighted=T, diag=F)
plot(g1, layout=layout, edge.width=E(g1)$weight*10, edge.color="black", vertex.label="", vertex.color="gray")

p.in=0.1
p.out=0.1

mat2=apply(memb.mat, c(1,2), function(x) if(x==TRUE) x=rnorm(1,mean=p.in,sd=0.1) else x=rnorm(1,mean=p.out,sd=0.1))
mat2[mat2<0.2]=0
mat2=(mat2+t(mat2))/2
g2=graph_from_adjacency_matrix(mat2, mode="undirected", weighted=T, diag=F)
#plot(g2, layout=layout, edge.width=E(g2)$weight*10, edge.color="black", vertex.label="", vertex.color="gray")


p.in=0.3
p.out=0.3
mat3=apply(memb.mat, c(1,2), function(x) if(x==TRUE) x=rnorm(1,mean=p.in,sd=0.5) else x=rnorm(1,mean=p.out,sd=0.5))
mat3[mat3<0.2]=0
mat3=(mat3+t(mat3))/2
g3=graph_from_adjacency_matrix(mat3, mode="undirected", weighted=T, diag=F)
#plot(g3, layout=layout, edge.width=E(g3)$weight*10, edge.color="black", vertex.label="", vertex.color="gray")

p.in=0.6
p.out=0.1
mat4=apply(memb.mat, c(1,2), function(x) if(x==TRUE) x=rnorm(1,mean=p.in,sd=0.1) else x=rnorm(1,mean=p.out,sd=0.1))
mat4[mat4<0.2]=0
mat4=(mat4+t(mat4))/2
g4=graph_from_adjacency_matrix(mat4, mode="undirected", weighted=T, diag=F)
#plot(g4, layout=layout, edge.width=E(g4)$weight*10, edge.color="black", vertex.label="", vertex.color="gray")


shps=c("triangle", "triangle", "triangle", "circle", "circle", "circle", "square", "square", "square", "diamond", "diamond", "diamond")
png("figs/Conceptual.png",width=13,height=5,units="in",res=300)
par(mfrow=c(1,3))
#p1=qgraph(mat1, color="gray", diag=F,fade=F, layout=layout, shape=shps, edge.color="black")
qgraph(mat2, color="gray", diag=F,fade=F, layout=layout, shape=shps, edge.color="black", maximum=1)
qgraph(mat3, color="gray", diag=F,fade=F, layout=layout, shape=shps, edge.color="black", maximum=1)
qgraph(mat4, color="gray", diag=F,fade=F, layout=layout, shape=shps, edge.color="black", maximum=1)

dev.off()

