#This is an edited version of a co-phyloplot script from Paola De Lima Ferreira


library(phytools)

setwd("C:/Users/kristine")
dir()


#Import Data - trees

#paralogs
tr1 <- read.tree("astral_tree_probabilities_renamed.tre")
plot(tr1, show.node.label = TRUE, use.edge.length = F)
title("With paralogs")
nodelabels()

tr1reroot <- reroot(tr1, node.number=98, position=NULL, interactive=FALSE)
plot(tr1reroot, show.node.label = TRUE, use.edge.length = F)
str(tr1reroot)
tr1reroot$tip.label
tr1reroot$edge.length[is.na(tr1reroot$edge.length)] <- 1
tr1reroot$edge.length <- rep(1, length(tr1reroot$edge.length)) 
tr1reroot$edge.length[is.na(tr1reroot$edge.length)] <- 1
write.tree(tr1reroot, file = "rrPhylogeny_paralogs.tre")

#no_paralogs
tr2 <- read.tree("Phylogeny_ceroxyloideae.nex")
plot(tr2, show.node.label = TRUE, use.edge.length = F)
nodelabels()


tr2reroot <- reroot(tr2, node.number=152, position=NULL, interactive=FALSE)
plot(tr2reroot, show.node.label = TRUE, use.edge.length = F)
title("Without paralogs")
str(tr2reroot)
tr2reroot$tip.label
tr2reroot$edge.length[is.na(tr2reroot$edge.length)] <- 1
tr2reroot$edge.length <- rep(1, length(tr2reroot$edge.length)) 
tr2reroot$edge.length[is.na(tr2reroot$edge.length)] <- 1
write.tree(tr2reroot, file = "rrPhylogeny_no_paralogs.tre")

#Check the tip names
tr1reroot$tip.label
tr2reroot$tip.label
match(tr1reroot$tip.label, tr2reroot$tip.label)

#Co-phyogenies

obj<-cophylo(tr1reroot,tr2reroot, use.edge.length	= FALSE)
pdf("2f_paralogs_vs_no_paralogs.pdf", height=60, width=40)
plot(obj,link.type="curved",link.lwd=1,link.lty="solid",cex=0.1, link.col=make.transparent("black", 0.8), fsize =2.8, label.offset = 0.75, use.edge.lenth = FALSE, link.col = (3))
nodelabels.cophylo(obj$trees[[1]]$node.label,1:obj$trees[[1]]$Nnode+
                     Ntip(obj$trees[[1]]), cex=2.5)
title("Multi-copy genes not separated", adj=0.1, line=-5, cex.main=5)
nodelabels.cophylo(obj$trees[[2]]$node.label,1:obj$trees[[2]]$Nnode+
                     Ntip(obj$trees[[2]]),which="right", cex=2.5)
title("Multi-copy genes separated", adj=0.85, line=-5, cex.main=5)



dev.off()

