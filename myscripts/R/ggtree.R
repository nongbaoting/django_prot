knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(ggtree)
library(treeio)
#input
# wd = "/dat1/nbt2/proj/21-prot/web/data/res/phylotree/b246ad08-86c8-4719-9802-3585a0436803"
# setwd(wd)
tree = read.newick("seq.msa.treefile")
id_name = read.table("id_names.txt", sep="\t")
names(id_name) = c("label","desc")

x = as_tibble(tree)
#join
left_join(x, id_name)  %>% unique(.) ->x
names(x)[c(4,5)] = c("desc", "label")
tree= as.phylo(x)
ggtree(tree, branch.length='none')+ geom_tiplab()  + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 600, 6, 6)) ->p

pdf("ggtree.pdf", width=12,height = 12)
#plot
p
dev.off()
