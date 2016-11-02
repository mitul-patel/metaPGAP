> library("gplots")
> library("ape")
> options(expressions = 100000) #https://stat.ethz.ch/pipermail/r-help/2004-January/044109.html
> tab <- read.table(file="intersection/pangenome_matrix_t0.tab", header=TRUE)
> mat_dat <- data.matrix(tab[,2:ncol(tab)])
> rnames <- tab[,1]
> rownames(mat_dat) <- rnames
> 
> if(100 < 100)
+ {
+    tmp_mat = mat_dat
+    diag(tmp_mat) = NA
+    rows <- (!apply( tmp_mat , 1 , function(x) any( x > 100 , na.rm=T) ) )
+    mat_dat <- mat_dat[rows, rows]
+ }
> 
> if(1 > 0){
+   svg("intersection/pangenome_matrix_t0_heatmap.svg", width=15, height=10, pointsize=)  
+   heatmap.2(mat_dat, cellnote=mat_dat, main="intersection/pangenome_matrix_t0.tab", notecol="black", density.info="none", trace="none", margins=c(18,18), lhei = c(1,5))
+   dev.off()
+ } else {
+   svg("intersection/pangenome_matrix_t0_heatmap.svg", width=15, height=10, pointsize=)  
+   heatmap.2(mat_dat, cellnote=mat_dat, main="intersection/pangenome_matrix_t0.tab", notecol="black", density.info="none", trace="none", margins=c(18,18), lhei = c(1,5), dendrogram = "row", Colv = FALSE)
+    dev.off()
+ }  
null device 
          1 
> if(0 > 0){
+   sim2dist <- function(x) 100 - x
+   bionj <- bionj(as.dist(apply(mat_dat, 1, sim2dist)))
+   write.tree(phy=bionj, file="intersection/pangenome_matrix_t0_BioNJ.ph")
+ }
> 
