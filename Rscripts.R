drawChordDiagram <- function (chordfile = "chord.lists.txt", ratio = 1, cex = 1)
{
	require( circlize )
	circos.clear()
	circos.par(start.degree = 270, message = F, points.overflow.warning = F)
	g = read.table(chordfile ,comment.char = "", check.names = FALSE,sep="\t")
	chord = g[,c(1:3,3)]
	chord[,4] = chord[,3] / 3 * ratio
	
	g1 = unique(g[,c(1,4)])
	node_color1 = as.character(g1[,2])
	names(node_color1) = as.character(g1[,1])

	g2 = unique(g[,2])
	node_color2 = rep("grey",length(g2))
	names(node_color2) = as.character(g2)
	node_color = c(node_color1, node_color2)

	gene_res = chordDiagram( chord ,annotationTrack = c("grid"),preAllocateTracks = 1,grid.border=T,grid.col=node_color,col=as.character(g[,5]),directional = TRUE,)
	circos.track(track.index = 1, panel.fun = function(x, y) { circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", cex=cex, niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA)
	for( i in seq_len(nrow(gene_res))){
		circos.rect(gene_res[i,"x1"], -uy(2,"mm"), gene_res[i,"x1"] - abs(gene_res[i, "value1"]), -uy(19, "mm"), col = as.character(g[i,6]), border=T,lwd=.1, sector.index = gene_res$rn[i], track.index=2)
	}
}


getChordDiagram <- function ( module, ratio = 1, tfnum = 50, targetnum = 15, cex = 1)
{
	if (!file.exists("results.regulator.tfs.txt")){
		cat ("The file 'results.regulator.tfs.txt' not found.\n")
		return("Exit")
	}

	if (!file.exists("getChordLists.pl")){
		cat ("The file 'getChordLists.pl' was not found. This file is required for extracting the TF-target gene pairs for the input module. It should be placed within the home folder of the EXPLICIT package.\n")
		return("Exit")
	}
	a <- read.table("results.regulator.tfs.txt",heade=T,sep="\t")
	if ( module %in% a[,1]){
		unlink( "chord.lists.txt" ) 
		perlcmd = paste("perl getChordLists.pl",module,tfnum,targetnum,sep=" ")
		system( perlcmd )
		if (!file.exists("chord.lists.txt")){
			cat ("Perl is required for the analysis. Make sure Perl is installed properly and it can be invoked from command line.\n")
			return("Exit")
		}
		drawChordDiagram( chordfile = "chord.lists.txt", ratio = ratio, cex = cex)
	}else{
		cat ("The module ",module," was not found in 'results.regulator.tfs.txt'.","\n",sep="")
	}
}


