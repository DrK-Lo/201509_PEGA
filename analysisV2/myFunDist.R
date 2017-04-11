  myFun <- function(inDist) {
    if (class(inDist) != "dist") stop("wrong input type")
    A <- attr(inDist, "Size")
    B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
    if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
    if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
    data.frame(
      row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
      col = rep(B[-length(B)], (length(B)-1):1),
      value = as.vector(inDist))
  }
  

MakeCorNetwork <- function(clustnum, filename, corthreshold, colordf){
  
  ### Based on correlation in associations (co-association)
  names(results_pine3)[env_cols]
  multidat <- abs(results_pine3[which(results_pine3$cluster==clustnum), env_cols])
  #str(multidat)
  rownames(multidat) <- results_pine3$gtcontig_pos[which(results_pine3$cluster==clustnum)]
  dim(multidat)
  
  tidy_cors <- t(multidat) %>% 
    correlate() %>% 
    stretch()
  #head(tidy_cors)
  
  nedges <- nrow(tidy_cors %>%
    filter(abs(r) > corthreshold))
  
  print(c("Number of edges:", nedges))
  
  graph_cors <- tidy_cors %>%
    filter(abs(r) > corthreshold) %>%
    graph_from_data_frame(directed = FALSE )
  
    #write_graph(tidy_cors %>%  filter(abs(r) > 0.75) %>% graph_from_data_frame(directed = FALSE ), file = paste("../results/NetworkCorr0",filename,".GML",sep=""), format="gml")
    
    #write.table(tidy_cors %>%  filter(abs(r) > 0.75) , file = paste("../results/NetworkCorr0",filename,".txt",sep=""), row.names=FALSE)
  
    # check node names
    V(graph_cors)$name
    toadd <- !(rownames(multidat) %in% V(graph_cors)$name)
    graph_cors2 <- add_vertices(graph_cors, nv = sum(toadd, na.rm=TRUE), name= rownames(multidat)[toadd])
    V(graph_cors2)$name
    sum(!(rownames(multidat) %in% V(graph_cors2)$name))
    graph_cors <- graph_cors2
    
    # line up names
    (myo <- match(V(graph_cors)$name, results_pine3$gtcontig_pos))
    
    contigname = results_pine3$gtcontig[myo]
    
    myocol <- colordf$col[match(contigname, colordf$gtcontig)]
    myoname <- colordf$contigID[match(contigname, colordf$gtcontig)]
    
    #(myocol <- results_pine3$clustercol[myo])
    #V(graph_cors)$name <- gsub("_", "\n",V(graph_cors)$name )
    #z <- cbind(myocol, V(graph_cors)$name, results_pine3$gtcontig_pos[myo])
    
   #print(table(z[,1])) #check against colors
  
   png(file = paste("../results/NetworkCorr2Hist", filename,".pdf", sep=""), width=15, height=15, units="in", res=300)
   hist(tidy_cors$r)
   dev.off()
    #head(E(graph_cors)$r)
#pdf(file = paste("../results/NetworkCorr", filename,".pdf", sep=""), width=15, height=15) #units="in", res=300)
    ggraph(graph_cors, layout="kk") +
    geom_edge_fan2(aes(edge_alpha = 0.01)) +
    guides(edge_alpha = "none", edge_width = "none") +
    #scale_edge_colour_gradientn(limits = c(-1, 1), 
    #                            colors = c("firebrick2", "dodgerblue2")) +
    geom_node_point(color = myocol, size = 10) +
    geom_node_text(aes(label = myoname, size=0.1)) +
    theme_graph() +
    labs(title = filename)
   
   ggsave(paste("../results/Network", filename,corthreshold, "corr.pdf", sep="_"), width=15, height=15, units="in")
   
    #dev.off()
} # end function