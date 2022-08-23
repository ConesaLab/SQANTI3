generatePDFreport = function() 
  {
  pdf(file=pdf.report.file, width = 7.5, height = 6.5)

  #cover
  grid.newpage()
  cover <- textGrob("SQANTI3 report",
                    gp=gpar(fontface="italic", fontsize=40, col="orangered"))
  grid.draw(cover)
  
  # TABLE 1: Number of isoforms in each structural category
  
  freqCat <- as.data.frame(table(data.class$structural_category))
  #freqCat$ranking = order(freqCat$Freq,decreasing = T)
  table1 <- tableGrob(freqCat, rows = NULL, cols = c("Category","Isoforms, count"))
  title1 <- textGrob("Transcript Classification\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
  gt1 <- gTree(children=gList(table1, title1))
  
  
  # TABLE 2: Number of Novel vs Known Genes
  freqCat = as.data.frame(table(isoPerGene$novelGene))
  table2 <- tableGrob(freqCat, rows = NULL, cols = c("Category","Genes, count"))
  title2 <- textGrob("Gene Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -4)
  gt2 <- gTree(children=gList(table2, title2))
  
  
  # TABLE 3: Junction Classification
  
  uniq_sj_count <- nrow(uniqJunc)
  
  freqCat <- as.data.frame(table(uniqJunc$SJ_type))
  freqCat$Var1 <- gsub(" ", "", freqCat$Var1)
  freqCat$Var1 <- gsub("\n", " ", freqCat$Var1)
  freqCat$Frac <- round(freqCat$Freq*100 / uniq_sj_count, 2)
  table2 <- tableGrob(freqCat, rows = NULL, cols = c("Category","SJs, count","Percent"))
  title2 <- textGrob("Splice Junction Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
  gt3 <- gTree(children=gList(table2, title2))
  
  
  # TABLE 4: Summary number of Unique Isoforms and Unique Genes
  nGenes = nrow(isoPerGene)
  nIso = nrow(data.class)
  sn = paste("Unique Genes: ", nGenes, "\n", "Unique Isoforms: ", nIso)
  gt4 <- textGrob(sn, gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  
  
  # Plot Table 1 and Table 2
  grid.arrange(gt4,gt2,gt1, layout_matrix = cbind(c(1,2),c(1,4)))
  grid.arrange(gt3)
  
  s <- textGrob("Gene Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p0)
  print(p7)
  print(p6)
  print(p.classByLen.a)
  print(p.classByLen.b)
  
  if (!all(is.na(data.class$iso_exp))){
    print(p10)
  }
  if (!all(is.na(data.class$FL))){
    print(p11)
  }
  
  # PLOT length of isoforms
  # p.length.all: length of all isoforms, regardless of category
  # p.length.cat: length of isoforms, by category
  # p.length.exon: length of isoforms, mono- vs mult-exon/ufrc/conesa/fpardopalacios/SQANTI_QDE/SQANTI3/melanoma_example/melanoma_chr13_tappAS_annot_from_SQANTI3.gff3
  # (optional) p.length.all.sample: length of all isoforms by sample
  print(p.length.all)
  print(p.length.cat)
  print(p.length.exon)
  if (length(FL_multisample_indices)>0) {
    print(p.length.all.sample)
    print(p.length.exon.sample)
  }
  
  # 2. general parameters by structual categories
  s <- textGrob("Structural Isoform Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p1)
  if (length(p1.s.list) > 0) {
    for (i in 1:length(p1.s.list)) {
      print(p1.s.list[i])
    }
  }
  print(p4)
  print(p4.s1)
  print(p4.s2)
  print(p4.s3)
  print(p5)
  print(p5.s1)
  print(p5.s2)
  print(p5.s3)
  print(pSTM)
  print(pSTM_perc)
  print(pSTM.s1)
  print(pSTM_perc.s1)
  print(pSTM.s2)
  print(pSTM_perc.s2)

  if (!all(is.na(data.class$iso_exp))){
    print(p8)
    if (!all(is.na(data.FSMISM$iso_exp))){
      print(p8.s1)
    }
    if (!all(is.na(data.NICNNC$iso_exp))){
      print(p8.s2)
    }
    if (!all(is.na(data.other$iso_exp))){
      print(p8.s3)
    }
  }
  if (!all(is.na(data.class$FL))){
    print(p9)
    if (!all(is.na(data.FSMISM$FL))){
      print(p9.s1)
    }
    if (!all(is.na(data.NICNNC$FL))){
      print(p9.s2)
    }
    if (!all(is.na(data.other$FL))){
      print(p9.s3)
    }
  }
  
  
  # (optional) p.FL_TPM_sample.by_cat
  # (optional) p.FL_TMP_sample.by_length
  #if (length(FL_multisample_indices)>0 & length(FL_multisample_indices)<4 ) {
  if (length(FL_multisample_indices)>0 & length(FL_multisample_indices)<0 ) {
    data.class$length_cat <- "<1kb"
    data.class[data.class$length>=1000,"length_cat"] <- "1-3kb"
    data.class[data.class$length>=3000&data.class$length<5000,"length_cat"] <- "3-5kb"
    data.class[data.class$length>=5000&data.class$length<10000,"length_cat"] <- "5-10kb"
    data.class[data.class$length>10000,"length_cat"] <- ">10kb"
    
    for (i in 1:(length(FL_TPM_multisample_names)-1)) {
      j1 <- FL_TPM_multisample_names[i]
      j1_log10 <- paste(FL_TPM_multisample_names[i], "_log10", sep='')
      n1 <- FL_multisample_names[i]
      for (i2 in (i+1):length(FL_multisample_names)) {
        j2 <- FL_TPM_multisample_names[i2]
        j2_log10 <- paste(FL_TPM_multisample_names[i2], "_log10", sep='')
        n2 <- FL_multisample_names[i2]
        
        print(paste("Printing FL TPM for sample", j1, "vs", j2, "...."))
        
        max_j1j2 <- floor(max(data.class[,j1_log10], data.class[,j2_log10])) + 1
        pearson <- round(cor(data.class[,j1], data.class[,j2], method="pearson"), 2)
        p.tmp <- ggplot(data.class, aes_string(j1_log10, j2_log10, color="structural_category")) +
          geom_point(alpha=0.3) +
          annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
          xlim(c(0, max_j1j2)) +
          ylim(c(0, max_j1j2)) +
          labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2)) +
          guides(fill=FALSE) +
          mytheme
        print(p.tmp)
        ggsave(paste("Rplot.",n1,"vs",n2,".by_cat.png",sep=''), p.tmp, width=8, height=6)
        
        data.class.gene <- group_by(data.class, by=associated_gene) %>% dplyr::summarise(n=dplyr::n(), sum1=sum(!!sym(j1)), sum2=sum(!!sym(j2)))
        pearson <- round(cor(data.class.gene$sum1, data.class.gene$sum2, method="pearson"), 2)
        p.tmp.gene <- ggplot(data.class.gene, aes(x=log10(sum1), y=log10(sum2))) +
          geom_point(alpha=0.3, color='orange') +
          annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
          xlim(c(0, max_j1j2)) +
          ylim(c(0, max_j1j2)) +
          xlab(j1) +
          ylab(j2) +
          labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, ", grouped by gene")) +
          guides(fill=FALSE) +
          mytheme
        print(p.tmp.gene)
        ggsave(paste("Rplot.",n1,"vs",n2,".summed_by_gene.png",sep=''), p.tmp.gene, width=8, height=6)
        
        
        data.class.tmp <- subset(data.class,length_cat=="<1kb")
        pearson <- round(cor(data.class.tmp[,j1], data.class.tmp[,j2], method="pearson"), 2)
        p.tmp.le1k <- ggplot(data.class.tmp, aes_string(j1_log10, j2_log10)) +
          geom_point(alpha=0.3, color='orange') +
          annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
          xlim(c(0, max_j1j2)) +
          ylim(c(0, max_j1j2)) +
          labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, "< 1kb only"))+
          guides(fill=FALSE) +
          mytheme
        print(p.tmp.le1k)
        ggsave(paste("Rplot.",n1,"vs",n2,".le1k.png",sep=''), p.tmp.le1k, width=8, height=6)
        
        data.class.tmp <- subset(data.class,length_cat=="1-3kb")
        pearson <- round(cor(data.class.tmp[,j1], data.class.tmp[,j2], method="pearson"), 2)
        p.tmp.1to3k <- ggplot(data.class.tmp, aes_string(j1_log10, j2_log10)) +
          geom_point(alpha=0.3, color='purple') +
          annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
          xlim(c(0, max_j1j2)) +
          ylim(c(0, max_j1j2)) +
          labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, "1-3kb only")) +
          guides(fill=FALSE) +
          mytheme
        print(p.tmp.1to3k)
        ggsave(paste("Rplot.",n1,"vs",n2,".1to3k.png",sep=''), p.tmp.1to3k, width=8, height=6)
        
        data.class.tmp <- subset(data.class,length_cat=="3-5kb")
        pearson <- round(cor(data.class.tmp[,j1], data.class.tmp[,j2], method="pearson"), 2)
        p.tmp.3to5k <- ggplot(data.class.tmp, aes_string(j1_log10, j2_log10)) +
          geom_point(alpha=0.3, color='royalblue4') +
          annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
          xlim(c(0, max_j1j2)) +
          ylim(c(0, max_j1j2)) +
          labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, "3-5kb only")) +
          guides(fill=FALSE) +
          mytheme
        print(p.tmp.3to5k)
        ggsave(paste("Rplot.",n1,"vs",n2,".3to5k.png",sep=''), p.tmp.3to5k, width=8, height=6)
        
        data.class.tmp <- subset(data.class,length_cat=="5-10kb")
        pearson <- round(cor(data.class.tmp[,j1], data.class.tmp[,j2], method="pearson"), 2)
        p.tmp.5to10k <- ggplot(data.class.tmp, aes_string(j1_log10, j2_log10)) +
          geom_point(alpha=0.3, color='hotpink4') +
          annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
          xlim(c(0, max_j1j2)) +
          ylim(c(0, max_j1j2)) +
          labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, "5-10 kb only")) +
          guides(fill=FALSE) +
          mytheme
        print(p.tmp.5to10k)
        ggsave(paste("Rplot.",n1,"vs",n2,".5to10k.png",sep=''), p.tmp.5to10k, width=8, height=6)
        
        data.class.tmp <- subset(data.class,length_cat==">10kb")
        pearson <- round(cor(data.class.tmp[,j1], data.class.tmp[,j2], method="pearson"), 2)
        p.tmp.ge10k <- ggplot(data.class.tmp, aes_string(j1_log10, j2_log10)) +
          geom_point(alpha=0.3, color='darkolivegreen4') +
          annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
          xlim(c(0, max_j1j2)) +
          ylim(c(0, max_j1j2)) +
          labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, ">10 kb only")) +
          guides(fill=FALSE) +
          mytheme
        print(p.tmp.ge10k)
        ggsave(paste("Rplot.",n1,"vs",n2,".ge10k.png",sep=''), p.tmp.ge10k, width=8, height=6)
        
        
        #grid.arrange(p.tmp.le1k, p.tmp.1to3k, p.tmp.3to5k, p.tmp.5to10k, p.tmp.ge10k, ncol=2)
        
      }
    }
  }
  
  #
  if (nrow(data.FSM) > 0 ) {
    print(p2)
    print(p3)
  }
  if (!all(is.na(data.class$gene_exp))){
    if (nrow(data.class[data.class$structural_category=="NNC",])!=0){
      print(p12)
    }
  }
  #if (!all(is.na(data.class$gene_exp))){
  #    if (nrow(data.class[data.class$structural_category=="NNC",])!=0 & nrow(data.class[data.class$structural_category=="FSM",])!=0 ){]
  #        print(p13)
  #        print(p13.c)
  #    }
  #}
  
  
  #3. splice junction
  
  s <- textGrob("Splice Junction Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p23.a)
  print(p23.b)
  #   print(p24)
  #   print(p25)
  #   print(p26)
  #
  #   if (!all(is.na(data.junction$total_coverage)) & !all(is.na(data.class$iso_exp))){
  #     print(pn1.2)
  #   }
  #
  
  if (!all(is.na(data.junction$total_coverage))) {
    print(pn4.a)
    print(pn4.b)
  }
  
  if (sum(data.junction$RTS_junction=='TRUE') > 0) {
    print(p29.a)
    print(p29.b)
  }
  
  
  s <- textGrob("Comparison With Annotated TSS and TTS", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  if (nrow(data.FSM) > 0) {
    print(p21.a)
    print(p21.b)
    print(p22.a)
    print(p22.b)
  }
  
  if (nrow(data.ISM) > 0) {
    print(p21.dist3.ISM.a)
    print(p21.dist3.ISM.b)
    print(p22.dist5.ISM.a)
    print(p22.dist5.ISM.b)
  }
  
  s <- textGrob("Comparison With Annotated TSS and TTS \nby Subcategories", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  if (!all(is.na(data.FSM$polyA_motif))){
    if (length(p21.FSM.list) > 0) {
      for (i in 1:length(p21.FSM.list)) {
        print(p21.FSM.list[i])
        print(p21.FSM.list.a[i])
      }
    }
  }
  
  if (nrow(data.ISM) > 0) {
    if (length(p21.ISM.list) > 0) {
      for (i in 1:length(p21.ISM.list)) {
        print(p21.ISM.list[i])
        print(p21.ISM.list.a[i])
      }
    }
  }
  
  if (!all(is.na(data.FSM$within_cage_peak))) {
    if (length(p22.FSM.list) > 0) {
      for (i in 1:length(p22.FSM.list)) {
        print(p22.FSM.list[i])
        print(p22.FSM.list.a[i])
      }
    }
    if (length(p22.ISM.list) > 0) {
      for (i in 1:length(p22.ISM.list)) {
        print(p22.ISM.list[i])
        print(p22.ISM.list.a[i])
      }
    }
  }

  s <- textGrob("PolyA Distance Analysis", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  
  #PolyA Distance Analysis by categories
  
  if (sum(!is.na(data.class$polyA_dist)) > 10) {
    print(p.polyA_dist)
    
    # PLOT polyA motif ranking, distance from 3' end
    
    table.polyA <- tableGrob(df.polyA, rows = NULL, cols = c("Category","Count","polyA\nDetected","%"))
    title.polyA <- textGrob("Number of polyA Motifs Detected", gp=gpar(fontface="italic", fontsize=15), vjust=-12)
    gt.polyA <- gTree(children=gList(table.polyA, title.polyA))
    
    
    table.polyA_freq <- tableGrob(df.polyA_freq, rows = NULL, cols = c("Motif", "Count", "%"))
    title.polyA_freq <- textGrob("Frequency of PolyA Motifs", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
    gt.polyA_freq <- gTree(children=gList(title.polyA_freq, table.polyA_freq))
    
    grid.arrange(gt.polyA, gt.polyA_freq, ncol=2)
  }
  
  #PolyA Distance Analysis by subcategories
  
  if (sum(!is.na(data.class$polyA_dist)) > 10) {
    print(p.polyA_dist_subcat)
    print(p.polyA_dist_subcat.s2)
    # PLOT polyA motif ranking, distance from 3' end by subcategory
    
    
    df.polyA_subcat <- tableGrob(df.polyA_subcat, rows = NULL, cols = c("Subcategory","Count","polyA\nDetected","%"))
    title.polyA <- textGrob("Number of polyA Motifs Detected", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
    gt.polyA <- gTree(children=gList(df.polyA_subcat, title.polyA))
    
    

    table.polyA_freq <- tableGrob(df.polyA_freq, rows = NULL, cols = c("Motif", "Count", "%"))
    title.polyA_freq <- textGrob("Frequency of PolyA Motifs", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
    gt.polyA_freq <- gTree(children=gList(title.polyA_freq, table.polyA_freq))
    
    grid.arrange(gt.polyA)
    par(mfrow=c(1,1))
    grid.arrange(gt.polyA_freq)
  }
  
  
  
  if (!all(is.na(data.class$dist_to_cage_peak))) {
    s <- textGrob("CAGE Distances Analysis", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
    grid.arrange(s)
    print(cage_hist_FSM)
    print(cage_hist_FSM_perc)
    if (length(cage.FSM.list)>0) {
      for (i in 1:length(cage.FSM.list)){
        print(cage.FSM.list[i])
        print(cage.FSM.list.a[i])
      }
    }
    print(cage_hist_ISM)
    print(cage_hist_ISM_perc)
    if (length(cage.ISM.list)>0) {
      for (i in 1:length(cage.ISM.list)){
        print(cage.ISM.list[i])
        print(cage.ISM.list.a[i])
      }
    }
    #print(cage_hist_ISM3frag)
    #print(cage_hist_ISM3frag_perc)
    #print(cage_hist_ISM5frag)
    #print(cage_hist_ISM5frag_perc)
    print(cage_hist_NIC)
    print(cage_hist_NIC_perc)
    if (length(cage.NIC.list)>0) {
      for (i in 1:length(cage.NIC.list)){
        print(cage.NIC.list[i])
        print(cage.NIC.list.a[i])
      }
    }
    print(cage_hist_NNC)
    print(cage_hist_NNC_perc)
    if (length(cage.NNC.list)>0) {
      for (i in 1:length(cage.NNC.list)){
        print(cage.NNC.list[i])
        print(cage.NNC.list.a[i])
      }
    }
  
    table.cage <- tableGrob(df.cage, rows = NULL, cols = c("Category","Count","CAGE\nDetected","%"))
    title.cage <- textGrob("Number of CAGE Detected", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
    gt.cage <- gTree(children=gList(table.cage, title.cage))
    
    
    
    table.cage_subc <- tableGrob(df.cage_subc, rows = NULL, cols = c("Subcategory","Count","CAGE\nDetected","%"))
    title.cage_subc <- textGrob("Number of CAGE Detected", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
    gt.cage_subc <- gTree(children=gList(table.cage_subc, title.cage_subc))
    
    
    grid.arrange(gt.cage)
    par(mfrow=c(1,1))
    grid.arrange(gt.cage_subc)
  }
  
  
  s <- textGrob("Redundancy Analysis", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  if (nrow(data.ISM) > 0 || nrow(data.FSM) > 0) {
    print(new.1.FSM)
    print(new.1.ISM)
    print(new.1.total)
    if (!all(is.na(data.class$dist_to_cage_peak))) {
      print(new.2.FSM)
      print(new.2.ISM)
      print(new.2.total)
    }
    if (!all(is.na(data.class$polyA_motif))) {
      print (new.3.FSM)
      print (new.3.ISM)
      print (new.3.total)
    }
    if (!all(is.na(data.class$polyA_motif)) & !all(is.na(data.class$dist_to_cage_peak))) {
      print(new.4.FSM)
      print(new.4.ISM)
      print(new.4.total)
    }
  }
  
  s <- textGrob("Intra-Priming Quality Check", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p30.s1)
  print(p30.s2)
  print(p30.s3)
  print(p31)
  print(p32)
  
  if (saturation.curves=='True'){
    if (!all(is.na(data.class$FL))) {
      s <- textGrob("Saturation curves", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
      grid.arrange(s)
      print(rar1[[1]])
      print(rar1[[2]])
      print(rar2[[1]])
      print(rar2[[2]])
      print(rar3[[1]])
      print(rar3[[2]])
      print(rar5[[1]])
      print(rar5[[2]])
    }
  }
  
  s <- textGrob("Features of Bad Quality", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p28.RTS)
  print(p28.SJ)
  if (n_t3.SJ>0 & n_t3.RTS>0 & !all(is.na(data.class$min_cov))) {
    print(p28.Cov)
  }
  if (n_t3.SJ>0 & n_t3.RTS>0 &!all(is.na(data.class$predicted_NMD))) {
    print(p28.NMD)
  }
  if (n_t3.SJ>0 & n_t3.RTS>0) {
    print(p28)
  }
  
  s <- textGrob("Features of Good Quality", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  if (!all(is.na(data.class$dist_to_cage_peak))) {
    print(p28.a.Cage)
  }
  print(p28.a.annot)
  
  if (!all(is.na(data.class$polyA_motif))) {
    print(p28.a.polyA)
  }
  print(p28.a.SJ)
  
  if (!all(is.na(data.class$min_cov))) {
    print(p28.a.Cov)
  }
  print(p28.a)
  
  if (!all(is.na(data.ratio$ratio_TSS))){
    print(p28.a.ratio.pdf)
  }
  
  dev.off()
  
  print("SQANTI3 report successfully generated!")
}
