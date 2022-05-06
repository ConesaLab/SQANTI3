### Fran's code for ploting saturation plots
plot.rarefaction <- function (rarefact, k = 0, sample = 1, depth.increase = 2, break.by = "none", 
                              specify = "all transcripts", what = c("sat", "diff")) {
  force.sat <- function(x) {
    for (t in 2:(length(x) - 1)) {
      if (x[(t + 1)] < x[t]) {
        x[(t + 1)] = x[t]
      }
    }
    x
  }
  myPalette = c("#6BAED6","#FC8D59","#78C679","coral2","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")
  if (break.by == "none") {
    all.sample <- length(rarefact$simulation)
    transcripts <- rarefact$simulation[[all.sample]]
    ticks <- rarefact$ticks[[all.sample]]
    gene.names <- names(transcripts[[length(ticks)]][[1]])
    selection <- data.frame(rep("all transcripts", length(gene.names)))
    rownames(selection) <- gene.names
  }
  if (break.by == "biotype") {
    selection = rarefact$biotype
  }
  if (break.by == "category") {
    selection = rarefact$category
  }
  for (o in 1:sample){
    #ticks_df=data.frame(row.names = levels(as.factor(rarefact$category[,1])))
    ticks_df=data.frame(Row.names=levels(as.factor(selection[,1])))
    ticks=rarefact$ticks[[o]]
    for (t in 1:length(rarefact$simulation[[o]])){
      reps_df=data.frame()
      cat_df=data.frame(matrix(ncol = 1, nrow = 0))
      colnames(cat_df) <- c("category")
      for (r in 1:length(rarefact$simulation[[o]][[t]])){
        tmp_r=as.data.frame(rarefact$simulation[[o]][[t]][[r]])
        reps_df=merge(reps_df, tmp_r, all=T, by.x = 0, by.y = "Var1")
        #reps_df=reps_df[,-which(colnames(reps_df)=="Row.names")]
        rownames(reps_df)=reps_df[,1]
        colnames(reps_df)[which(colnames(reps_df)=="Freq")] <- paste("rep",r, sep = "_")
        if (any(colnames(reps_df)=="ID")){reps_df=reps_df[,-which(colnames(reps_df)=="ID")]}
        colnames(reps_df)[which(colnames(reps_df)=="Row.names")] <- "ID"
        #tmp_c=aggregate(cbind(count=as.character(paste("rep",r, sep = "_"))) ~ category , data=reps_df, FUN=function(x) length(which(x>=k)))
      }
      reps_df=merge(reps_df, selection, by=0, all=T)
      colnames(reps_df)[which(colnames(reps_df)=="featureData(input)$Category")] <- "category"
      for (c in grep(x=colnames(reps_df), patter="rep")){
        counts_cat=aggregate(cbind(count=reps_df[,c]) ~ category , data=reps_df, FUN = function(x) length(which(x>=k)) )
        colnames(counts_cat)[which(colnames(counts_cat)=="count")] <- colnames(reps_df)[c]
        cat_df=merge(cat_df, counts_cat, all=T, by = "category")
      }
      rownames(cat_df)=cat_df$category
      #cat_df=cat_df[, - which(colnames(cat_df)=="category")]
      cat_df$avg=apply(cat_df[, - which(colnames(cat_df)=="category")],1,mean)
      ticks_df=merge(ticks_df, cat_df[,c("category","avg")] , by.x="Row.names", by.y = "category", all = T)
      colnames(ticks_df)[which(colnames(ticks_df)=="avg")] = paste("tick",t,sep="_")
     
    }
  } 
  tticks_df=t(ticks_df)
  colnames(tticks_df)=tticks_df["Row.names",]
  tticks_df=tticks_df[-1,]
  #cbind(tticks_df, ticks=ticks)
  tticks_df=mutate_all(data.frame(cbind(tticks_df, ticks=ticks)),function(x){as.numeric(as.character(x))})
  more.ticks <- as.data.frame(round(tail(ticks, n = 1) * seq(0.1, depth.increase, 0.1), 0))
  colnames(more.ticks) = "ticks"
  seq_depth=max(tticks_df$ticks)/10^6
  pred_df=data.frame(ticks=more.ticks)
  for (c in 1:(length(colnames(tticks_df))-1)){
    fm1 = loess(tticks_df[,colnames(tticks_df)[c]] ~ ticks, tticks_df, control = loess.control(surface = "direct"), degree=2)
    prediction = cbind(predict(fm1, newdata = more.ticks), more.ticks)
    colnames(prediction)[1]=colnames(tticks_df)[c]
    pred_df=merge(pred_df,prediction, by="ticks")
  }
  tticks_df=rbind(tticks_df, pred_df[(length(tticks_df[,1])+1):length(pred_df[,1]),])
  tticks_df=as.matrix(tticks_df)
  tticks_df[which(is.na(tticks_df))]<-0
  tticks_df=as.data.frame(tticks_df)
  sat_tticks_df=data.frame(apply(tticks_df,2,function(y){force.sat(y)}))
  sat_tticks_df$ticks=sat_tticks_df$ticks/10^6
  all_ticks=sat_tticks_df$ticks/10^6
  sat_tticks_df[,- which(colnames(sat_tticks_df)=="ticks")]=sat_tticks_df[,- which(colnames(sat_tticks_df)=="ticks")]/1000
  increment=sat_tticks_df[(2:length(sat_tticks_df$ticks)),] - sat_tticks_df[(1:(length(sat_tticks_df$ticks)-1)),]
  increment_corr=cbind(increment[,- which(colnames(increment)=="ticks")]/increment$ticks, ticks=sat_tticks_df$ticks[-length(sat_tticks_df$ticks)])
  
  xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic.Genomic",  "Antisense", "Fusion","Intergenic", "Genic.Intron")
  
  saturation_df=reshape::melt(sat_tticks_df, id.vars="ticks")
  saturation_df$variable=factor(saturation_df$variable,
         levels = xaxislabelsF1, 
         ordered=TRUE)
  
  increase_df=reshape::melt(increment_corr, id.vars="ticks")
  increase_df$variable=factor(increase_df$variable,
                              levels = xaxislabelsF1, 
                              ordered=TRUE)
  sat_plot=ggplot(saturation_df, aes(x=ticks, y=value, fill=variable, color=variable)) +
    geom_line(linetype = "dashed", show.legend = FALSE)+
    scale_color_manual(values=myPalette) +
    geom_point(aes(color=variable), show.legend = TRUE) +
     #theme(legend.position="none") +
    mytheme +
    ylab("Number of isoforms detected (in thousands)") +
    xlab("Sequencing depth \n (millions of FL sequences)") +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=10))+
    labs(title = "Saturation plot per structural category\n\n",
         subtitle = paste("Minimum number of FL-reads required to call an isoform = ",k, "\n\n", sep=" ") ) +
    geom_vline(xintercept = seq_depth)
  inc_plot=ggplot(increase_df, aes(x=ticks, y=value, fill=variable)) +
    geom_bar(stat="identity", position = "dodge")+
    scale_fill_manual(values=myPalette) +
    mytheme +
    ylab("Absolute increment of isoforms detected (in thousands)") +
    xlab("Sequencing depth \n (millions of FL sequences)") +
    theme(legend.position="bottom", legend.title=element_blank()) +
    theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
    labs(title = "Increments of detected isoforms per structural category\n\n",
         subtitle = paste("Minimum number of FL-reads required to call an isoform = ",k, "\n\n", sep=" ") ) +
    geom_vline(xintercept = seq_depth)
  list(sat_plot,inc_plot)
}
