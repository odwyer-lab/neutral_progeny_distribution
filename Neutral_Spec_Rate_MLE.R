packages <- c("hypergeo", "grid", "gridBase", "ggplot2", "plyr", "reshape")
for (i in packages){
  if(!require(i, character.only = T)){
    install.packages(i, dependencies = T)
    require(i, character.only = T)
  }
}
cumulasym<-function(b,d,K, thresh=1){
  cumul<-Re(exp((-1 + 2*K)*log(2)-log(b)+K*log(b*d)+(1-2*K)*log(b + d)+lchoose(1/2,K)+log(hypergeo(1,K-1/2,1+K,4*b*d/(b+d)^2))-((-1 + 2*thresh)*log(2)-log(b)+thresh*log(b*d)+(1-2*thresh)*log(b + d)+lchoose(1/2,thresh)+log(hypergeo(1,thresh-1/2,1+thresh,4*b*d/(b+d)^2)))))
  return(cumul)
}  
nu_MLE<-function(progenybyclass){
  #assume that the progeny distribution is given in the form of a two column matrix, where the first column is a vector of total progeny, and the second is a vector of numbers of species with the number of progeny in the first column
  # this returns a maximum likelihood estimator for the combination of neutral parameters 1-b/d
  names(progenybyclass)<-c("births","species")
 return(1-((sum(progenybyclass[,1]*progenybyclass[,2])-sum(progenybyclass[,2]))/sum(progenybyclass[,1]*progenybyclass[,2])))
}

plot_prog_vs_neutral<-function(progenybyclass,returnneutral=TRUE){
  #assume that the progeny distribution is given in the form of a two column matrix, where the first column is a vector of total progeny, and the second is a vector of numbers of species with the number of progeny in the first column
  # this computes a maximum likelihood estimator for the combination of neutral parameters 1-b/d, and then plots the data alongside a neutral progeny dbn with the MLE estimate
colorprog<-c("black","dark gray")
labelvalues<-c("Observed","Neutral Large T Theory\n (ML estimate)")
#estimate nu via MLE
nu_est<-nu_MLE(progenybyclass)
names(progenybyclass)<-c("births","species")
#remove species with zero births from the data
temp<-progenybyclass[which(progenybyclass[,1]>0),]
#compute cumulative version of the same distribution for plotting
cumul<-rev(cumsum(temp[nrow(temp):1,]$species)/sum(temp[,2]))
data_cumulative_prog<-data.frame(births=rev(temp[nrow(temp):1,]$births),prob=cumul)
#compute corresponding neutral prediction for the progeny dbn, using the MLE of speciation rate
neutral_nzs_cumulative_prog<-data.frame(births=data_cumulative_prog$births,theory=cumulasym(1-nu_est,1,data_cumulative_prog$births,1))
dfs<-list(data_cumulative_prog,neutral_nzs_cumulative_prog)
mergeddata<-join_all(dfs,type="full")
plotdata<-melt(mergeddata,id="births")
plotdata<-plotdata[!is.na(plotdata$value),]
progplot<-ggplot(plotdata, aes(x=births,y=value,group=variable,color=variable)) + theme_bw()+scale_linetype_manual(values=c(1,1),labels=labelvalues) +geom_line(aes(y=value))+scale_y_log10()+scale_x_log10()+xlab("births")+ylab("Probability P(k) of number births >= k")+ scale_colour_manual(values = colorprog,labels=labelvalues)+theme(legend.position="none")
plot.new()            
vps <- baseViewports()
pushViewport(vps$figure)
vp1 <-plotViewport(c(0,0,0,0.5))
print(progplot,vp=vp1)                    
popViewport()
  # return neutral progeny dbn if option selected
  if(returnneutral){
  return(neutral_nzs_cumulative_prog)
   }
}
