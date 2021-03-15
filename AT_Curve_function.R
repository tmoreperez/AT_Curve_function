#Thanks for checking out this code for fitting temperature-assimilation curves 
#Be sure to the read the Read.me.txt file.
#Below, the "ATcurve" function can be used to calculate  Topt, Popt, "Omega", and Tmax. These parameters are calculated using 
#two separate models. Model "1" is used to estimate Topt, Popt, "Omega", and model "2" is used to estimate a second set of estimates of Topt, Popt, plus Tmax.
#The model you want to use can be indicated as "1", "2", or both ("3").
#For this function, define the the columns with the 1) leaf temperature & 2) Assimilation variables, then the 3) id/unique factor that you 
#want to estimate the heat tolerances for 4)  indicate if you want to make plots, 6) the number of bootsrap iterations for estimating heat tolerances,
# and 7) the model(s) you want to use.

setwd("/Users/timothyperez/Google Drive/Git_stuff/Git_projects/AT_curve_function/")
atdat=read.csv("Sample_AT.data.csv")
set.seed(23)
#Tleaf=atdat$Tleaf; A=atdat$A; id=atdat$species; plot.est=T; boots; mods=3

ATcurve=function( Tleaf, Assimilation, id, plot.est, boots, mods){
  l1=list(plot.est=plot.est, boots=boots, mods=mods)
  attach(l1)
  ATdf=data.frame(Tleaf=Tleaf, A=Assimilation, id=id)
  return(do.call("rbind", by(ATdf, list(ATdf$id),  function(df){
    Tleaf=df[,which(colnames(df)=="Tleaf")]
    A=df[,which(colnames(df)=="A")]
    id=df[,which(colnames(df)=="id")]
    if(mods==1|mods==3){
      P.opt=T.opt=T.omeg=c()#empty vectors that will hold estimates for Topt and Popt for each iteration of bootstrapping/resampling
      pred.assim=matrix(NA,75,boots)
      for(k in 1:boots){
        srows <- sample(1:length(Tleaf),length(Tleaf),TRUE)#randomly selcts rows from our subset data=
        #If else to catch errors when nls model is parameters are unable to converge
        if(class(try(nls(A[srows]~ (Pop*exp(-((Tleaf[srows] - Top)/(omeg))^2)), start = list(Pop=15,  Top=28, omeg=18)), silent=T)[[1]])=="nlsModel")
        {fit.Tleaf <- nls(A[srows]~ (Pop*exp(-((Tleaf[srows] - Top)/(omeg))^2)),  start = list(Pop=15,  Top=28, omeg=18))
        
        tleafx=seq(from=min(na.omit(Tleaf[srows])), to=max(na.omit(Tleaf[srows])), length.out = 75)
        
        predA=(coef(fit.Tleaf)[[1]]*exp(-((tleafx - coef(fit.Tleaf)[[2]])/(coef(fit.Tleaf)[[3]]))^2))
        df3=data.frame(Tleafs=tleafx, A=predA)
        pred.assim[,k]=predA #These are he predictions from the leaf temperatre model...in other words these are the assimilation rates
        T.opt[k]=coef(fit.Tleaf)[2] #Here we extract the leaf temperature where assiilation is predictd to be highest
        T.omeg[k]=coef(fit.Tleaf)[3]#Here we extract the leaf temperature where assiilation is predictd to be highest
        P.opt[k]=max(df3$A)
        }else{(class(try(nls(A[srows]~ (Pop*exp(-((Tleaf[srows] - Top)/(omeg))^2)),  start = list(Pop=15,  Top=28, omeg=18)), silent=T)[[1]])=="list")
          pred.assim[,k]=NA
          T.opt[k]=NA
          P.opt[k]=NA
          T.omeg[k]=NA}
      }
      #Calculates the 95% confidence interval and the mean for our bootstrapped data
      #A.boot=t(apply(pred.assim, 1, function(x){quantile(x,c(0.025,0.975),na.rm=T)}))
      A.mean=t(apply(pred.assim, 1, function(x){mean(na.omit(x))}))}
    
    if(mods==2|mods==3){
      P.opt2=T.opt2=T.max=T.min=c()#empty vectors that will hold estimates for Topt and Popt for each iteration of bootstrapping/resampling
      pred.assim2=matrix(NA,75,boots)
      #v bootstrapping starts
      for(k in 1:boots){
        srows <- sample(1:length(Tleaf),length(Tleaf),TRUE)#randomly selcts rows from our subset data
        if(class(try(nls(A[srows]~ (b*(Tleaf[srows] - Tmin)*( 1-( exp(cc*(Tleaf[srows]-Tmax)))))^2, start = list(b=.1, cc=.1, Tmax=48, Tmin=13)), silent=T)[[1]])=="nlsModel")
        {fit.Tleaf2 <- nls(A[srows]~ (b*(Tleaf[srows] - Tmin)*( 1-( exp(cc*(Tleaf[srows]-Tmax)))))^2, start = list(b=.1, cc=.1, Tmax=48, Tmin=13))
        tleafx=seq(from=min(na.omit(Tleaf[srows])), to=max(na.omit(Tleaf[srows])), length.out = 75)
        
        pred.assim2.1=(coef(fit.Tleaf2)[[1]]*(seq(from=min(na.omit(Tleaf[srows])), to=max(na.omit(Tleaf[srows])), length.out = 75) - coef(fit.Tleaf2)[[4]])*( 1-( exp(coef(fit.Tleaf2)[[2]]*(seq(from=min(na.omit(Tleaf[srows])), to=max(na.omit(Tleaf[srows])), length.out = 75)-coef(fit.Tleaf2)[[3]])))))^2
        df32=data.frame(Tleafs=tleafx, A=pred.assim2.1)
        pred.assim2[,k]=pred.assim2.1 #These are he predictions from the leaf temperatre model...in other words these are the assimilation rates
        T.opt2[k]=df32[which(df32$A==max(df32$A)),]$Tleafs#Here we extract the leaf temperature where assiilation is predictd to be highest
        T.max[k]=coef(fit.Tleaf2)[3]#Here we extract the leaf temperature where assiilation is predictd to be highest
        P.opt2[k]=max(df32$A)
        T.min[k]=coef(fit.Tleaf2)[4]
        #aftop=df3[which(df3$Tleafs>30),]
        #T.max[k]=which(abs(aftop$A-1)==min(abs(aftop$A-1)))
        #lines(pred.df$GDs, pred.GD, col=alpha(df2$kolor, 0.025))
        }else{(class(try(nls(A~ (b*(Tleaf - Tmin)*( 1-( exp(cc*(Tleaf-Tmax)))))^2, data=df2, start = list(b=.1, cc=.1, Tmax=48, Tmin=13)), silent=T)[[1]])=="list")
          pred.assim2[,k]=NA
          T.opt2[k]=NA
          P.opt2[k]=NA
          T.max[k]=NA
        }}
      #Calculates the 95% confidence interval and the mean for our bootstrapped data
      #A2.boot=t(apply(pred.assim2, 1, function(x){quantile(x,c(0.025,0.975),na.rm=T)}))
      A2.mean=t(apply(pred.assim2, 1, function(x){mean(na.omit(x))}))}
    
    if(plot.est==T){
      plot(Tleaf,A, col="black", ylim=c(0,30), pch=5, xlim=c(25,46), 
           ylab=expression("Photosynthesis ("~mu~"mol m"^-2~"s"^-1~")"),
           xlab=expression("Leaf temperature ("~degree~ "C)") )
      text(25,28,pos=4, paste(unique(paste(id))), font=4,cex=.7)
      
      if(mods==1|mods==3){ 
        lines(seq(from=min(na.omit(Tleaf)), to=max(na.omit(Tleaf)), length.out = 75), A.mean, col="black", lwd=1.5)
      }
      if(mods==2|mods==3){
        lines(seq(from=min(na.omit(Tleaf)), to=max(na.omit(Tleaf)), length.out = 75), A2.mean, col="gray", lwd=1.5)
      }}
    
    mean.Topt=NA
    mean.Popt=NA
    mean.Tomeg=NA
    if(mods==1|mods==3){
      mean.Topt=mean(na.omit(T.opt)) 
      mean.Popt=mean(na.omit(P.opt))
      mean.Tomeg=mean(na.omit(T.omeg))}
    
    mean.Topt2=NA
    mean.Popt2=NA
    mean.Tmax=NA
    if(mods==2|mods==3){
      mean.Topt2=mean(na.omit(T.opt2)) 
      mean.Popt2=mean(na.omit(P.opt2))
      mean.Tmax=mean(na.omit(T.max))}
    
    return(data.frame(id=paste(unique(id)), mean.Topt=mean.Topt, mean.Popt=mean.Popt, mean.Tomeg=mean.Tomeg,
                      mean.Topt2=mean.Topt2, mean.Popt2=mean.Popt2, mean.Tmax=mean.Tmax))
  }) ))
  detach(l1)}

#When plotting model 2, I've been noticing that at high (and sometimes low) temperatures, the model tends to increase
#I think this is an artifact of the bootstrapping procedure when the model is re-fit. This increase in modeled A at high temperature
#occurs when temperatures higher than those used to fit the model are used to predict new values of A. 
ATcurve(Tleaf=atdat$Tleaf, A=atdat$A, id=atdat$species, plot.est=T, boots=100, mods=2)