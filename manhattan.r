library(scales)

manhattan<-function(x,cutoff=1,type="",Pcol="p",CHRcol="chr",POScol="pos",RScol="rs",tophit=F,MHC=F,realpos=F)
{
	chrlength<-c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
	chrsum<-c(0,cumsum(chrlength))
	if (type=="plink")
	{
		Pcol<-"P"
		CHRcol<-"CHR"
		RScol<-"SNP"
		POScol<-FALSE
		realpos<-FALSE
	}
	if (type=="gemma")
	{
		Pcol<-"p_wald"
		CHRcol<-"chr"
		RScol<-"rs"
		POScol<-"ps"
	}
	dat<-x[x[,Pcol]<cutoff,]
	print(dim(dat))
	chr.uniq<-unique(dat[,CHRcol])
	cols<-rep(c("black","blue"),11)[match(dat[,CHRcol],chr.uniq)]
	if (realpos==F)
	{
		xpos<-1:nrow(dat)
		chrmid<-tapply(1:nrow(dat),dat[,CHRcol],mean)
		
	} else 
	{
		xpos<-apply(dat,1,function(d){
                        chrsum[as.numeric(as.character(d[CHRcol]))]+as.numeric(as.character(d[POScol]))
                })
		chrmid<-sapply(2:23,function(i){(chrsum[i-1]+chrsum[i])/2})
		
	}
	plot(xpos,-log10(dat[,Pcol]),col=cols,pch=20,xaxt="n",xlab="Chromosome",ylab="-log10(pval)")
	axis(side=1,at=chrmid,labels=1:22)
	if (tophit==T)
	{	
		mtext(at=xpos[order(dat[,Pcol])[1]],side=3,line=1,text=dat[order(dat[,Pcol])[1],RScol])
	}
	if (MHC==T)
	{
		if (POScol!=FALSE)
		{
			idx<-(which(dat[,CHRcol]==6 & dat[,POScol]>28477797 & dat[,POScol]<33448354))
			points(xpos[idx],-log10(dat[idx,Pcol]) ,col="red")

		}
	}	
}


qq<-function(x,type=NULL,Pcol="P",CHRcol="chr",RScol="rs"){
	if (type=="plink")
	{
		Pcol<-"P"
		CHRcol<-"CHR"
		RScol<-"SNP"
		POScol<-FALSE
	}
	if (type=="gemma")
	{
		Pcol<-"p_wald"
		CHRcol<-"chr"
		RScol<-"rs"
		POScol<-"ps"
	}

	
	observed <- sort(x[,Pcol])
	lobs <- -(log10(observed))

	expected <- c(1:length(observed)) 
	lexp <- -(log10(expected / (length(expected)+1)))



	plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
	points(lexp, lobs, pch=23, cex=.4, bg="black",type="l") 

}
