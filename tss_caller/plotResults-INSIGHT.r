library("Hmisc",lib.loc="/usr/projects/INSIGHT/webtool/Rlibs/");

# I/O FILES
args          <- commandArgs(trailingOnly=TRUE);
tableFileName <- args[1];
plotName      <- args[2];
pValArgs      <- args[3];

# GENERAL SETTINGS
params=expression(paste("E[",P[w],"] / kbp"),rho,paste("E[",D[p],"] / kbp"));
paramCols=c(5,1,3);
lrtStatCols=c(18,16,17)
pValCutoffs=rep(pValArgs,3);     # CUTOFFS FOR P-VALUE
pValDF=c(1,3,1);                   # DEGREES OF FREEDOM FOR P-VAL CALCULATIONS
pValLabels=c("w","*","p");
barColors=c("grey95","grey80","grey50");
barWidths=c(0.8,1,0.8);
heightScaling=c(0.25,1,0.25);
barSpace=c(0,1);
plotHeight=6;      # INCHES
plotWidthProp=0.6; # per set proportion between width and height

raw_table<- read.table(tableFileName);

if(raw_table[1,1] == "dataID" && raw_table[1,3] == "rho") {
   raw_table=raw_table[2:nrow(raw_table),];
}
numSets<-nrow(raw_table);

if(numSets == 1) {
   numSets = 2;
   raw_table_tmp <- raw_table;
   raw_table<-matrix(0.0,2,ncol(raw_table_tmp));
   raw_table[1,]<-as.matrix(raw_table_tmp[1,]);
   raw_table[2,1:2]<-"";
   raw_table[2,2+paramCols]<--0.1;
   #plotWidthProp=plotWidthProp*1.5;
}

primaryLabels   <- raw_table[,1];
secondaryLabels <- raw_table[,2];

data_table      <- matrix(0.0, numSets,ncol(raw_table)-3);
data_table[,]   <- as.numeric(as.matrix(raw_table[,4:ncol(raw_table)-1]));

pdf(plotName,width=2+plotHeight*plotWidthProp*numSets,height=plotHeight);
par(mar=c(6,6.2,1,6.7));

# PLOT RESULTS
height=t(data_table[,paramCols])*heightScaling;
stds=t(data_table[,paramCols+1])*heightScaling;
noStds=which(stds<0)
stds[noStds] = height[noStds];
errBarColors=matrix("grey20",1,length(stds));
errBarColors[noStds] = "red";
ylim=c(0,max(height+stds,height[1:6]+stds[1:6]+0.4));

barPositions = 
   barplot ( beside=TRUE , height=height , width=barWidths , space=barSpace ,
             col=barColors , ylim=ylim , axes=FALSE, xpd=FALSE);

   errbar ( add=TRUE , x=barPositions , y=height , yplus=height+stds , # yplus=pmin(height+stds,1.0) ,
            yminus=pmax(height-stds,0.0) , lwd=2 , cap=0.1*barWidths/numSets ,
            pch='' , errbar.col=errBarColors, xpd=FALSE,);

# CLEAR LINE AT BASE OF EMPTY SECOND RESULT IF THERE IS ONLY ONE RESULT TO PLOT
if(numSets == 2 && primaryLabels[2] == "") {
   segments(x0=barPositions[1,2]-0.5, x1=barPositions[3,2]+0.5,y0=-0.0,y1=-0.0,col="white",lwd=20);
}


# LABEL EACH RESULT
mtext(side=1 , at = colMeans(barPositions) , line=2.5, xpd=TRUE, text=primaryLabels , cex=2);
mtext(side=1 , at = colMeans(barPositions) , line=4, xpd=TRUE, text=secondaryLabels , col="grey30", cex=1.5);

# LEFT AXIS FOR rho
axis(side=2,cex.axis=2.0,las=2);
mtext(side=2,line=4.5,text=params[2], cex=2.5);

# RIGHT AXIS FOR E[Dp] and E[Pw]
rightAxisLim=ylim[2]/heightScaling[1];
labels=seq(0,rightAxisLim,0.5);
axis(side=4,at=labels*heightScaling[1],labels=labels,cex.axis=2.0,las=2);
mtext(side=4,line=5.5,text=params[c(1,3)], at=ylim[2]*c(0.2,0.75), cex=2.5);

# SIGNIFICANCE LABELS
sigLabels=matrix("",3,numSets);
for(i in 1:3) {
   pvals=1-pchisq(df=pValDF[i],q=data_table[,lrtStatCols[i]]);
   sigLabels[i,which(pvals<pValCutoffs[i])] = pValLabels[i];
   print(pvals)
}
text(x = barPositions[2,], -0.03, xpd=TRUE, labels=sigLabels[2,],cex=4,col="darkred");
text(x = barPositions[c(1,3),], -0.03, xpd=TRUE, labels=sigLabels[c(1,3),],cex=2,col="darkred");


legend(legend=params,fill=barColors,x='topleft',inset=-0.00,bty='n',horiz=FALSE,cex=2, y.intersp=1);
sigLabels = c(
   substitute(paste(' ',gamma,' > 0    (p < ',x,')'),list(x=pValCutoffs[1])),
   substitute(paste(' ',rho,  ' > 0    (p < ',x,')'),list(x=pValCutoffs[2])),
   substitute(paste(' ',eta,  ' > 0    (p < ',x,')'),list(x=pValCutoffs[3])),
   expression(rho));
sigLabel2 = substitute(paste(' (',gamma,' > 0 ; p < ',x,')'),list(x=pValCutoffs[1]));
sigLabel3 = substitute(paste(' (',gamma,' > 0 ; p < ',x,')'),list(x=pValCutoffs[1]));
legend(  legend=sigLabels[1:3],pch=pValLabels,col="darkred", y.intersp=1.2,
         x=3.5-numSets*0.1,y=ylim[2],bty='n',horiz=FALSE,pt.cex=c(2,3,2),cex=1.8);



dev.off();

