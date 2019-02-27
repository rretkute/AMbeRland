# setwd("/Users/renataretkute/Documents/GitHub/AMbeRland")
#library(factoextra)
#library(FactoMineR)


###################################################
### Shiny app for visualisation of MIC densities
###################################################


require(shiny)
require("RColorBrewer")
library(Rtsne)


#  Data upload and analysis
data<-read.table("data.txt", header=T)

antibiotics<-as.character(sort(unique(data$Antibiotic)))
countries<- as.character(sort(unique(data$Country)))
years<-seq(2004, 2019, by=1)
val.log2.MIC<-seq(-8,8)

ylim<-c(min(val.log2.MIC),max(val.log2.MIC))
xlim<-c(min(years),max(years))

ids<-unique(data$ID)
isol<-matrix(0, nrow=length(ids), ncol=length(antibiotics))
ii<-0
for(i in 1:length(ids)){
	id<-ids[i]
	MICs<-as.numeric(sapply(1:length(antibiotics), function(a) data$MIC[data$ID==id & data$Antibiotic == antibiotics[a]]))
	if(length(which(is.na(MICs)))==0){
		ii<-ii+1
		isol[ii,] <- log(MICs, base=2)
	}
}
isol<-data.frame(isol[1:ii,])
colnames(isol)<-antibiotics
#head(isol)

# PCA analysis
res.pca <- prcomp(isol)
prj.pca<-data.frame(res.pca$x[,c(1,2)])

# t-SNE analysis
tsne<-Rtsne(isol, check_duplicates =F, pca=F, perplexity=50)


# Functions

H <- function(x) as.numeric(x>0)

legend.col <- function(col, lev){  
	opar <- par  
	n <- length(col) 
	 bx <- par("usr")  
	 box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000, bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50) 
	 box.cy <- c(bx[3], bx[3]) 
	 box.sy <- (bx[4] - bx[3]) / n  
	 xx <- rep(box.cx, each = 2)  
	 par(xpd = TRUE) 
	 for(i in 1:n){  
	 	yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i - 1))) 
	 	polygon(xx, yy, col = col[i], border = col[i])  
	 	} 
	 par(new = TRUE) 
	 plot(0, 0, type = "n", ylim = c(min(lev), max(lev)), yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE) 
	 axis(side = 4, las = 2, tick = FALSE, line = .25) 
	 par <- opar 
	 }




##########################
#  shiny interface
##########################

ui = fluidPage(

	#The title
	titlePanel("Analysis and visualisation of antimicrobial resistance"),
	
	
	#The sidebar for parameter input
	sidebarPanel(
	
	radioButtons("analysis", "Type of plot:", c("Annual distribution" = "distr", "PCA visualisation of log2(MIC)" = "viz1", "t-SNE visualisation of log2(MIC)" = "viz2")),
	
	      br(),

  	 selectInput("cntr", "Choose country (for distribution, Worldwide for visualisation):", choices = c("Worldwide", countries)),
    	 selectInput("antb", "Choose antibiotic:", choices = antibiotics)
	),
	#Main panel for figures and equations
	mainPanel(        
	# Output     
		plotOutput(outputId = "plot"))
) #End of ui()


server = function(input, output) {
 
  #Plot: renderPlot to be passed to UI 
  output$plot = renderPlot({
  	if(input$analysis=="distr"){ 
  		# Calculate distribution
  		fr.sampl<-matrix(0, ncol=length(val.log2.MIC), nrow=length(years))
	for(i in 1:length(years)){
		if(input$cntr=="Worldwide") {
				wh<-as.character(data$MIC[which(data$Antibiotic==input$antb & data$Year==years[i])])
			} else {
				wh<-na.omit(as.character(data$MIC[which(data$Antibiotic==input$antb &  data$Country==input$cntr & data$Year==years[i])]))}
	wh<-round(log(suppressWarnings(as.numeric(wh)),base=2))
	fr.sampl[i,]<-sapply(1:length(val.log2.MIC), function(a) length(which(wh==val.log2.MIC[a])))
	if(sum(fr.sampl[i,])>0) fr.sampl[i,]<-fr.sampl[i,]/sum(fr.sampl[i,])
	}
	toti=c(sapply(1:nrow(fr.sampl), function(a) fr.sampl[a,]))
	toty=c(sapply(1:length(years), function(a) rep(years[a],length(val.log2.MIC)))) 
	totm=rep(val.log2.MIC,length(years))  
	tot=cbind(toty,totm,toti)  

	# Plot 
	plot(toty, totm, cex=H(toti)*3, pch=16, col='gray', xlab="Year", ylab="log2(MIC)", xlim=xlim, ylim=ylim, main='', asp=1) 
	points(toty, totm, cex=toti*3, pch=16, col='black')  
  	} 
  	if(input$analysis=="viz1"){ 
 		wh.ant<-which(antibiotics==input$antb)
 		if(length(wh.ant)>0){
 	colours <- factor(isol[,wh.ant], levels=val.log2.MIC)
	colours <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(unique(colours)))[factor(colours)]
 			
 			# Plot
 				plot(prj.pca$PC1, prj.pca$PC2,col=colours, pch=1, asp=1, xlab="PC1", ylab="PC2", main="")
	legend.col(col = rev(brewer.pal(9, "RdBu")), lev = seq(-8,8)) 
 	  }	
 	 }
 	if(input$analysis=="viz2"){ 
 		wh.ant<-which(antibiotics==input$antb)
 		if(length(wh.ant)>0){
 	colours <- factor(isol[,wh.ant], levels=val.log2.MIC)
	colours <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(unique(colours)))[factor(colours)]
 			
 			# Plot  
 				plot(tsne$Y, col=colours, pch=1, asp=1, xlab="Dim1", ylab="Dim2", main="")
legend.col(col = rev(brewer.pal(9, "RdBu")), lev = seq(-8,8))
 	  }	
 	 }
   })
 } #End of server()

shinyApp(ui, server)
