corrWheel <- function( Data, Threshold, Groups=NULL, plotCex=1, plotMain='' ) {
# corrWheel Graphs a "correlation wheel," which is a graphic with names of each factor arranged in a circular fashion with lines drawn between each pair of factors that have an absolute correlation value >= some threshold.  A correlation wheel is useful for choosing variables for analysis that have acceptible levels of correlation between them (i.e., are not connected by a line on the wheel).
#
# ARGUMENTS
# Data [required]
# data frame with data to be analyzed... pairwise Pearson correlations are calculated for each possible pair of columns
#
# Group [optional]
# data frame with first column the name of a field in Data and the second a numeric or character value indicating which group the value should be considered part of... this is used to color-code names of factors in the corr wheel to aid visual grouping of like factors
#
# Threshold
# numeric value bewteen 0 and 1 above which absolute value of a Pearson correlation coefficient is considered "too corrrelated" and so indicate the connection of the respective factors with a line
#
# plotCex
# cex for plot (character size)
#
# plotMain
# character, text for main title of plot
#
# VALUES
# (graph)
#
# BAUHAUS
# 
#
# EXAMPLE
# corrWheel( Data=data.frame(A=runif(100), B=1:100*runif(100), C=1:100, D=sort(runif(100)) ), Groups=data.frame(Name=c('A', 'B', 'C', 'D'), Group=c(1, 1, 2, 2)), Threshold=0.6, plotCex=1, plotMain='TITLE' )
#
# SOURCE	source('C:/Work/Code/Utility - Correlation Wheel.r')
#
# LICENSE
# This document is copyright ©2012 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/
#
# AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
# DATE		2012-01
# REVISIONS 

#############################
## libraries and functions ##
#############################

if (!is.null(Groups)) require(colorspace)

##########################
## initialize variables ##
##########################

lineWidthMult <- 2 # 1.4 # multiply line width by this factor first
lineWidthPow <- 3 # 2.5 # then take line width to this power

lineColors <- rev( heat.colors(n=10, alpha=0.4) )
lineColorsLegend <- rev( heat.colors(n=10) )

##########
## MAIN ##
##########

circPos <- seq( 0, 2 * pi - 2 * pi / ncol(Data), length.out=ncol(Data) ) # get circular coordinates for each factor

# set up plot
par( cex=plotCex, bty='n', pty='s', cex.main=1, bg='black', col='white')
plot(x=c(-1, 1), y=c(-1,1), col='black', col.main='white', xaxt='n', yaxt='n', xlab=NA, ylab=NA, main=paste(plotMain, ' thresholded at ', Threshold, sep='')) # blank plot
	
for (countFact1 in 1:(length(names(Data)) - 1) ) { # for each pair of factors, link with line

	for (countFact2 in 2:length(names(Data))) { # for second factor
	
		thisCor <- abs( cor(x=Data[ , names(Data)[countFact1]], y=Data[ , names(Data)[countFact2]], method='pearson', use='pairwise.complete.obs') )
	
		if ( thisCor >= Threshold ) { # draw line if too correlated

			lines(
				x=c( sin(circPos[countFact1]), sin(circPos[countFact2]) ),
				y=c( cos(circPos[countFact1]), cos(circPos[countFact2]) ),
				lty='solid',
				lwd=( lineWidthMult * thisCor )^lineWidthPow,
				col=lineColors[ ceiling(10* thisCor^2) ]
			)
				
		}
		
	} # for second factor
	
} # for each pair of factors, link with line

## plot names of each factor

if (is.null(Groups)) { # get text colors by group

		labelColor <- 'white' # all black because no groups
	
} else {
		
		# each group a different color
		groupNames <- unique( as.character(Groups[ , 2]) ) # get group names
		groupColors <- rainbow_hcl(n=length(groupNames), start=30, end=300) # number of colors equal to group names
		labelColor <- rep(NA, length(names(Data)))
		for (countGroups in 1:(length(groupNames))) labelColor[as.character(Groups[ , 2])==groupNames[countGroups]] <- groupColors[countGroups]

}

# plot factor names
text(
	x=sin(circPos),
	y=cos(circPos),
	labels=names(Data),
	adj=c(0.5, 0.5),
	col=labelColor,
	xpd=NA
	
)

legend(
	x=0,
	y=-1.1,
	legend=paste( '±', round( seq(from=1, to=Threshold, length.out=4), digits=2), sep='' ),
	lwd=( seq(from=1, to=Threshold, length.out=4) * lineWidthMult ) ^ lineWidthPow,
	col=lineColorsLegend[ round(10 * seq(from=1, to=Threshold, length.out=4)^2) ],
	xpd=NA,
	ncol=4,
	adj=0.5,
	bty='n',
	cex=1,
	xjust=0.5
)

}

