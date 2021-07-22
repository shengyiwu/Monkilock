#######################################
##### Simple misc functions         ###
#######################################

stdize <- function(x, ...) {
        (x - mean(x,...)) / sd(x, ...)
}

sem <- function(x) { 
        if(length(not.na(x)) < 2) { NA }
        else { sd( not.na(x) ) / sqrt(length( not.na(x) )) }
}

# Create function for distance formula
dist <- function(x1, y1, x2, y2){ sqrt((x2-x1)^2 + (y2-y1)^2) }


not.na <- function(x) { x[!is.na(x)] }

# aggregate by two functions (ie mean and sem)
aggregate2 <- function(x, by, FUN1=mean, FUN2=sem, names=c("x", "mean", "stderr") ) {
        a <- aggregate(x, by=by, FUN=FUN1)
        a2 <- aggregate(x, by=by, FUN=FUN2)
        a$x2 <- a2$x # the rows should be the same since we aggregated by the same thing
        names(a) <- names
        a
}

errbar2 <- function(x,y,yl,yu,col="black",...) {
	par(fg=col)
	errbar(x,y,yl,yu,...)
	par(fg="black")
}

mybin <- function(x, bins) {
	as.numeric(as.character(cut2(x, seq(min(x, na.rm=T), max(x, na.rm=T), (max(x, na.rm=T)-min(x, na.rm=T) + 0.001)/bins), levels.mean=T)))
	}

#######################################
##### Aggregate binomial responses  ###
#######################################

aggregate.binomial <- function(outcome, by, na.rm=T, names=NULL, outcome.name="x", stats=TRUE) {
	
	# k out of n
	k <- aggregate(outcome, by=by, function(x) { sum(x,na.rm=na.rm) })
	n <- aggregate(outcome, by=by, function(x) { length(x[!is.na(x)]) } )
	
	if( !is.null(names) ) { names(k) <- c(names, "x") }
	
	# because we aggregated by the same thing, we can just map the rows
	k$n <- n$x
	k$p <- NA
	k$lower <- NA
	k$upper <- NA
	k$p.value <- NA

	if(stats) {
		
		for(i in 1:dim(n)[1]) {
			
			if(k$n[i] > 0) {
			t <- binom.test(k$x[i], k$n[i])
			k$p[i] <- t$estimate
			k$lower[i] <- t$conf.int[1]
			k$upper[i] <- t$conf.int[2]		
			k$p.value[i] <- t$p.value
			}
		}
	}
	else { # don't run stats, just compute the prob 
		k$p <- k$x / k$n
	}

	# set what we made "x" to be outcome.name
	names(k)[length(names)+1] <-  outcome.name
	
	k
	
}

###########################
##### Plot gam function ###
###########################

plotgam <- function(x, y, sp=15, xlab="", ylab="", ylim=c(0,1), xlim=NULL, main="", cex=2, cex.lab=2, cex.axis=2, lwd=2, col="#f36d90", lty=2, linecol = "black", title) {
	
		par(mar=c(4.5, 4.5, 2, 2) + .02)
       
       # deal with NAs
       keep <- (!is.na(x)) & (!is.na(y))
	x <- x[ keep ]
	y <- y[ keep ]
	
       g <- gam( y ~ 1 + s(x), data=data.frame(x=x,y=y), sp=sp) 

	xes <- seq(min(x), max(x), 0.01)
	newdata <- data.frame(x=xes)
	
	pred <- predict(g, newdata=newdata, type="terms", se.fit=T)
	pred$fit <- pred$fit + coef(g)[1] # apparently gam doesn't include the intercept in its predictions? 
	
	# for some reason, subject gets put before s(seq_item)
	plot(newdata$x, pred$fit[,1], type="l", lwd=lwd, col=col, ylim=ylim, xlim=xlim, ylab=ylab, xlab=xlab, cex.lab=cex.lab, cex.axis=cex.axis, main=main, cex=cex, lty=lty)
 	lines(newdata$x, pred$fit[,1] + pred$se.fit[,1], lwd=lwd, lty=lty, col=linecol)
 	lines(newdata$x, pred$fit[,1] - pred$se.fit[,1], lwd=lwd, lty=lty, col=linecol)
 	title(title, cex.main = 2,   font.main= 2, col.main= "black", adj=0)
 	rug( jitter(x) )
}


plotgam.controlled <- function(f, v, d, sp=10, xlab="", ylab="", ylim=c(0,1), xlim=NULL, main="", cex=2, cex.lab=2, cex.axis=2, lwd=2, col="#f36d90", lty=2, title) {
	
	par(mar=c(4.5, 4.5, 2, 2) + .02)
	
	o <- order(d[,v]) # sort the order so that 
	d <- d[o,]
	
	# make the term to add to our formula
	addf <- update(f, as.formula(paste("~.+ s(", v, ")", sep="")))
	
	g <- gam( addf, data=d, sp=sp) 
	
	pred <- predict(g, type="terms", newdata=d, se.fit=T)
	
	py <- pred$fit[,ncol(pred$fit)]
	pye <- pred$se.fit[,ncol(pred$se.fit)]
	
	# for some reason, subject gets put before s(seq_item)
	plot(d[,v], py, type="l", lwd=lwd, col=col, ylim=ylim, xlim=xlim,
ylab=ylab, xlab=xlab, cex.lab=cex.lab, cex.axis=cex.axis, main=main,
cex=cex, lty=lty)
	title(title, cex.main = 2,   font.main= 2, col.main= "black", adj=0)
 	lines(d[,v], py+pye, lwd=lwd, lty=2)
 	lines(d[,v], py-pye, lwd=lwd, lty=2)
 	rug( jitter(d[,v]) )
}