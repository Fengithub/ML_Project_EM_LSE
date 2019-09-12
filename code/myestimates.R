# source('Seojin/myestimates.R')

my.range <- function(x, xlist, ylist) {
	index = sum(x > xlist)
	early.x = max(xlist[index], 0, na.rm = T)
	late.x = max(xlist[c(index + 1)], x, na.rm = T)
	early.y = min(ylist[index], 1, na.rm = T)
	late.y = min(ylist[c(index + 1)], early.y, na.rm = T)
	if (late.x == early.x) {
		if (early.y == 0) return(sort(unique(ylist))[2])
		return(early.y)
		}
	new.y = early.y + (late.y - early.y)*(1-(late.x - x)/(late.x-early.x))
	return(new.y)
}

my.range.hazard <- function(x, xlist, ylist) {
	# if (x>=130) return(0)
	index = sum(x > xlist)
	early.x = max(xlist[index], 0, na.rm = T)
	late.x = max(xlist[c(index + 1)], x, na.rm = T)
	if(index == 0) {
		early.y = ylist[which.min(xlist)]
	} else {
		early.y = ylist[index]
	}
	late.y = max(ylist[c(index + 1)], min(ylist), na.rm = T)
	if (late.x == early.x) return(early.y)
	#while(late.y == 0) late.y = sort(unique(ylist))[2]; ylist = sort(ylist)[-1]
	new.y = early.y + (late.y - early.y)*(1-(late.x - x)/(late.x-early.x))
	return(new.y)
}

myestimates <- function(x, time, estimates, type = "survival") {
	x <- as.numeric(x)
	time <- as.numeric(time)
	
	estimates <- as.numeric(estimates)
	n <- length(x)
	if (sum(x %in% time) == n) return(estimates[match(x, time)]) 
	if (type == "survival") {
		y <- sapply(x, my.range, xlist = time, ylist = estimates)
		} else if (type == "hazard") {
			y <- sapply(x, my.range.hazard, xlist = time, ylist = estimates) 
		} else {
			return("error")
		}
	return(as.numeric(y))
}
