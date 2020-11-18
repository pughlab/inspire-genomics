#############
# "S-plot" function for multiple groups
# Updated: Nov 18, 2020
# Created by : Cindy Yang
# Generate ordered, grouped, scatter plots 
# Updated to include log-scale option for y-axis
# Example: 
# source (~/R/testScores.R)
# main.df <- read.csv ("~/data/Fig1.csv", stringsAsFactors = FALSE)
# testScores (dat = main.df, xcat = "COHORT", ycat = "TMB_n72", log_y = TRUE, 
#             cex.lbs = 0.8, grid = T, pt.col = COLS_CANCER_COHORT)
############

# Helper function to calculate log10 axis tick positions
#https://stackoverflow.com/questions/47890742/logarithmic-scale-plot-in-r
log10Tck <- function(side, type){
  lim <- switch(side, 
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceiling(lim[2])
  return(switch(type, 
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}

##Splot function
testScores <- function (dat, xcat, ycat, orderByMedian = FALSE,
                        cex.pts = 1, cex.lbs = 1, box = FALSE, pt.col = "black", 
                        log_y = FALSE, grid = FALSE){
  # plot boxplot between categories
  par(mar=c(5,5,2,3))
  # remove any entry with missing data
  dat <- dat[!is.na(dat[,ycat]),]
  box.stats <- boxplot (formula (paste (ycat, xcat, sep = "~")), 
                        data = dat, xlab = xcat, ylab = ycat, plot = FALSE)
  group.order <- box.stats$names
  if(orderByMedian) group.order <- box.stats$names[order(box.stats$stats[3,], decreasing = FALSE)]
  dat[,xcat] <- factor(dat[,xcat], levels = group.order, order = TRUE)
  dat <- dat[order(dat[,xcat], dat[,ycat], decreasing = FALSE),]
  print ("ordered data")
  print(dim(dat))
  print(box.stats$n)
  # Calculate the order of the sigmoid plots
  calculateX <- function (n, mp = seq(1,length(n)), width = 0.4){
    min_x <- mp -width/2
    max_x <- mp +width/2
    return(unlist(lapply(mp, function (x) return(seq (min_x[x], max_x[x], length.out = n[x])))))
  }
  print(length(calculateX (box.stats$n)))
  dat$x.coord <- calculateX (box.stats$n)
  if(orderByMedian) dat$x.coord <- calculateX (box.stats$n[order(box.stats$stats[3,], decreasing = FALSE)])
  medians <- box.stats$stats[3,]
  print(medians)
  if (orderByMedian) medians <- sort (box.stats$stats[3,], decreasing = FALSE)
  # To recolor a specific group, need to add another column of colors to be used for colors
  # TODO:
  # dat$x.cols <- factor(as.character(dat[,xcat]), levels = names(group.cols), labels = group.cols)
  # Boxplot and sigmoid scatter
  plot(x = dat$x.coord, 
       y = dat[,ycat],  
       log = ifelse (log_y, yes = "y", no = ""),
       xlim = c(0,length(group.order)+1), 
       xlab = xcat, ylab = ycat, xaxt = "n",
       yaxt =  ifelse (log_y, yes = "n", no = "s"))
  
  if (log_y){
    # plot y-axis ticks (major and minor)
    axis(2, at=log10Tck("y",'major'), tcl= -0.5, las = 2, cex.axis = 0.80) # left
    axis(2, at=log10Tck("y",'minor'), tcl= -0.3, labels=NA) # left
  }

  # Add boxplot outline if required
  if (box){
    tryCatch(expr = 
            {boxplot (formula (paste (ycat, xcat, sep = "~")), data = dat, 
             log = ifelse (log_y, yes = "y", no = ""),
             xlab = xcat, ylab = ycat, ylim = ylimit,
             add = TRUE, boxwex = 0.5, las = 2, border = "grey",
             boxlwd = 2, outline = FALSE,
             cex.axis = cex.lbs)},
            error = function(e) {
              message("Cannnot plot boxplot, check if data has zeros")
              return (NA)
            })
  }
  # Plot points
  points(x = dat$x.coord, 
         y = dat[,ycat], 
         ylog = log_y,
         pch = 16, col = pt.col, cex = cex.pts)
  
  # Plot medians
  segments (x0 = seq(1,length(box.stats$n))-0.2,
            x1 = seq(1,length(box.stats$n))+0.2,
            y0 = medians, col = "red", lwd = 2)
  # Plot x axes labels
  n.labs <- box.stats$n
  if (orderByMedian) n.labs <- box.stats$n[order(box.stats$stats[3,])]
  axis (side = 1, at = seq(1, length(box.stats$n), by = 1), labels = group.order, las = 2, cex.axis = cex.lbs)
  axis (side = 3, at = seq (1, length(box.stats$n), by = 1), labels = n.labs, las = 1, cex.axis = cex.lbs)
  if (grid){
    grid (NULL,NULL, lty = 6, col = "cornsilk2")
  }
  # do wilcox test or kruskal test
  if (length(box.stats$n) ==2){
    stats <- wilcox.test (formula (paste (ycat, xcat, sep = "~")), data = dat)
  }else{
    stats <- kruskal.test (formula (paste (ycat, xcat, sep = "~")), data = dat)
  }
  p.val <- ifelse (test = round (stats$p.value, digits = 3) < 0.001,
                   yes = "p < 0.001",
                   no = paste0 ("p=",round (stats$p.value, digits = 3)))
  legend ("topright", legend = p.val, bty = "n", cex = cex.lbs)
}