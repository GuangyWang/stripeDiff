set.seed(123)
suppressMessages(require(bcp))
suppressMessages(require(plotly))
suppressMessages(require(proxy))
suppressMessages(require(zoo))
suppressMessages(require(changepoint))


scores <- function (x, type = c("z", "t", "chisq", "iqr", "mad"), prob = NA, 
                    lim = NA) 
{
  if (is.matrix(x)) 
    apply(x, 2, scores, type = type, prob = prob, lim = lim)
  else if (is.data.frame(x)) 
    as.data.frame(sapply(x, scores, type = type, prob = prob, 
                         lim = lim))
  else {
    n <- length(x)
    s <- match.arg(type)
    ty <- switch(s, z = 0, t = 1, chisq = 2, iqr = 3, mad = 4)
    if (ty == 0) {
      res <- (x - mean(x))/sd(x)
      if (is.na(prob)) 
        res
      else {
        if (prob == 1) 
          pnorm(res)
        else if (prob == 0) 
          abs(res) > (n - 1)/sqrt(n)
        else abs(res) > qnorm(prob)
      }
    }
    else if (ty == 1) {
      t <- (x - mean(x))/sd(x)
      res <- (t * sqrt(n - 2))/sqrt(n - 1 - t^2)
      if (is.na(prob)) 
        res
      else {
        if (prob == 1) 
          pt(res, n - 2)
        else if (prob == 0) 
          abs(res) > (n - 1)/sqrt(n)
        else abs(res) > qt(prob, n - 2)
      }
    }
    else if (ty == 2) {
      res <- (x - mean(x))^2/var(x)
      if (is.na(prob)) 
        res
      else {
        if (prob == 1) 
          pchisq(res, 1)
        else abs(res) > qchisq(prob, 1)
      }
    }
    else if (ty == 3) {
      res <- x
      Q1 <- quantile(x, 0.25)
      Q3 <- quantile(x, 0.75)
      res[x >= Q1 & res <= Q3] <- 0
      res[x < Q1] <- (res[x < Q1] - Q1)/IQR(x)
      res[x > Q3] <- (res[x > Q3] - Q3)/IQR(x)
      if (is.na(lim)) 
        res
      else abs(res) > lim
    }
    else if (ty == 4) {
      res <- (x - median(x))/mad(x)
      if (is.na(prob)) 
        res
      else {
        if (prob == 1) 
          pnorm(res)
        else if (prob == 0) 
          abs(res) > (n - 1)/sqrt(n)
        else abs(res) > qnorm(prob)
      }
    }
  }
}

argmax <- function(x, y, w=1, ...) {
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

localPeak <- function(outputName5, outputName7,  x, y, w=2, span=0.1, peakPvalue) {
  peaks <- argmax(x, y, w=w, span=span)
  diffScore = scores(peaks$y.hat, prob = 1)
  
  if(outputName5!='NA'){
    pdf(outputName5)
    plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep=""), xlim=c(0,200))
    lines(x, peaks$y.hat,  lwd=2) #$
    y.min <- min(y)
    sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, peaks$y.hat[i]), col="Red", lty=2))
    points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25)
    points(x[peaks$i][diffScore[peaks$i]>peakPvalue], peaks$y.hat[peaks$i][diffScore[peaks$i]>peakPvalue], col="Blue", pch=19, cex=1.25)
    dev.off()
  }
  
  if(outputName7!='NA'){
    pdf(outputName7)
    plot(density(log(y+1)))
    dev.off()
  }
  
  return(list(x[peaks$i][diffScore[peaks$i]>peakPvalue], peaks$y.hat))
}

stripeCalling <- function(contactMap, upLimit, outputPath, outputName1, outputName2, startPoint, endPoint){
  # input: Hi-C matrix
  # output: sum of differetial Hi-C matrix on each line
  contactMap = as.matrix(contactMap)
  contactMap[contactMap > upLimit] = upLimit
  contactMapUpperTri = contactMap
  contactMapUpperTri[lower.tri(contactMapUpperTri, diag = TRUE)] = 0
  diffMap = diff(contactMapUpperTri)
  vectorDiffMap = as.vector(diffMap)
  quantileCutoff = quantile(vectorDiffMap, 0.9)     # set a cutoff for the differential value to be considered
  diffMap[abs(diffMap) < quantileCutoff] = 0
  binDiffMap = diffMap
  binDiffMap[binDiffMap<0] = -1
  binDiffMap[binDiffMap>0] = 1
  
  if(startPoint!=0&endPoint!=0){
    contactMapUpperTriMean = colMeans(t(binDiffMap[startPoint:endPoint, startPoint:endPoint]))
  }
  else{
    contactMapUpperTriMean = colMeans(t(binDiffMap))
  }
  
  
  bcpMean <- bcp(contactMapUpperTriMean,mcmc = 10000)
  
  if(outputName1!='NA'){
    if(startPoint!=0&endPoint!=0){
      s <- subplot(
        plot_ly(y = contactMapUpperTriMean, type = "bar"), 
        plot_ly(z = t(diff(contactMapUpperTri)[startPoint:endPoint, startPoint:endPoint]), type = "heatmap", 
                zauto = F, zmin = -30, zmax = 20),
        nrows = 2, heights = c(0.5, 0.5), 
        shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
      )}
    else{
      s <- subplot(
        plot_ly(y = contactMapUpperTriMean, type = "bar"), 
        plot_ly(z = t(diff(contactMapUpperTri)), type = "heatmap", zauto = F, zmin = -30, zmax = 20),
        nrows = 2, heights = c(0.5, 0.5), 
        shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
      )
    }
    orca(s, outputName1)
  }
  
  if(outputName2!='NA'){
    pdf(outputName2, width = 4.5, height = 4)
    plot(bcpMean)
    dev.off()
  }
  return(bcpMean)
}

oneSideData <- function(bcpMean){
  # input: sum of diff Hi-C contact on each line
  # output: up peaks and down peaks
  rawData = bcpMean$data[,2]
  posteriorMeanLowProb = bcpMean$data[bcpMean$posterior.prob<0.1, 2]
  posteriorMeanLowProb = posteriorMeanLowProb[!is.na(posteriorMeanLowProb)]
  dens = density(posteriorMeanLowProb)
  sampleData= dens$x[dens$y>(max(dens$y)/2)]
  meanLowProb = mean(posteriorMeanLowProb)
  
  n1 = length(bcpMean$data[bcpMean$data[,2] < meanLowProb, 2])
  if(length(sampleData[sampleData<meanLowProb]) ==0){
    sampleDownSide = sample(meanLowProb, n1, replace = TRUE)
  }
  else{
    sampleDownSide = sample(sampleData[sampleData<meanLowProb], n1, replace = TRUE)
  }
  
  n2 = length(bcpMean$data[bcpMean$data[,2] >= meanLowProb, 2])
  if(length(sampleData[sampleData >= meanLowProb])==0){
    sampleUpSide = sample(meanLowProb, n2, replace = TRUE)
  }
  else{
    sampleUpSide = sample(sampleData[sampleData >= meanLowProb], n2, replace = TRUE)
  }
  
  nMissValue = length(which(rawData==0))
  
  upSideData = rawData
  upSideData[upSideData<meanLowProb] = sampleDownSide
  upSideData[which(rawData==0)-1] = sample(sampleData, nMissValue, replace = TRUE)
  upSideData[which(rawData==0)+1] = sample(sampleData, nMissValue, replace = TRUE)
  
  downSideData = rawData
  downSideData[downSideData>=meanLowProb] = sampleUpSide
  downSideData[which(rawData==0)-1] = sample(sampleData, nMissValue, replace = TRUE)
  downSideData[which(rawData==0)+1] = sample(sampleData, nMissValue, replace = TRUE)
  
  ########### use estimated TAD size to detect the down peaks ##############
  #localPeak(outputName5, 1:length(upSideData), upSideData)  ########## local peak ##########
  #localPeak(outputName6, 1:length(downSideData), -downSideData)  ########## local peak ##########
  
  bcpUpSide = bcp(upSideData)
  bcpDownSide = bcp(downSideData)
  
  return(list(bcpUpSide, bcpDownSide))
}

matchPeaks <- function(upSide, downSide, p, map, outputName5, outputName6, outputName7){
  upLocalPeak = localPeak(outputName5, outputName7, x= upSide$data[,1],y= upSide$posterior.mean, peakPvalue = 0.9)  # local peak 
  downLocalPeak = localPeak(outputName6, outputName7, x = downSide$data[,1], y = -downSide$posterior.mean, peakPvalue = 0.8)  # local peak 
  
  upSideData <- upSide$data
  downSideData <- downSide$data
  upPeak = upSideData[upSide$posterior.prob>p, ]   # get up peaks
  downPeak = downSideData[downSide$posterior.prob>p, ]   # get down peaks
  
  if( !all(is.na(upLocalPeak[[1]])) & !all(is.na(downLocalPeak[[1]])) & is.matrix(upPeak)){
    # cluster up-peaks and estimate the cutoff for distance between adjacent up peaks and down peaks
    boundaryDist = proxy::dist(upPeak[-nrow(upPeak),1], method = "euclidean")
    boundaryDist = as.matrix(boundaryDist)
    diagonal = boundaryDist[row(boundaryDist) == col(boundaryDist)+1]
    
    if(length(diagonal)>2){
      cluster = kmeans(diagonal, 2)
      meanTadSize = max(cluster$centers)   # mean value for TAD size
    }
    else{
      meanTadSize = nrow(map)
    }
    
    stripe = data.frame()
    for(i in 1:length(upLocalPeak[[1]])){
      dis.tmp = min(abs(upLocalPeak[[1]][i]-downLocalPeak[[1]]))
      if(dis.tmp<0.3*meanTadSize){
        index = which.min(abs(upLocalPeak[[1]][i]-downLocalPeak[[1]]))
        stripe = rbind(stripe, c(upLocalPeak[[1]][i], downLocalPeak[[1]][index]))
      }
    }
    
    if(nrow(stripe)>0){
      
      colMin = data.frame(apply(stripe, 1, FUN=min)-5)
      colMin[,2] = 1
      colMax = data.frame(apply(stripe, 1, FUN=max)+5)
      colMax[,2] = ncol(map)
      
      stripe[,3] = apply(colMin, 1, FUN = max)
      stripe[,4] = apply(colMax, 1, FUN = min)
      
      #stripe = addFC(stripe)
      colnames(stripe) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
      stripe[nrow(stripe)+1,] <- NA
      stripe = stripe[stripe[,1]<stripe[,2],]    # up peak need to be left of down peak
      return(list(stripe, 0.1*meanTadSize, upLocalPeak[[2]], downLocalPeak[[2]]))
    }
    else{
      stripe = data.frame(matrix(ncol = 4, nrow = 0))
      colnames(stripe) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
      stripe[1,] <- NA
      return(list(stripe, 1, upLocalPeak[[2]], downLocalPeak[[2]]))
    }
    
  }
  else{
    stripe = data.frame(matrix(ncol = 4, nrow = 0))
    colnames(stripe) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
    stripe[1,] <- NA
    return(list(stripe, 1, upLocalPeak[[2]], downLocalPeak[[2]]))
  }
}

reduceStripe <- function(stripe, matchCutoff){
  if(nrow(stripe)>2){
    stripeReduce = data.frame(matrix(ncol = 4, nrow = 0))
    colnames(stripeReduce) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
    
    flag = FALSE
    stripeTmp = data.frame(matrix(ncol = 4, nrow = 0))
    colnames(stripeTmp) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
    stripeTmp[1, ] = stripe[1, ]
    
    for(i in 2:(nrow(stripe)-1)){
      if(abs(stripe[i,1]-stripe[(i-1),2]) < matchCutoff){
        if(flag){
          stripeTmp[1,] <- c(stripeTmp[1,1], stripe[i,2], stripeTmp[1,3], stripe[i,4])
          colnames(stripeTmp) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
        }
        else{
          stripeTmp[1,] <- c(stripe[(i-1),1], stripe[i,2], stripe[(i-1),3], stripe[i,4])
          colnames(stripeTmp) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
        }
        flag = TRUE
      }
      else{
        stripeReduce = rbind(stripeReduce, stripeTmp)
        stripeTmp = data.frame(matrix(ncol = 4, nrow = 0))
        colnames(stripeTmp) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
        stripeTmp[1, ] = stripe[i, ]
        flag = FALSE
      }
    }
    stripeReduce = rbind(stripeReduce, stripeTmp)
    stripeReduce[nrow(stripeReduce)+1,]<-NA
    return(stripeReduce)
  }
  else{
    return(stripe)
  }
}

getStripe <- function(contactMap, upLimit, outputPath='NA', outputName1='NA', 
                     outputName2='NA', outputName3='NA', outputName4='NA',outputName5='NA', outputName6='NA', outputName7='NA', 
                     startPoint=0, endPoint=0, p=0.9){
  # contactMap: a matrix for Hi-C contact map
  # upLimit: a number of up limit of Hi-C contact
  # outputPath: the output file path 
  # outputName1: the file name for heatmap of differential contact map
  # outputName2: the file name for bcp plot for all peaks
  # outputName3: the file name for bcp plot for up peaks
  # outputName4: the file name for bcp plot for down peaks
  # outputName5: the file name for up local peak
  # outputName6: the file name for down local peak
  # outputName7: the file name for density of up local peak
  # startPoint: the start point for the sub-contact map
  # endPoint: the end point for the sub-contact map
  # p: p-value cutoff for peak calling
  # output: a list of stripe, probability of up peak, probability of down peak
  
  if(outputPath!='NA'){
    setwd(outputPath)
  }

  bcpLeft = stripeCalling(contactMap, upLimit, outputPath, outputName1, outputName2, startPoint, endPoint)
  oneSideLeft = oneSideData(bcpLeft)
  upSideLeft = oneSideLeft[[1]]
  downSideLeft = oneSideLeft[[2]]
  matchPeaks.temp = matchPeaks(upSideLeft, downSideLeft, p, contactMap, outputName5, outputName6, outputName7)
  stripeLeft = matchPeaks.temp[[1]]
  matchCutoff = matchPeaks.temp[[2]]
  
  stripeReduce = reduceStripe(stripeLeft, matchCutoff)
  
  if(outputName3!='NA'){
    pdf(outputName3,width = 4.5, height = 3.5)
    plot(upSideLeft)
    dev.off()
  }
  
  if(outputName4!='NA'){
    pdf(outputName4,width = 4.5, height = 3.5)
    plot(downSideLeft)
    dev.off()
  }
  return(list(stripeReduce, upSideLeft, downSideLeft, matchPeaks.temp[[3]], matchPeaks.temp[[4]])) # output: strip, up bcp, down bcp, up local peak, down loacal peak
}

compareStripe <- function(list1, list2){
  # list1: a list for the first stripe generating by getStripe function
  # list2: a list for the second stripe generating by getStripe function
  # output: a list of result for the first and second stripe with differential p value
  
  foldChange <- function(stripe, upData, downData, upLocal, downLocal){
    # input: stripe: stripe of the first sample
    # upData: up local peak for the second sample
    # downData: down local peak for the second sample
    stripe2 = head(stripe, -1)
    getUp <- function(location){
      return(max(upLocal[location[3]:location[4]]))
    }
    getDown <- function(location){
      return(min(downLocal[location[3]:location[4]]))
    }
    
    getPvalue <- function(location){
      return(max(upData$posterior.prob[location[3]:location[4]]))
    }
    
    stripe2[,ncol(stripe2) + 1] = apply(stripe2, 1, FUN=getUp)
    stripe2[,ncol(stripe2) + 1] =  apply(stripe2, 1, FUN=getDown)
    stripe2[,ncol(stripe2) + 1] = log(abs(stripe2[,ncol(stripe2)-1])/apply(data.frame(-stripe2[,ncol(stripe2)],0.0001),1,FUN=max))
    #stripe2[,ncol(stripe2) + 1] = log(abs(stripe2[,ncol(stripe2)-1])/abs(stripe2[,ncol(stripe2)]))
    stripe2[,ncol(stripe2) + 1] =  apply(stripe2, 1, FUN=getPvalue)
    stripe2[(nrow(stripe2)+1),] <- NA
    return((stripe2))
  }
  
  stripe1 = foldChange(list1[[1]], list1[[2]], list1[[3]], list1[[4]], -list1[[5]])
  stripe1 = foldChange(stripe1, list2[[2]], list2[[3]], list2[[4]], -list2[[5]])
  stripe2 = foldChange(list2[[1]], list2[[2]], list2[[3]], list2[[4]], -list2[[5]])
  stripe2 = foldChange(stripe2, list1[[2]], list1[[3]], list1[[4]], -list1[[5]])
  
  stripe1$pvalue = pnorm(stripe1[,7]-stripe1[,11],0,2*0.5538227^2)
  colnames(stripe1) <- c('upPeak.loc', 'downPeak.loc', 'leftEdge', 'rightEdge', 'upPeak.sample1', 
                        'downPeak.sample1', 'logFoldChange.sample1', 'stripe.pValue.sample1', 'upPeak.sample2', 
                        'downPeak.sample2', 'logFoldChange.sample2', 'stripe.pValue.sample2', 'diffStripe.pValue')
  
  stripe2$pvalue = pnorm(stripe2[,7]-stripe2[,11],0,2*0.5538227^2)
  colnames(stripe2) <- c('upPeak.loc', 'downPeak.loc', 'leftEdge', 'rightEdge', 'upPeak.sample2', 
                        'downPeak.sample2', 'logFoldChange.sample2', 'stripe.pValue.sample2', 'upPeak.sample1', 
                        'downPeak.sample1', 'logFoldChange.sample1', 'stripe.pValue.sample1', 'diffStripe.pValue')
  return(list(stripe1, stripe2))
}

stripeLocation <- function(stripe, contactMap, upLimit, j, p, outputName5 = 'NA', outputName6 = 'NA', outputName7='NA'){
  # stripeMap: contact map for the pre-stripe
  # contactMap: contact map
  # start: start for the pre-stripe
  # end: end for the pre-stripe
  # p: p value cutoff for a break point
  # output: a vector for the stripe
  
  contactMap = as.matrix(contactMap)
  contactMap[contactMap > upLimit] = upLimit
  
  stripe.left = min(stripe[j,1], stripe[j,2])
  stripe.right = max(stripe[j,1], stripe[j,2])
  
  start = max(stripe.left - max((stripe.right - stripe.left)*3, 10), 1)
  end = min(stripe.right + max((stripe.right - stripe.left)*3, 10), nrow(contactMap))
  subMap = contactMap[start:end, end:nrow(contactMap)]
  stripeRollMean = t(rollmean(t(subMap), end-start))
  stripeRollMean.diff = diff(stripeRollMean)
  stripeList = c()
  
  if(length(stripeRollMean.diff) > 0){
    for(i in 1:(ncol(subMap)-(end-start)+1)){
      bcp.temp = bcp(stripeRollMean.diff[,i])
      oneSide.temp = oneSideData(bcp.temp)
      peaks = matchPeaks(oneSide.temp[[1]], oneSide.temp[[2]], 0.8, contactMap, outputName5 = 'NA', outputName6 = 'NA', outputName7='NA')
      point = reduceStripe(peaks[[1]], peaks[[2]])
      if(nrow(point)>1){
        pointCenter = (point[,1]+point[,2])/2
        stripeList = c(stripeList, head(pointCenter, -1))
      }else{
        stripeList = c(stripeList, 0)
      }
    }
    m.binseg <- cpt.mean(stripeList, test.stat = "Normal", method = "AMOC")
    change.data = m.binseg@data.set
    change.cpts = c(0, m.binseg@cpts)
    stripeEnd = 0
    
    for(t in 1:(length(change.cpts)-1)){
      data.tmp = change.data[(change.cpts[t]+1):(change.cpts[t+1])]
      m = length(data.tmp[data.tmp>0])
      n = length(data.tmp)
      if(m/n > 0.5){
        stripeEnd = change.cpts[t+1]
      }
    }
    stripeLength = end+stripeEnd
  }else{
    stripeLength = end
  }
  return(stripeLength)
}

diffStripe <- function(contactMap1, contactMap2, upLimit1, upLimit2){
  # contactMap1: the first matrix of the Hi-C contact Map
  # contactMap2: the second matrix of the Hi-C contact Map
  # upLimit: up limitation of contact map
  # output: two lists of differential strip: list(contactMap1ToContactMap2, contactMap2ToContactMap1)

  stripe1 = getStripe(contactMap1, upLimit1)
  stripe2 = getStripe(contactMap2, upLimit2)

  compare1 = compareStripe(stripe1, stripe2)
  compare1[[1]]$direction = 'left'
  compare1[[2]]$direction = 'left'
  
  contactMap1.reverse = t(contactMap1[nrow(contactMap1):1, ncol(contactMap1):1])
  contactMap2.reverse = t(contactMap2[nrow(contactMap2):1, ncol(contactMap2):1])
  
  stripe1 = getStripe(contactMap=contactMap1.reverse, upLimit=upLimit1)
  stripe2 = getStripe(contactMap=contactMap2.reverse, upLimit=upLimit2)

  compare2 = compareStripe(stripe1, stripe2)
  compare2[[1]][,1:4] = nrow(contactMap1) - compare2[[1]][,1:4]
  compare2[[2]][,1:4] = nrow(contactMap2) - compare2[[2]][,1:4]
  compare2[[1]]$direction = 'right'
  compare2[[2]]$direction = 'right'
  
  pair1ToPair2 = rbind(compare1[[1]], compare2[[1]])
  pair1ToPair2 = pair1ToPair2[order(pair1ToPair2$upPeak.loc),]
  
  pair2ToPair1 = rbind(compare1[[2]], compare2[[2]])
  pair2ToPair1 = pair2ToPair1[order(pair2ToPair1$upPeak.loc),]
  
  addLength<-function(pair, contactMap, upLimit){
    contactMap.reverse = t(contactMap[nrow(contactMap):1, ncol(contactMap):1])
    pair$stripeLength = 0
    for(i in 1:(nrow(pair)-1)){
      if(pair[i,]$diffStripe.pValue < 0.05){
        if(pair[i,]$direction == 'left'){
          stripeL = stripeLocation(pair, contactMap, upLimit, i, 0.8)
          pair[i,]$stripeLength = stripeL
        }else{
          stripeL = stripeLocation(pair, contactMap.reverse, upLimit, i, 0.8)
          pair[i,]$stripeLength = nrow(contactMap.reverse) - stripeL  ## change
        }
      }
    }
    return(pair)
  }
  
  pair1ToPair2 = addLength(head(pair1ToPair2, -1), contactMap1, upLimit1)
  pair2ToPair1 = addLength(head(pair2ToPair1, -1), contactMap2, upLimit2)
  
  return(list(head(pair1ToPair2, -1), head(pair2ToPair1, -1)))
}


suppressMessages(library("optparse"))
options(warn=-1)

option_list = list(
  make_option(c("-f", "--inputFile"), type="character", 
              help="Two input Hi-C contact maps (separated by comma)", metavar="character"),
  make_option(c("-p", "--parameter"), type="character", 
              help="Two parameters for Hi-C contact maps (separated by comma)", metavar="character"),
  make_option(c("-o", "--output"), type="character", 
              help="Output directory name", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

inputFile = opt$inputFile
maps = unlist(strsplit(inputFile, ","))
parameter = opt$parameter
parameters = unlist(strsplit(parameter, ","))
hicMap1 = read.table(maps[1])
hicMap2 = read.table(maps[2])
out = diffStrap(hicMap1, hicMap2, as.numeric(parameters[1]), as.numeric(parameters[2]))
dir.create(file.path(getwd(), opt$output), showWarnings = FALSE)
write.table(out[[1]], file.path(getwd(), opt$output, '1.txt'), sep='\t', row.names = F)
write.table(out[[2]], file.path(getwd(), opt$output, '2.txt'), sep='\t', row.names = F)
