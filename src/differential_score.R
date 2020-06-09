stripeP <- function(stripe){
  u1 = stripe$upPeak.sample1
  u2 = stripe$upPeak.sample2
  m = min(u1, u2)
  if(m<0){
    u1 = u1-1.1*m
    u2 = u2-1.1*m
  }
  x1 = log(u1)-log(u2)
  m1 = mean(x1)
  sd1 = sd(x1)
  
  d1 = stripe$downPeak.sample1
  d2 = stripe$downPeak.sample2
  m = max(d1,d2)
  if(m>0){
    d1 = d1-1.1*m
    d2 = d2-1.1*m
  }
  x2 = log(-d1)-log(-d2)
  m2 = mean(x2)
  sd2 = sd(x2)
  
  x1.norm = (x1-m1)/sd1
  x2.norm = (x2-m2)/sd2
  chiq = x1.norm^2+x2.norm^2
  stripe$chiq = chiq
  stripe$p.chip = 1-pchisq(chiq,1)
  
  stripe$up.tscore = (x1-m1)/(sd1)
  stripe$p.up = 2*pt(-abs(stripe$up.tscore),df=nrow(stripe)-1)
  stripe$down.tscore = (x2-m2)/(sd2)
  stripe$p.down = 2*pt(-abs(stripe$down.tscore),df=nrow(stripe)-1)
  
  return(stripe)
}

suppressMessages(library("optparse"))
options(warn=-1)

option_list = list(
  make_option(c("-f", "--inputFile"), type="character",  # need to change
              help="Input stripe file", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

inputFile = opt$inputFile
stripe = read.table(inputFile)
stripe.diff = stripeP(stripe)

write.table(stripe.diff, file.path(getwd(), opt$output, '1.txt'), sep='\t', row.names = F) # need to change
