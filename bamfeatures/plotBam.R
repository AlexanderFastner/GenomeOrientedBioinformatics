#read data in
#RPKM distribution Y across genes X
#1.all reads
#2.PCRindex = 0 reads
#add legend


# load Data
install.packages("data.table")
library("data.table")
install.packages("ggplot2")
library("ggplot2")
#do math in java, then output only necessary data into file
#RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
all<- fread("C:\\Users\\alex\\gobi2021\\bamfeatures\\nookaew_cm\\RPKMall")

nonpcr<-fread("C:\\Users\\alex\\gobi2021\\bamfeatures\\nookaew_cm\\RPKMPCRO")

all$type<-"allreads"
nonpcr$type<-"nonpcr"

plotting <-rbind(all, nonpcr)

p<-ggplot(plotting,aes(V1,fill = type)) + geom_density(alpha = 0.2) + scale_x_log10() + xlab("RPKM (log-scale)") + ggtitle("nookaew_cm")

png("nookaew_cm.png")
print(p)
dev.off()




