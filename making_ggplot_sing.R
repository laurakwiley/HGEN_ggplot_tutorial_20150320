library(ggplot2)
library(dplyr)

##################
# Data Processing
##################
man.data<-original
# Read in Data
#man.data<-tbl_df(read.table(file="plink.assoc.txt",header=TRUE))
man.data<-tbl_df(read.csv(file="data.csv"))
man.data<-man.data %.% na.omit()
man.data<-man.data[sample(1:nrow(man.data),40000,replace=FALSE),]
man.data$CHR<-as.factor(man.data$CHR)


# Sort by CHR,BP
man.data<-arrange(man.data, CHR, BP)
  
# Add in index position for plot
man.data$position<-c(1:nrow(man.data))

##################
# Manhattan Plot
##################

m.plot<-ggplot(man.data, aes(x=position, y=-log10(P))) + geom_point(aes(color=man.data$CHR))
m.plot

# Add horizontal line
m.plot <- m.plot + geom_hline(y=-log10(1e-5),color="green")
m.plot

# Remove legend
m.plot <- m.plot + theme(legend.position="none")
m.plot

## Change colors of datapoints
cols<-rep(c("black","blue"),11)
m.plot <- m.plot + scale_colour_manual(values=cols)
m.plot

# Make x-axis labels by chromosome
# Need tick position and label name
axis_labels<-man.data %.% group_by(CHR) %.% summarize(mean(position))
names(axis_labels)<-c("CHR","position")

m.plot<-m.plot + scale_x_continuous(breaks = axis_labels$position, labels = axis_labels$CHR)
m.plot


# Prettify label positions
axis_labels$CHR<-as.numeric(axis_labels$CHR)
axis_labels$label<-as.character(axis_labels$CHR)
axis_labels[axis_labels$CHR%%2==0,]$label <- paste("\n",axis_labels[axis_labels$CHR%%2==0,]$CHR,sep="")
m.plot<-m.plot + scale_x_continuous(breaks = axis_labels$position, labels = axis_labels$label)
m.plot


# Prettify Axis labels
m.plot <- m.plot + xlab("Chromosome") + ylab(expression(-log[10](p)))
m.plot



## Annotation
# Select out points we want
annotate<-man.data[man.data$P<1e-5,]
m.plot <- m.plot + geom_point(data=annotate,colour="darkgreen")
m.plot

m.plot + geom_text(data=annotate,aes(label=SNP))

m.plot <- 
  m.plot + geom_text(data=annotate,colour="darkgreen",aes(label=SNP,vjust=-.2,hjust=1.1))
m.plot


##################
# Mito Manhattan 
##################