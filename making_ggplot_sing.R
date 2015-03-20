library("ggplot2")
library("dplyr")
library("magrittr")

##################
# Data Processing
##################

# Read in Data
#man.data<-tbl_df(read.table(file="plink.assoc.txt",header=TRUE))
man.data <- read.csv(file="HGEN_ggplot_tutorial_20150320/data.csv") %>% tbl_df()
man.data %<>% 
  na.omit() %>% 
  sample_n(40000) %>% 
  mutate(CHR=as.factor(CHR))


# Sort by CHR,BP
man.data %<>% arrange(CHR, BP)
  
# Add in index position for plot
man.data %>% mutate(position=1:length(SNP))

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
axis_labels <- man.data %.% group_by(CHR) %.% summarize(mean(position))
names(axis_labels) <- c("CHR","position")

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

mitodata <- read.table(file="HGEN_ggplot_tutorial_20150320/mitodata.txt", header=T) %>% tbl_df()

# Sets gene names for bp ranges
addgenelabel <- function(bp,gene) { 
  gene <- ifelse(bp < 577,gene <- "Control-Region", 
               ifelse(bp < 648,gene <- "tRNA", 
               ifelse(bp < 1602,gene <- "rRNA", 
               ifelse(bp < 1671,gene <- "tRNA", 
               ifelse(bp < 3230,gene <- "rRNA", 
               ifelse(bp < 3305,gene <- "tRNA", 
               ifelse(bp < 3307,gene <- "Non-Coding", 
               ifelse(bp < 4263,gene<- "ND1", 
               ifelse(bp < 4332,gene <- "tRNA", 
               ifelse(bp < 4401,gene <- "tRNA", 
               ifelse(bp < 4402,gene <- "Non-Coding", 
               ifelse(bp < 4470,gene <- "tRNA", 
               ifelse(bp < 5512,gene <- "ND2", 
               ifelse(bp < 5580,gene <- "tRNA", 
               ifelse(bp < 5587,gene <- "Non-Coding", 
               ifelse(bp < 5656,gene <- "tRNA", 
               ifelse(bp < 5657,gene <- "Non-Coding", 
               ifelse(bp < 5730,gene <- "tRNA", 
               ifelse(bp < 5826,gene <- "tRNA", 
               ifelse(bp < 5892,gene <- "tRNA", 
               ifelse(bp < 5904,gene <- "Non-Coding", 
               ifelse(bp < 7446,gene <- "CO1", 
               ifelse(bp < 7515,gene <- "tRNA", 
               ifelse(bp < 7518,gene <- "Non-Coding", 
               ifelse(bp < 7586,gene <- "tRNA", 
               ifelse(bp < 8270,gene <- "CO2", 
               ifelse(bp < 8295,gene <- "Non-Coding", 
               ifelse(bp < 8365,gene <- "tRNA", 
               ifelse(bp < 8366,gene <- "Non-Coding", 
               ifelse(bp < 8573,gene <- "ATP8", 
               ifelse(bp < 9208,gene <- "ATP6", 
               ifelse(bp < 9991,gene <- "CO3", 
               ifelse(bp < 10059,gene <- "tRNA", 
               ifelse(bp < 10405,gene <- "ND3", 
               ifelse(bp < 10470,gene <- "tRNA", 
               ifelse(bp < 10767,gene <- "ND4L", 
               ifelse(bp < 12138,gene <- "ND4", 
               ifelse(bp < 12207,gene <- "tRNA", 
               ifelse(bp < 12266,gene <- "tRNA", 
               ifelse(bp < 12337,gene <- "tRNA", 
               ifelse(bp < 14149,gene <- "ND5", 
               ifelse(bp < 14674,gene <- "ND6", 
               ifelse(bp < 14743,gene <- "tRNA", 
               ifelse(bp < 14747,gene <- "Non-Coding", 
               ifelse(bp < 15888,gene <- "CYB", 
               ifelse(bp < 15954,gene <- "tRNA", 
               ifelse(bp < 15956,gene <- "Non-Coding", 
               ifelse(bp < 16024,gene <- "tRNA", 
               ifelse(bp < 17000,gene <- "Control-Region")))))))))))))))))))))))))))))))))))))))))))))))))}

# Add gene names to each SNP
mitodata %<>% mutate(gene=addgenelabel(bp,gene))

# Display internal structure of mitodata
mitodata %>% str()


mitodata %<>% 
  mutate(neglogp = -log10(p)) %>% # Creates and stores negative log p as a new variable
  mutate(neglogpline = -log10(0.05)) %>% # Adds sig threshold line at -log(0.05)
  mutate(extraline = -3) # Adds -3 label to y axis

# Set colors for each gene
colours <- c("Control-Region" <- "lightblue4", "tRNA" <- "magenta4", "rRNA" <- "mediumaquamarine", "Non-Coding" <- "sienna4", "ND1" <- "magenta", "ND2" <- "mediumblue", "CO1" <- "olivedrab", "CO2" <- "orange2", "ATP8" <- "orchid4", "ATP6" <- "red3", "CO3" <- "royalblue2", "ND3" <- "palegreen4", "ND4L" <- "grey0", "ND4" <- "pink4", "ND5" <- "yellow4", "ND6" <- "steelblue4", "CYB" <- "tan","red")

# Create gene boundaries and lines
visibleboundaries <- c(1,576,1601,3229,4262,5511,7445,8269,9207,9990,10404,10766,12137,14148,14673,15887)

bdries <- data.frame(x = visibleboundaries,y=-.5)

bdries$gene <- addgenelabel(bdries$x,bdries$gene)

lines <- data.frame(x = seq(0,16567,by=1),y = 0)

lines$gene <- addgenelabel(lines$x,lines$gene)

# Plot everything and GO
ggplot(mitodata, aes(x = bp,y = neglogp,color = gene)) +
  geom_point() + 
  coord_polar(direction = -1) +
  geom_line(aes(x,1.30,color = "red"),data = lines) +
  #facet_grid(.~pheno) +
  geom_line(aes(y=extraline)) +
  geom_point(aes(x,y,color = gene),data=lines) +
  scale_colour_manual(values = colours,"Genes",
                      breaks = c("Control-Region","tRNA","rRNA","Non-Coding","ND1","ND2","CO1","CO2","ATP8","ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB"),
                      labels = c("Control Region","tRNA","rRNA","Non-Coding","ND1","ND2","CO1","CO2","ATP8","ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB")) +
  xlab("Mitochondrial Base-Pair Location") +
  ylab("-log(p-value)") +
  ggtitle("Negative Log P-value of Mitochondrial Hits") +
  layer(geom="text",mapping =aes(x,y,label = x),data = bdries,size=2.5)

ggsave("solarplot.png", w=6, h=6, dpi=110)

