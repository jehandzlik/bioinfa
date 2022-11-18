library("readxl")
library("xlsx")
####### Subsetting dataframes to the desired gene lists ####### 
#setwd("R:/Kirchhofflab/homes/Kirchhofflabspace/Project-World-Trade-Center_U01/Oncoprint/")

#setwd("/Volumes/Kirchhofflabspace/Project-World-Trade-Center_U01/Oncoprint/")

####### NYU Driver Genes Top 50 ####### 
#NYU_Driver_50 = read.table("Gene Lists/Gene_List_WTC_U01_pt_Driver_top50.txt", header = T,stringsAsFactors=FALSE, sep = "\t")
NYU_Driver_50 = read.table("/Users/handzj01/Desktop/World-Trade-Center-U01-211112/Oncoprint/Input/Gene_List/Gene_List_TCGA_pt_Melanoma_top50.txt", header = TRUE,stringsAsFactors=FALSE, sep = "\t")

## NYU Patients


####### Prepare input

## Concatenate separate files and create the
## data


### Cleaning and making FULL dataset, not subset to just NYU PACT Genes 
#### Removing unknown, combining all frameshift/nonframeshift
input.dir <- "/Users/handzj01/Desktop/Mutations_for_WTC_U01/FirstFilter/"
setwd(input.dir)

data <- data.frame()
for(file in list.files(input.dir)){
  
  if(length(grep("_Mutec.xlsx",file)) == 1){
    
    input <- read_excel(file)
    prt <- input[,c(3,13,14,16)]
    colnames(prt) <- c("Study_ID", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene")
    data <- rbind(data, prt)
    
  }
  
  if(length(grep("_MAF.xlsx",file)) == 1){
    
    input <- read_excel(file)
    prt <- input[,c(2,10,11,13)]
    colnames(prt) <- c("Study_ID", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene")
    data <- rbind(data, prt)
    
  }

}

##################
### OR
##################

### Create data based on collpased files

input <- read.table("/Users/handzj01/Desktop/BigPurpleSFTP/Ziyan/2022-11-WES/WES_Results/VCF-GATK-HC-annot.all.txt",
                    comment.char = "",
                    check.names = FALSE, sep = "\t", header = T)

data <- input[,c(2,10,11,13)]
colnames(data) <- c("Study_ID", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene")





#Prelim_pts_NYU <- read.delim("../WTC_Annotated_VCFs/WTC_U01_all_gene_dataframe.txt", header = T, stringsAsFactors = F, sep="\t")

#data=Prelim_pts_NYU
#subset to only relevant cols to avoid confusion
# Study ID, Func.refGene, Gene.refGene, ExonicFunc.refGene
#data=data[,c(1,7,8,10)]


data[which(data$Func.refGene == "splicing"),]
data$ExonicFunc.refGene=as.character(factor(data$ExonicFunc.refGene))
data$Func.refGene=as.character(factor(data$Func.refGene))
data$ExonicFun.new=ifelse(data$Func.refGene=="splicing","splicing",data$ExonicFunc.refGene) #keep splicing
#remove missing
data_sub=data[data$ExonicFun.new!=".",]
# Making no synonymous SNVs dataset only
data_sub_nosyn=data_sub[data_sub$ExonicFun.new!="synonymous SNV",]
data_sub_nounk=data_sub_nosyn[data_sub_nosyn$ExonicFun.new!="unknown",]
input=data_sub_nounk

# Recoding to remove spaces, make uniform uppercase 
input$ExonicFun.new[input$ExonicFun.new=="nonsynonymous SNV"] <- "NONSYN_SNV"
input$ExonicFun.new[input$ExonicFun.new=="splicing"] <- "SPLICE"
input$ExonicFun.new[input$ExonicFun.new=="stopgain"] <- "STOPGAIN"
input$ExonicFun.new[input$ExonicFun.new=="stoploss"] <- "STOPLOSS"
#input$ExonicFun.new[input$ExonicFun.new=="unknown"] <- "UNKNOWN"
input$ExonicFun.new[input$ExonicFun.new %in% c("frameshift substitution", "frameshift deletion", "frameshift insertion")] <- "FRAME"
input$ExonicFun.new[input$ExonicFun.new %in% c("nonframeshift substitution", "nonframeshift deletion", "nonframeshift insertion")] <- "NONFRAME"

gene_list=NYU_Driver_50

sample=unique(input$Study_ID)
for (i in 1:length(unique(input$Study_ID))){
  
  temp=input[input$Study_ID==sample[[i]],] #subset to one patient at a time
  temp2=aggregate(ExonicFun.new~ Gene.refGene, temp, paste,collapse=";") # collapsing duplicate genes and print exonic function separated by ";"
  gene_list[[i+1]]=with(temp2, ExonicFun.new[match(gene_list$Genes,Gene.refGene)]) #use gene.ref as key and match to print exonic function
  
}

colnames(gene_list)=c("Genes",paste0("Sample.",as.character(sample)))

write.table(gene_list,file="/Users/handzj01/Desktop/World-Trade-Center-U01-211112/Oncoprint/Input/NYU_Driver_50_first_filter.txt",col.names = T,row.names = F,quote = F,sep = "\t")

gene_list[which(gene_list$Genes == "NRAS"),]
