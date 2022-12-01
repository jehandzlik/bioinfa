library(ComplexHeatmap)

#################################################
## Script for preparing input for OncoPrint and 
## generating oncoprints 
#################################################

## Filtering the annotated VCF files
filterInput <- function(input){
  
  ## Take only 0 or . from gnomAD_exome_ALL AND 0 or . from gnomAD_genome_ALL OR 
  ## gnomAD_exome_ALL == 1 AND gnomAD_genome_ALL == 1
  inputA <- input[which( ((input$gnomAD_exome_ALL == "." & input$gnomAD_genome_ALL == 0) |
                            (input$gnomAD_exome_ALL == 0 & input$gnomAD_genome_ALL == ".")) | 
                           input$gnomAD_exome_ALL == 1 & input$gnomAD_genome_ALL == 1),]
  
  ## For the ones that have . AND . in gnomAD_genome_ALL AND gnomAD_exome_ALL and 
  ## use two sorts: 
  ## 1) for Kaviar_AF we take everything <= 5*10-5,
  ## 2) for the Kaviar_AF that have . we only take those that have
  ## cosmic94 
  
  input.int <- input[which(input$gnomAD_exome_ALL == "." & input$gnomAD_genome_ALL == "."),]
  
  inputB <- input.int[which(input.int$Kaviar_AF != "."), ]
  inputB <- inputB[which(as.numeric(inputB$Kaviar_AF) <= 5*10^-5 ), ]
  
  inputC <- input.int[which(input.int$Kaviar_AF == "." & input.int$cosmic94 != "."),]
  
  input <- rbind(inputA, inputB, inputC)
  
  ## Take everything in range [0.05 - 0.68] from FREQ
  input <- input[which(input$FREQ <= 0.68 & input$FREQ >= 0.05),]
  
  ## Exclude chromosomes X and Y
  input <- input[which(input$CHR != "chrX" | input$CHR != "chrY"),]
  
  return(input)
  
}

#################################################
## Prepare input for Oncroprint
#################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) ## you might need to change that if you're not using RStudio

data <- read.delim("VCF-LoFreq-annot.all.txt", header = T, stringsAsFactors = F, sep="\t")

data <- filterInput(data) ## Adjust to your own needs. Also make sure that 
## Mutect2 results are not filtered here or use different filter for those.

#subset to only relevant cols to avoid confusion
# Study ID, Func.refGene, Gene.refGene, ExonicFunc.refGene
data=data[,c(2,10,11,13)] ## adjust to your input!

colnames(data) <- c("Study_ID", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene")

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

gene_list = read.table("gene_list.txt",
                       header = TRUE,stringsAsFactors=FALSE, sep = "\t")

sample=unique(input$Study_ID)
for (i in 1:length(unique(input$Study_ID))){
  
  temp=input[input$Study_ID==sample[[i]],] #subset to one patient at a time
  temp2=aggregate(ExonicFun.new~ Gene.refGene, temp, paste,collapse=";") # collapsing duplicate genes and print exonic function separated by ";"
  gene_list[[i+1]]=with(temp2, ExonicFun.new[match(gene_list$Genes,Gene.refGene)]) #use gene.ref as key and match to print exonic function
  
}

colnames(gene_list)=c("Genes",paste0("Sample_",sample))

mat <- gene_list

#write.table(mat,file="oncoprint.input.txt",col.names = T,row.names = F,quote = F,sep = "\t")

#################################################
## Call Oncroprint
#################################################

mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]

##### Graphics for alterations

col = c("NONSYN_SNV" = "#726991", 
        "SPLICE" = "#F3A358", "STOPGAIN" = "#37E1BF", "STOPLOSS" = "#D8D32C",
        "FRAME" = "red",  "NONFRAME" = "blue")

#### Defining alteration colors and widths

alter_fun = list(
  background = alter_graphic("rect", width=0.9, height=0.9, fill = "#CCCCCC"),   
  NONSYN_SNV = alter_graphic("rect", width=0.9, height=0.9, fill = col["NONSYN_SNV"]),
  SPLICE = alter_graphic("rect", width=0.9, height=0.6, fill = col["SPLICE"]),
  STOPGAIN = alter_graphic("rect", width=0.9, height=0.4, fill = col["STOPGAIN"]),
  STOPLOSS = alter_graphic("rect", width=0.9, height=0.3, fill = col["STOPLOSS"]),
  FRAME = alter_graphic("rect", width=0.9, height=0.2, fill = col["FRAME"]),
  NONFRAME = alter_graphic("rect", width=0.9, height=0.2, fill = col["NONFRAME"]))

#### Add title and Legends
column_title = "OncoPrint for profile of top \nTCGA Driver Genes  \nin Melanoma patients" ## change the title of the plot
heatmap_legend_param = list(title = "Alterations", at = c("NONSYN_SNV", "SPLICE","STOPGAIN", "STOPLOSS","FRAME","NONFRAME"), 
                            labels = c("nonsynonymous SNV", "splicing", "stopgain", "stoploss",
                                       "frameshift (sub, ins, del)","nonframeshift (sub, ins, del)"))


col.ord <- match(sort(colnames(mat)), colnames(mat))


pdf("Oncoprint.pdf", height = 12, width = 10)
oncoPrint(mat,
          column_order = col.ord,
          show_column_names = TRUE,
          alter_fun = alter_fun, col = col,
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)

dev.off()
