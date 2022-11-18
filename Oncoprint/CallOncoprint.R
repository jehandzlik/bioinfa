# Oncoprints for WTC U01 Melanoma Grant 

#setwd("R:/Kirchhofflab/homes/Kirchhofflabspace/Project-World-Trade-Center_U01/Oncoprint/")

setwd("/Users/handzj01/Desktop/World-Trade-Center-U01-211112/Oncoprint/Output/")


library(ComplexHeatmap)

####### Oncoprint of NYU Driver Genes in NYU Pts ####### 
###### Top 50 Genes (NYU NYU)  ########
#mat = read.table("/Users/handzj01/Desktop/World-Trade-Center-U01-211112/Oncoprint/Input/WTC_U01_pt_NYU_Driver_50_df_forOncoprint.txt", header = TRUE,stringsAsFactors=FALSE, sep = "\t")
mat = read.table("/Users/handzj01/Desktop/World-Trade-Center-U01-211112/Oncoprint/Input/NYU_Driver_50_first_filter.txt", header = TRUE,stringsAsFactors=FALSE, sep = "\t")

mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]

##### Define colors for Alterations 

# Synonymous SNV = SYN_SNV, "#D8D32C"
# Nonsynonymous SNV = NONSYN_SNV,  "#008000"
# Splicing = SPLICE, "#87456B"
# Stopgain = STOPGAIN, "#37E1BF"
# Stoploss = STOPLOSS,  "#726991"
# Frameshift Substitution = FRAME_SUB,  "red"
# Nonframeshift Substitution = NONFRAME_SUB,   "blue"
# Frameshift deletion = FRAME_DEL, "#D8D32C"
# Nonframeshift deletion = NONFRAME_DEL, "#8b7774"
# Frameshift insertion = FRAME_INS, "#55c575"
# Unknown = UNKNOWN, "#A00303"

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

column_title = "OncoPrint for profile of top \nNYU PACT Driver Genes (n=50) \nin WTC U01 Melanoma (n=9 Samples)"
heatmap_legend_param = list(title = "Alterations", at = c("NONSYN_SNV", "SPLICE","STOPGAIN", "STOPLOSS","FRAME","NONFRAME"), 
                            labels = c("nonsynonymous SNV", "splicing", "stopgain", "stoploss",
                                       "frameshift (sub, ins, del)","nonframeshift (sub, ins, del)"))

pdf("newSamplesFirstFilter.pdf", height = 15, width = 10)
oncoPrint(mat,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)
dev.off()
