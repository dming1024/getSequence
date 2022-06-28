#
#获取融合基因，断点处，上下游的序列信息
#examples
#Rscript getFusionsSeq.R ENST00000373886,ENST00000358487,60553375,123239535,500
#args1: 基因1转录本ID
#args2: 基因2转录本ID
#args3: 基因1断点位置
#args4: 基因2断点位置
#args5: 断点上下游长度

vars <- commandArgs(trailingOnly = TRUE)
split.vars <- unlist(strsplit(vars, ","))
g1 = as.character(split.vars[1])
g2 = as.character(split.vars[2])
b1 = as.numeric(split.vars[3])
b2 = as.numeric(split.vars[4])
wth = as.numeric(split.vars[5])
m = load("db/EnsemblGene.Transcripts.RData")#染色体位置，start-end, exon数目，exon位置信息: 21139

getSeq = function(ensemblID, breakpoint, wth) {
  breakpoint = as.numeric(breakpoint)
  tmp = EnsemblGene.Structures[EnsemblGene.Structures$EnsemblGene == ensemblID, ]
  N = which(EnsemblGene.Structures$EnsemblGene == ensemblID)
  
  exons = list()
  exons_end = as.numeric(unlist(strsplit(as.character(tmp$exonEnd), ",")))
  exons_start = as.numeric(unlist(strsplit(as.character(tmp$exonStart), ",")))
  #这个是正义链获取坐标的形式
  if (tmp$Strand == "+") {
    #获取各exons坐标
    for (i in 1:as.numeric(as.character(tmp$exonCount))) {
      exons[[i]] = c(exons_start[i], exons_end[i])
      
      if (exons[[i]][1] - 1 < breakpoint &
          exons[[i]][2] + 1 > breakpoint) {
        #在该exons上的长度
        exon_length = as.numeric(breakpoint) - exons[[i]][1] + 1
        #在那个exons
        exon_location = i
      }
      
    }
    #计算断点前长度
    if (exon_location == 1) {
      sequence_all = exon_length
    } else{
      exon_location = exon_location - 1
      sequence_length = 0
      for (i in 1:exon_location) {
        sequence_length = sequence_length + (exons[[i]][2] - exons[[i]][1] + 1)
      }
      sequence_all = sequence_length + exon_length
    }
  }
  
  if (tmp$Strand == "-") {
    #获取各exons坐标
    exons_start_x = rev(exons_end)
    exons_end_x = rev(exons_start)
    for (i in 1:as.numeric(as.character(tmp$exonCount))) {
      exons[[i]] = c(exons_start_x[i], exons_end_x[i])
      
      if (exons[[i]][1] + 1 > breakpoint &
          exons[[i]][2] - 1 < breakpoint) {
        #在该exons上的长度
        exon_length = exons[[i]][1] - breakpoint
        #在那个exons
        exon_location = i
      }
      
    }
    #计算断点前长度
    if (exon_location == 1) {
      sequence_all = exon_length
    } else{
      exon_location = exon_location - 1
      sequence_length = 0
      for (i in 1:exon_location) {
        sequence_length = sequence_length + (exons[[i]][1] - exons[[i]][2] + 1)
      }
      sequence_all = sequence_length + exon_length
    }
  }
  
  if (wth > 0) {
    sequences_select = substr(sequences[N], sequence_all - wth, sequence_all)
  }
  if (wth < 0) {
    sequences_select = substr(sequences[N], sequence_all, sequence_all - wth)
  }
  
  return(sequences_select)
}

g1seq = getSeq(g1, b1, wth)
g2seq = getSeq(g2, b2, -wth)
fusions = paste(g1seq, g2seq, sep = "--")
writeLines(paste0(">",g1,"_",g2,"\n",fusions), "fusionSeq.txt")

#examples
#Rscript getFusionsSeq.R ENST00000373886,ENST00000358487,60553375,123239535,500
#
