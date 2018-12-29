library(aphid) #The 'aphid' package for analysis with profile hidden Markov models
library(ape) #The 'ape' package for Analyses of Phylogenetics and Evolution
library(seqinr) #The 'seqinr' package for Biological Sequences Retrieval and Analysis
library(magicfor) #Functions to Obtain Results from for Loops
library(readxl) #read xls files
library(caret) #Classification and Regression Training


          ##Model Building##

#read seed alignment sequences from FASTA file - sequences unaligned
trainglobin = read.fasta(file = "globin-training.fasta", 
                         seqtype = "AA", set.attributes = FALSE)

#store sequences as AAbin object 
seqs = as.AAbin.list(trainglobin)

#Multiple Sequence Alignment (MSA) using Hidden Markov Models
algn = align(seqs, type = "global", residues = "AMINO")

#Derive profile from MSA
globinPHMM = derivePHMM(algn, residues = "AA")

#Export model
writePHMM(globinPHMM, file = "Globin.hmm")

       ##Test Model##

#Test a single protein sequence if it is Globin 

proteinToCheck = as.AAbin("MVLSPADKTNVKAAWGKVGAHAGEYGAEAL
                          ERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNAL
                          SALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR")

Viterbi(globinPHMM, proteinToCheck)

test = read.fasta(file = "testing-main.fasta", 
                 seqtype = "AA", set.attributes = FALSE)

magic_for(print,silent = TRUE)

for (i in test){
  res = Viterbi(globinPHMM, i)
  print(res[["score"]])
} 

scoreddt = magic_result_as_dataframe()

#read predicted dataset with column Score and Class empty to be filled
globin.predicted = read_excel("globin-predicted.xls")

#fill Score column
globin.predicted$Score = scoreddt 

#fill the Class column with yes or no 
globin.predicted$Class = ifelse(globin.predicted$Score>0,"Yes","No")

#read actual dataset with column Class filled with actual data. Globin in this case is yes rest are no
globin.actual <- read_excel("globin-actual.xls")

confusionMatrix(globin.predicted$Class, globin.actual$Class) 

     ##Plot Model##

plotg = read.fasta(file = "Globin-members.fasta", 
                   seqtype = "AA", set.attributes = FALSE)

magic_for(print,silent = TRUE)

for (i in plotg){
  res = Viterbi(globinPHMM, i)
  print(res[["score"]])} 

globinsp = magic_result_as_dataframe() 

write.table(globinsp, "globinplotscores.txt", sep="\t") #export scores as text for plotting

plotting = read_excel("globinplot.xls")

Len = plotting$Length

scores = plotting$Score

plot.default(Len, scores, xlab = "Protein Sequence Length", ylab = "Log-odds Score", main = "Globin Family")