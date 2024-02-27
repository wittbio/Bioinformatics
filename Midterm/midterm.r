#Loading in the packages Biostrings, msa, and seqinr
library(Biostrings)
library(msa)
library(seqinr)
#Setting the working directory
setwd("C:/Users/phoen/Documents/GitHub/Bioinformatics/Midterm")
#Reading the fasta file for our sequences
sequences <- readDNAStringSet("sequences.fasta")
#Printing the sequences
print(sequences)
#Aligning the sequence with ClustalW
alignment <- msa(sequences, "ClustalW")
#Converting the alignment to sequinr format
distance_alignment <- msaConvert(alignment, type = "seqinr::alignment")
#Performing a distance alignment
dm <- dist.alignment(distance_alignment)
#Printing distance alignment
dm
#Showed that the sequence with the greatest difference was #6
#Observing the differences between 6 and 7. This demonstrates that there is a deletion and a substitution in #6
sequences$Homo_sapiens_6
sequences$Homo_sapiens_7
#The gene is Homo sapiens hbb gene for beta globin. The accession number is LC121775.
#Putting the sequence for Homo_sapiens_6 into a DNAString
seq6 <- DNAString(sequences$Homo_sapiens_6)
#Translating the DNAString into a protein file using Biostrings
seq6t <- Biostrings::translate(seq6)
#Printing the newly translated protein
print(seq6t)
#Writing the translated protein to a .fasta file
write.fasta(sequences=seq6t,names = names(seq6t),file.out="seq6t.fasta")
#I used blastp to search using this sequence. The closest result was "hemoglobin subunit beta [Homo sapiens]"
#The accession number is KAI2558340
#The diseases associated with this gene are sickle-cell disease and beta thalassemia
#I ran the sequence through UniProt and according to the variant viewer this variant has beta Thalassemia.
#The variant ID is VAR_002856

