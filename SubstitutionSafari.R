############ SUBSTITUTION SAFARI - SUBSTITUTION SPECTRA FOR EVERYONE ##############

## Last Update - 25/01/2022 ##
## mrs72 / Max Stammnitz ##
## maxrupsta@gmail.com ##
## University of Cambridge  ##

library(vcfR)
library(Biostrings)
library(stringr)
library(scales)

## 1. Functions
substitution.spectrum <- function(x, normalised){
  
  # a. Isolate triplets form Platypus VCF
  x[,'TRIPLET'] <- as.character(subseq(x = reference[as.character(x[,'CHROM'])], 
                                       start = as.numeric(x[,'POS']) - 1, 
                                       end = as.numeric(x[,'POS']) + 1))
  
  # b. Isolate base changes
  context <- x[,'TRIPLET']
  context.changes <- matrix(NA, nrow = length(context), ncol = 3)
  colnames(context.changes) <- c("REF", "ALT-5'", "ALT-3'")
  context.changes[,1] <- context
  context.changes[,2:3] <- str_split_fixed(context.changes[,1],"",3)[,c(1,3)]
  context.changes[,2] <- paste(context.changes[,2],as.character(x[,'ALT']),context.changes[,3],sep="")
  context.changes <- context.changes[,c(1,2)]
  colnames(context.changes) <- c("REF", "ALT")
  
  # Pyrimidine-context substitutions
  equivalent.mut <- matrix(c('ACA>AAA', 'TGT>TTT', # AC - AA, or TG - TT ###  C>A or G>T
                             'ACC>AAC', 'GGT>GTT',
                             'ACG>AAG', 'CGT>CTT',
                             'ACT>AAT', 'AGT>ATT',
                             'CCA>CAA', 'TGG>TTG', # CC - CA, or GG - GT
                             'CCC>CAC', 'GGG>GTG',
                             'CCG>CAG', 'CGG>CTG',
                             'CCT>CAT', 'AGG>ATG',
                             'GCA>GAA', 'TGC>TTC', # GC - GA, or CG - CT
                             'GCC>GAC', 'GGC>GTC',
                             'GCG>GAG', 'CGC>CTC',
                             'GCT>GAT', 'AGC>ATC',
                             'TCA>TAA', 'TGA>TTA', # TC - TA, or AG - AT
                             'TCC>TAC', 'GGA>GTA',
                             'TCG>TAG', 'CGA>CTA',
                             'TCT>TAT', 'AGA>ATA',
                             'ACA>AGA', 'TGT>TCT', # AC - AG, or TG - TC ###  C>G or G>C
                             'ACC>AGC', 'GGT>GCT',
                             'ACG>AGG', 'CGT>CCT',
                             'ACT>AGT', 'AGT>ACT',
                             'CCA>CGA', 'TGG>TCG', # CC - CG, or GG - GC
                             'CCC>CGC', 'GGG>GCG',
                             'CCG>CGG', 'CGG>CCG',
                             'CCT>CGT', 'AGG>ACG',
                             'GCA>GGA', 'TGC>TCC', # GC - GG, or CG - CC
                             'GCC>GGC', 'GGC>GCC',
                             'GCG>GGG', 'CGC>CCC',
                             'GCT>GGT', 'AGC>ACC',
                             'TCA>TGA', 'TGA>TCA', # TC - TG, or AG - AC
                             'TCC>TGC', 'GGA>GCA',
                             'TCG>TGG', 'CGA>CCA',
                             'TCT>TGT', 'AGA>ACA',
                             'ACA>ATA', 'TGT>TAT', # AC - AT, or TG - TA ###  C>T or G>A
                             'ACC>ATC', 'GGT>GAT',
                             'ACG>ATG', 'CGT>CAT',
                             'ACT>ATT', 'AGT>AAT',
                             'CCA>CTA', 'TGG>TAG', # CC - CT, or GG - GA
                             'CCC>CTC', 'GGG>GAG',
                             'CCG>CTG', 'CGG>CAG',
                             'CCT>CTT', 'AGG>AAG',
                             'GCA>GTA', 'TGC>TAC', # GC - GT, or CG - CA
                             'GCC>GTC', 'GGC>GAC',
                             'GCG>GTG', 'CGC>CAC',
                             'GCT>GTT', 'AGC>AAC',
                             'TCA>TTA', 'TGA>TAA', # TC - TT, or AG - AA
                             'TCC>TTC', 'GGA>GAA',
                             'TCG>TTG', 'CGA>CAA',
                             'TCT>TTT', 'AGA>AAA',
                             'ATA>AAA', 'TAT>TTT', # AT - AA, or TA - TT ###  T>A or A>T
                             'ATC>AAC', 'GAT>GTT',
                             'ATG>AAG', 'CAT>CTT',
                             'ATT>AAT', 'AAT>ATT',
                             'CTA>CAA', 'TAG>TTG', # CT - CA, or GA - GT
                             'CTC>CAC', 'GAG>GTG',
                             'CTG>CAG', 'CAG>CTG',
                             'CTT>CAT', 'AAG>ATG',
                             'GTA>GAA', 'TAC>TTC', # GT - GA, or CA - CT
                             'GTC>GAC', 'GAC>GTC',
                             'GTG>GAG', 'CAC>CTC',
                             'GTT>GAT', 'AAC>ATC',
                             'TTA>TAA', 'TAA>TTA', # TT - TA, or AA - AT
                             'TTC>TAC', 'GAA>GTA',
                             'TTG>TAG', 'CAA>CTA',
                             'TTT>TAT', 'AAA>ATA',
                             'ATA>ACA', 'TAT>TGT', # AT - AC, or TA - TG ### T>C or A>G
                             'ATC>ACC', 'GAT>GGT',
                             'ATG>ACG', 'CAT>CGT',
                             'ATT>ACT', 'AAT>AGT',
                             'CTA>CCA', 'TAG>TGG', # CT - CC, or GA - GG
                             'CTC>CCC', 'GAG>GGG',
                             'CTG>CCG', 'CAG>CGG',
                             'CTT>CCT', 'AAG>AGG',
                             'GTA>GCA', 'TAC>TGC', # GT - GC, or CA - CG
                             'GTC>GCC', 'GAC>GGC',
                             'GTG>GCG', 'CAC>CGC',
                             'GTT>GCT', 'AAC>AGC',
                             'TTA>TCA', 'TAA>TGA', # TT - TC, or AA - AG
                             'TTC>TCC', 'GAA>GGA',
                             'TTG>TCG', 'CAA>CGA',
                             'TTT>TCT', 'AAA>AGA',
                             'ATA>AGA', 'TAT>TCT', # AT - AG, or TA - TC ### T>G or A>C
                             'ATC>AGC', 'GAT>GCT',
                             'ATG>AGG', 'CAT>CCT',
                             'ATT>AGT', 'AAT>ACT',
                             'CTA>CGA', 'TAG>TCG', # CT - CG, or GA - GC
                             'CTC>CGC', 'GAG>GCG',
                             'CTG>CGG', 'CAG>CCG',
                             'CTT>CGT', 'AAG>ACG',
                             'GTA>GGA', 'TAC>TCC', # GT - GG, or CA - CC
                             'GTC>GGC', 'GAC>GCC',
                             'GTG>GGG', 'CAC>CCC',
                             'GTT>GGT', 'AAC>ACC',
                             'TTA>TGA', 'TAA>TCA', # TT - TG, or AA - AC
                             'TTC>TGC', 'GAA>GCA',
                             'TTG>TGG', 'CAA>CCA',
                             'TTT>TGT', 'AAA>ACA'),
                           nrow=96, ncol=2, byrow=TRUE)
  
  # c. Define variants from input table and convert them to pyrimidine context
  triplett_snvs <- paste0(context.changes[,1], '>', context.changes[,2])
  ind.convert <- match(triplett_snvs,equivalent.mut[,2])
  triplett_snvs[!is.na(ind.convert)] <- equivalent.mut[ind.convert[!is.na(ind.convert)],1]
  
  # d. Count triplets
  counts <- table(triplett_snvs)
  consensus.mut <- equivalent.mut[,1]
  if(length(counts)<96){
    add.names <- consensus.mut[which(is.na(match(consensus.mut,names(counts)))==T)]
    names.takeover <- names(counts)
    counts <- c(counts,rep(0,length(add.names)))
    names(counts) <- c(names.takeover,add.names)
  }
  counts <- counts[match(consensus.mut,names(counts))]
  
  # e. normalisation to genome triplet counts
  if(normalised == T){
    
    reference_trinucleotides.div <- c(rep(as.numeric(reference_trinucleotides[,5])[1:16],3),
                                      rep(as.numeric(reference_trinucleotides[,5])[17:32],3))
    counts.norm <- c(counts/reference_trinucleotides.div)/sum(c(counts/reference_trinucleotides.div))
    names(counts) <- names(counts.norm) <- consensus.mut
    out <- list("counts" = counts, "counts.normalised" = counts.norm)
  
  } else if (normalised == F){
    
    names(counts) <- names(counts.norm) <- consensus.mut
    out <- list("counts" = counts, "counts.normalised" = counts.norm)
    
  }
  
  # f. Output
  return(out)
}
plot.substition.spectrum <- function(x, title, peak.colour){
  
  ## setup plot background colours and blocks
  mar.default <- c(2,4,2,2) + 0.1
  par(mar = mar.default + c(3, 11, 3, -2))
  colors = c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
  mut <- c('C > A','C > G','C > T','T > A','T > C','T > G')
  y.top = 0.35
  borders = c(0, 115.4*1/6, 115.4*2/6, 115.4*3/6, 115.4*4/6, 115.4*5/6, 115.4)
  alphas = rep(0.2,6)
  plot(1, type="n", ylim=c(0,y.top), xlim=c(0, 115.4), xlab="", ylab="", axes=F,
       main = title,
       cex.main = 4)
  for (i in 1:6) {
    rect(xleft=borders[i], xright=borders[i+1], ybottom=0, ytop=y.top, col=alpha(colors[i], alphas[i]), border="white")
    rect(xleft=borders[i], xright=borders[i+1], ybottom=y.top-c(y.top*0.12), ytop=y.top, col=colors[i], border="white")
    text(x=(borders[i]+borders[i+1])/2, y=y.top-c(y.top*0.12)/2, labels=mut[i], cex=5, col="white")
  }
  
  ## add spectrum
  out <- barplot(x/sum(x),
                 col = rep(peak.colour,96),
                 border = NA,
                 ylab = "",
                 yaxt = 'n',
                 ylim = c(0, 0.35),
                 add = T,
                 names.arg = NA)
  
  ## add title
  title(ylab="Substitutions [%]", line = 9, cex.lab = 6)
  
  # add axes
  axis(side = 1, 
       las = 2,
       at = out[,1],
       labels = c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A",
                  "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C",
                  "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G",
                  "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T",
                  "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A",
                  "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C",
                  "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G",
                  "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
                  "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A",
                  "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", "A[T>A]C",
                  "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G",
                  "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T",
                  "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "A[T>C]A",
                  "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C",
                  "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G",
                  "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
                  "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A",
                  "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C",
                  "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G",
                  "T[T>G]T"),
      tick = F, line = -1.6, font = 1, family = 'mono', cex.axis = 1.2)
  axis(side = 2, las = 2, cex.axis = 4, col = 'black', col.axis = 'black',
       at = c(0, 0.1, 0.2, 0.3), labels = c('0', '10', '20', '30'), lwd = 5, pos=-2, hadj = 1.4)
  
}

## 2. Sequence of steps

### 2.a Load reference fasta
reference <- readDNAStringSet('my_reference.fa.gz')

### 2.b Load/preprocess the four essential substition VCF columns:
# 1.) CHROMOSOME -> 'CHROM'
# 2.) POSITION -> 'POS'
# 3.) REFERENCE SEQUENCE -> 'REF'
# 4.) ALTERNATE SEQUENCE -> 'ALT
input.vcf <- read.vcfR('my_substitutions.vcf.gz')
input.vcf <- cbind(input.vcf@fix[,c('CHROM', 'POS', 'REF', 'ALT')], 'TRIPLET'= NA)

### 2.c Load referene tri-nucleotide counts for normalisation
reference_trinucleotides <- read.table('my_reference_trinucleotides.txt', header = T)

### 2.c Generate the substitution spectrum
spectrum.out <- substitution.spectrum(x = input.vcf, normalised = T)

### 2.d Plot the substitution spectrum
pdf("my_first_substitution_spectrum.pdf", height = 12, width = 24)
plot.substition.spectrum(x = spectrum.out$counts.normalised, 
                         title = 'My first substitution spectrum!', 
                         peak.colour = 'cornflowerblue')
dev.off()
