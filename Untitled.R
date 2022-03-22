setwd("~/Desktop/NSP13/Enrich Analysis")

### Data Gathering ###

GO_BP <- read.delim("GO_Biological_Process_2021_table.txt")
GO_CC <- read.delim("GO_Cellular_Component_2021_table.txt")
GO_MF <- read.delim("GO_Molecular_Function_2021_table.txt")



### triming colems ###

GO_BP <- GO_BP [, -5:-8]
GO_CC <- GO_CC [, -5:-8]
GO_MF <- GO_MF [, -5:-8]



### constracting Rich factor from Overlap Columes ###

#############################################################

#GO_BP$Overlap <- as.numeric(GO_BP$Overlap)

RF_BP <- data.frame(a = sub("/.*", "", GO_BP$Overlap), b = sub(".*/","",GO_BP$Overlap))

RF_BP$a <- as.numeric(RF_BP$a)
RF_BP$b <- as.numeric(RF_BP$b)

RF_BP <- (RF_BP$a / RF_BP$b)
RF_BP <- as.data.frame(RF_BP)
RF_BP <- data.frame(RF = RF_BP, GeneNumber = sub("/.*", "", GO_BP$Overlap))

GO_BP <- cbind.data.frame( GO_BP, RF_BP)



################################################################


RF_CC <- data.frame(a = sub("/.*", "", GO_CC$Overlap), b = sub(".*/","",GO_CC$Overlap ))

RF_CC$a <- as.numeric(RF_CC$a)
RF_CC$b <- as.numeric(RF_CC$b)

RF_CC <- (RF_CC$a / RF_CC$b)
RF_CC <- as.data.frame(RF_CC)
RF_CC <- data.frame(RF = RF_CC,  GeneNumber = sub("/.*", "", GO_CC$Overlap))


GO_CC <- cbind.data.frame( GO_CC, RF_CC)

#GO_CC <- GO_CC [, -1]

################################################################


RF_MF <- data.frame(a = sub("/.*", "", GO_MF$Overlap), b = sub(".*/","",GO_MF$Overlap ))


RF_MF$a <- as.numeric(RF_MF$a)
RF_MF$b <- as.numeric(RF_MF$b)

RF_MF <- (RF_MF$a / RF_MF$b)
RF_MF <- as.data.frame(RF_MF)
RF_MF <- data.frame(RF = RF_MF, GeneNumber = sub("/.*", "", GO_MF$Overlap))

GO_MF <- cbind.data.frame( GO_MF, RF_MF)

#GO_MF <- GO_MF [, -1]

################################################################





### adjusting P valuse ###
#GO_Biological_process <- data.frame(Term = GO_Biological_process$Term, adjPValue = GO_Biological_process$Adjusted.P.value, PValue = GO_Biological_process$P.value ,RF = GO_Biological_process$Overlap, Genes = GO_Biological_process$Genes )
GO_BP <- subset(GO_BP, Adjusted.P.value < 0.05 )
GO_CC <- subset(GO_CC, Adjusted.P.value < 0.05 )
GO_MF <- subset(GO_MF, Adjusted.P.value < 0.05 )

dim(GO_BP)
dim(GO_CC)
dim(GO_MF)


### TOP 10 Trimming ###

GO_BP <- GO_BP [-12:-108, ]
GO_CC <- GO_CC [-12:-13,]
GO_MF <- GO_MF [-12:-15, ]

dim(GO_BP)
dim(GO_CC)
dim(GO_MF)


#install.packages('greta')
#RF <- calculate(GO_Biological_process$Overlap)
### Creating GO Enrich Plot Utlizing ggplot2 Library ###

library(ggplot2)
pdf("Results/GO_BP.pdf", width = 16, height = 12)
ggplot(GO_BP, aes(x = RF_BP, y = Term)) +geom_point(aes(color= Adjusted.P.value, size=GeneNumber)) +scale_color_gradientn(colors = rainbow(5)) +labs(x = "Rich Factor", y = NULL, color = "Adjusted.P.value", size="Gene Number") + theme(axis.title = element_text(face= "bold") ,axis.text = element_text(face = "bold"))
dev.off()

pdf("Results/GO_CC.pdf", width = 16, height = 12)
ggplot(GO_CC, aes(x = RF_CC, y = Term)) +geom_point(aes(color= Adjusted.P.value, size=GeneNumber)) +scale_color_gradientn(colors = rainbow(5)) +labs(x = "Rich Factor", y = NULL, color = "Adjusted.P.value", size="Gene Number") + theme(axis.title = element_text(face= "bold") ,axis.text = element_text(face = "bold"))
dev.off()


pdf("Results/GO_MF.pdf", width = 16, height = 12)
ggplot(GO_MF, aes(x = RF_MF, y = Term)) +geom_point(aes(color= Adjusted.P.value, size=GeneNumber)) +scale_color_gradientn(colors = rainbow(5)) +labs(x = "Rich Factor", y = NULL, color = "Adjusted.P.value", size="Gene Number") + theme(axis.title = element_text(face= "bold") ,axis.text = element_text(face = "bold"))
dev.off()

# NO RUN #

