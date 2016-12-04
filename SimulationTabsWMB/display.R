library(reshape)
library(ggplot2)
saveRDS(wmise,"ABMmise.rds")
saveRDS(omise,"OUmise.rds")
saveRDS(gmise,"GBMmise.rds")
saveRDS(cmise,"CIRmise.rds")
#TABLE
#out1 <- list(wmise[[1]], omise[[1]])
#out2 <- list(gmise[[1]], cmise[[1]])
#compvar <- function(mat){
 #   mat[,ncol(mat)] <- mat[,1] - mat[,2]
 #   return (mat)
#}
#getout <- function(out, outfile){
#    out <- lapply(out,compvar)
#    out<-do.call(cbind,out)
#    num <- function(x){
#        return (paste("\\num[scientific-notation=true,round-precision=3,round-mode=figures]{",x, "} &"))
#    }
#    out <- apply(out, c(1,2), num)
#    out[1,6] <- gsub("&", "\\\\\\\\", out[1,6])
#    out[2,6] <- gsub("&", "\\\\\\\\", out[2,6])
#    out[3,6] <- gsub("&", "\\\\\\\\", out[3,6])
#    cat(t(out),file=outfile)
#}
#getout(out1,"out1.txt")
#getout(out2,"out2.txt")

#FIGURES

#X<-data.frame(wmise[[2]])
#dl <- melt(X, id = "Times")
#pa <- ggplot(data=dl, aes(x=Times, y=value, colour=variable)) + geom_line()+ theme_bw() + theme(text=element_text(size=30), plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank() ,panel.border = element_blank())+ theme(axis.line = element_line(color = 'black'), legend.position="bottom", legend.title=element_blank())+ggtitle("Fourier Process") + ylab("Spot variance")
#X<-data.frame(omise[[2]])
#dl <- melt(X, id = "Times")
#po <- ggplot(data=dl, aes(x=Times, y=value, colour=variable)) + geom_line()+ theme_bw() + theme(text=element_text(size=30),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank() ,panel.border = element_blank())+ theme(axis.line = element_line(color = 'black'),legend.position="bottom", legend.title=element_blank())+ 
#ggtitle("Ornstein-Uhlenbeck (OU)") + ylab("Spot (variance)")
X<-data.frame(gmise[[2]])
dl <- melt(X, id = "Times")
pg <- ggplot(data=dl, aes(x=Times, y=value, colour=variable)) + geom_line()+ theme_bw() + theme(text=element_text(size=25),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank() ,panel.border = element_blank())+ theme(axis.line = element_line(color = 'black'),legend.position="bottom", legend.title=element_blank())+ ggtitle("Geometric Brownian Motion (GBM)") + ylab("Spot (variance)")
#X<-data.frame(cmise[[2]])
#dl <- melt(X, id = "Times")
#pc <- ggplot(data=dl, aes(x=Times, y=value, colour=variable)) + geom_line()+ theme_bw() + theme(text=element_text(size=30),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank() ,panel.border = element_blank())+ theme(axis.line = element_line(color = 'black'),legend.position="bottom", legend.title=element_blank())+ggtitle("Cox-Ingersoll-Ross (CIR)")+ylab("Spot (variance)")
#
#pdf("pw.pdf")
#print(pa)
#dev.off()
#pdf("po.pdf")
#print(po)
#dev.off()
pdf("pg.pdf")
print(pg)
dev.off()
#pdf("pc.pdf")
#print(pc)
#dev.off()
