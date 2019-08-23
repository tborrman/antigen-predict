library(ggplot2)

df <- read.table("hausmann_predict_table.txt", header=TRUE, sep="\t")

png("hausmann_predict_plot_A6.png", width=2800, height=2400, res=600)
par(mar=c(5,5,4,2) + 0.1)
p1 <- ggplot(df, aes(x=A6_Rosetta_dG_bind, y= netMHC_affinity, col=A6_lysis)) +
  scale_colour_gradient( low="dodgerblue", high="red") +
  geom_text(aes(label=peptide), size=2) + xlim(0,135) +
  theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x= "Rosetta dG bind", y= "NetMHC Affinity (nM)", col= "% lysis", title="A6 Hausmann predicted peptides")
p1
dev.off()

png("hausmann_predict_plot_B7.png", width=2800, height=2400, res=600)
par(mar=c(5,5,4,2) + 0.1)
p2 <- ggplot(df, aes(x=B7_Rosetta_dG_bind, y= netMHC_affinity, col=B7_lysis)) +
  scale_colour_gradient( low="dodgerblue", high="red") +
  geom_text(aes(label=peptide), size=2)  + xlim(0,190) +
  theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.border = element_blank(), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black")) +
  labs(x= "Rosetta dG bind", y= "NetMHC Affinity (nM)", col= "% lysis", title="B7 Hausmann predicted peptides")
p2
dev.off()

