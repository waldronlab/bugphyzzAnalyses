

sig1 <- sigs_by_subset_cat_s$mouth$decreased$`Smoking-behavior`
sig2 <- sigs_by_subset_cat_s$mouth$increased$`Smoking-behavior`


a <- dbEn2(control = sig1, case = sig2, term_list = bpSigs_s)
b <- dbHt(a)

rowData(a) |> as.data.frame() |> View()





# all ---------------------------------------------------------------------

sigs_by_subset_cat_s

(0/15) - (11/40)

mean(s1 %in% bpSigs_s$`bugphyzz:shape|coccus`) -
    mean(s2 %in% bpSigs_s$`bugphyzz:shape|coccus`)



x <- myOutput$matrix["shape|coccus",, drop = FALSE]



hist(calcEffectSize(s1, s2, myOutput$matrix))


wilcox.test(c(0, 0, 0, 0 ,0), c(1, 1, 1, 0, 0
                                ))
group2 <- rep(0, 15)
group1 <- c(rep(1, 11), rep(0, 29))

mean(g2) - mean(g1)
