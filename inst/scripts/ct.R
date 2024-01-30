
bsdb_sigs <- bsdb_sigs[sort(map_int(bsdb_sigs, length), decreasing = TRUE)]
background_sigs <- names(bsdb_sigs)
bsdb_sigs[[1]]

v_sigs <- bsdb_sigs[grep("vagina", names(bsdb_sigs), value = TRUE)]
b_sigs <- background_sigs[grep("vagina", names(bsdb_sigs), value = TRUE)]

v_sig <- bsdb_sigs$vagina_genus_UP_endometriosis
bck_sig <- background_sigs$vagina_genus_UP_endometriosis
bpsig <- bp_sigs$`bugphyzz:aerophilicity|facultatively anaerobic|genus`
bpsig2 <- bp_sigs$`bugphyzz:aerophilicity|aerobic|genus`



ct <- .contingencyTable(set = v_sig, reference = bck_sig, sig = bpsig2)
ct[1]

ct[1] + ct[2]

data.frame(
    bp_and_bsdb = ct[1],
    bp_and_ref = ct[2],
    no_bp_
)

tibble::tribble(
    ~ set, ~ n,
    "in bug sig and in bsdb sig" = ct[1],
    "in bug sig and not in bsdb sig" = ct[2]
)

sum(v_sig %in% bpsig2)
sum(bck_sig %in% bpsig2)


ct[1] # taxa in bsdb sig annotatied with fac -- bsdb_annotated
ct[1] + ct[2] # taxa in the bsdb sig + taxa in the bck annotatitade with fac -- background_annotated
ct[1] + ct[3] # Total size of bsdb sig -- bsdb_sig_size
sum(ct) # Total size of the background -- background_size









