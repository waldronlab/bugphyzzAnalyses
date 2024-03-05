
dat <- filtered_results$BSDB |>
    select(-combination) |>
    as_tibble()

dat |>
    filter(condition == "multiple myeloma") |>
    filter(direction == "increased") |>
    filter(grepl("fac", sig_name))

idx <- names(bsdb_sigs_sp) |>
    grep("bsdb:296/.*UP$", x = _)

sigList <- bsdb_sigs_sp[idx]
bpSig <- bp_sigs_sp$`bugphyzz:aerophilicity|facultatively anaerobic|species`

ti <- intersect(bpSig, sigList[[1]])
names(ti) <- taxizedb::taxid2name(ti, db = "ncbi")
ti <- ti[sort(names(ti))]
ti |> names()

taxizedb::classification(ti, db = "ncbi")
