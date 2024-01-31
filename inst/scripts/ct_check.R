target_set <- bsdb_sigs$feces_species_DOWN_psychosis
back_set <- background_sigs$feces_species_DOWN_psychosis
an <- bp_sigs$`bugphyzz:habitat|free-living|species`



ct <- .contingencyTable(set = target_set, reference = back_set, sig = an)

## CT checked
# sum(target_set %in% an)
# as.integer(target_set %in% an)
# sum(back_set %in% an)

ct <- ct + 0.5
epitools::oddsratio.wald(x = ct + 0.5)$measure[2,]
(or <- (ct[1] / ct[2]) / (ct[3] / ct[4]))
(lower_ci <- exp(log(or) - 1.96 * sqrt(sum(1 / (ct + 0.5) ))))
(upper_ci <- exp(log(or) + 1.96 * sqrt(sum(1 / (ct + 0.5) ))))

(p_value <- stats::fisher.test(ct, alternative = 'g')$p.value)


