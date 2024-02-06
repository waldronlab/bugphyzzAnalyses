library(onlineFDR)
input <- mainDF |>
    mutate(
        p_value = case_when(
            p_value > 1 ~ 1,
            p_value < 0 ~ 0,
            TRUE ~ p_value
        )
    ) |>
    mutate(
        id = paste0(`BSDB ID`, "---", rank, "---", bp_sig),
        pval = p_value,
        # batch_name = paste0(`BSDB ID`, "---", rank)
        batch_name = paste0(`BSDB ID`)
    ) |>
    relocate(id, pval, batch_name)

b <- 1:length(unique(input$batch_name))
set.seed(1234)
names(b) <- sample(unique(input$batch_name))

input$batch <- b[match(input$batch_name, names(b))]
input <- input |>
    relocate(batch, .after = pval) |>
    arrange(batch)

result <-  BatchBH(input, alpha = 0.1)
