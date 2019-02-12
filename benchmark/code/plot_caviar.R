plot_caviar <- function(x,
                        grid_nrow = NULL, 
                        grid_ncol = NULL, 
                        label_size = 2,
                        top_rank = 5,
                        lim_prob = c(0, 1.5),
                        ...)
{
  plot_snp(x, label_size, top_rank, lim_prob, ...)
}

plot_snp <- function(x, label_size, top_rank, lim_prob, ...)
{
  ptab <- x$snp

  ptab <- head(ptab, top_rank)

  ptab <- mutate(ptab,
    label = paste0(snp, "\n", 
      "P = ", round(snp_prob, 2),
      "; ", "P(set) = ", round(snp_prob_set, 2)))

  ggplot(ptab, aes(snp_prob, rank)) +
    geom_vline(xintercept = 1, linetype = 3) + 
    geom_point() + 
    geom_segment(aes(xend = snp_prob, yend = rank, x = 0)) + 
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.025, size = label_size) + 
    xlim(lim_prob) + 
    scale_y_continuous(limits  = c(top_rank + 0.5, 0.5), trans = "reverse")
}

png(plot_file)
for (r in 1:length(result)) {
    print(plot_caviar(result[[r]], top_rank = top_rank))
}
dev.off()
