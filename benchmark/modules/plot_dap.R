

plot_dap <- function(x,
                     grid_nrow = 2, 
                     grid_ncol = 1, 
                     label_size = 2,
                     top_rank = 5,
                     lim_prob = c(0, 1.2),
                     ...)
{
  label_size_config = label_size
  label_size_snp = label_size
  top_rank_config = top_rank
  top_rank_snp = top_rank
  lim_prob_config = lim_prob
  lim_prob_snp = lim_prob
    
  p2 <- plot_set(x,  
    top_rank = top_rank_config, 
    label_size = label_size_config, 
    lim_prob = lim_prob_config, ...)
  p3 <- plot_snp(x, 
    top_rank = top_rank_snp,
    label_size = label_size_snp, 
    lim_prob = lim_prob_snp, ...)
  
  plot_grid(p2, p3,  labels = "AUTO", nrow = grid_nrow, ncol = grid_ncol)
}


plot_set <- function(x, lim_prob, label_size, top_rank, ...)
{
  ptab <- x$set

  ptab <- head(ptab, top_rank)

  ptab <- mutate(ptab,
    label = paste0(snp, "\n", 
      "P = ", round(cluster_prob, 2),
      "; ", "avg(r^2) = ", round(cluster_avg_r2, 2)))

  ggplot(ptab, aes(cluster_prob, cluster)) + 
    geom_vline(xintercept = 1, linetype = 3) + 
    geom_point() + 
    geom_segment(aes(xend = cluster_prob, yend = cluster, x = 0)) + 
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.025, size = label_size) + 
    xlim(lim_prob) + 
    scale_y_continuous(limits  = c(min(top_rank, nrow(ptab)) + 0.5, 0.5), trans = "reverse")
}


plot_snp <- function(x, lim_prob, label_size, top_rank, ...)
{
  ptab <- x$snp
  
  ptab <- head(ptab, top_rank)

  ptab <- mutate(ptab,
    rank = seq(1, n()), 
    label = paste0(snp, "\n", 
      "P = ", round(snp_prob, 2),
      "; ", "log10(BF) = ", round(snp_log10bf, 2)))

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
    print(plot_dap(result[[r]], top_rank = top_rank))
}
dev.off()
