
plot_finemap <- function(x,
                         grid_nrow = NULL, 
                         grid_ncol = NULL, 
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
  lim_prob_ncausal = lim_prob   
    
  p1 <- plot_ncausal(x, 
    lim_prob = lim_prob_ncausal, ...)
  p2 <- plot_set(x,  
    top_rank = top_rank_config, 
    label_size = label_size_config, 
    lim_prob = lim_prob_config, ...)
  p3 <- plot_snp(x, 
    top_rank = top_rank_snp,
    label_size = label_size_snp, 
    lim_prob = lim_prob_snp, ...)
  
  plot_grid(p1, p2, p3,  labels = "AUTO", nrow = grid_nrow, ncol = grid_ncol)
}


plot_ncausal <- function(x, lim_prob, ...)
{
  ptab <- x$ncausal
  
  sum_prop_zero <- filter(ptab, ncausal_num == 0)[["prob"]]  %>% sum
  if(sum_prop_zero == 0) {
    ptab <- filter(ptab, ncausal_num != 0)
  }
  
  ptab <- mutate(ptab, 
    ncausal_num = factor(ncausal_num, levels = sort(unique(ncausal_num), 
                                                    decreasing = TRUE)),
    type = factor(type, levels = c("prior", "post")))
    
  p <- ggplot(ptab, aes(ncausal_num, ncausal_prob, fill = type)) + 
    geom_hline(yintercept = 1, linetype = 3) + 
    geom_bar(stat = "identity", position = "dodge") + 
    coord_flip() + theme(legend.position = "top") + 
    scale_fill_manual(values = c("grey50", "orange")) +
    ylim(lim_prob)
    
  return(p)
}

plot_set <- function(x, lim_prob, label_size, top_rank, ...)
{
  ptab <- x$set

  ptab <- head(ptab, top_rank)

  ptab <- mutate(ptab,
    label = paste0(config, "\n", 
      "P = ", round(config_prob, 2),
      "; ", "log10(BF) = ", round(config_log10bf, 2)))

  ggplot(ptab, aes(config_prob, rank)) + 
    geom_vline(xintercept = 1, linetype = 3) + 
    geom_point() + 
    geom_segment(aes(xend = config_prob, yend = rank, x = 0)) + 
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.025, size = label_size) + 
    xlim(lim_prob) + 
    scale_y_continuous(limits  = c(top_rank + 0.5, 0.5), trans = "reverse")
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
    print(plot_finemap(result[[r]], top_rank = top_rank))
}
dev.off()
