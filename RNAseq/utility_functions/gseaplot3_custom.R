library(RColorBrewer)

custom_gseaplot2 <- function(x, geneSetID) {
  statistics <- x %>% as.data.frame() %>% filter(Description == geneSetID)
  
  gsdata <- enrichplot:::gsInfo(x, geneSetID)
  
  # statistics label
  signif_label <- ifelse(nrow(x@result) > 1,
                         paste0("NES: ", round(statistics$NES,2),"\np-value: ", format(statistics$pvalue, digits=2, scipen = 0), "\nq-value: ",  format(statistics$qvalues, digits=2, scipen = 0)),
                         paste0("NES: ", round(statistics$NES,2),"\np-value: ", format(statistics$pvalue, digits=2, scipen = 0))
  )
  
  
  enrichmentScore <- x@result[geneSetID, "enrichmentScore"]
  
  es.df <- data.frame(es = which.min(abs(gsdata$runningScore - enrichmentScore)))
  # es.df <- data.frame(es = which.min(abs(p$data$runningScore - enrichmentScore)))
  
  # setup theme
  p1 <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_void(8) +
    coord_cartesian(clip = "off") +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.33, lineend="butt"),
          axis.ticks.length=unit(.07, "cm"),
          axis.line = element_line(color = "black", linewidth = 0.33, lineend="square"),
          axis.text = element_text(size = 8, color = "black"),
          axis.title.x = element_text(color = "black", size = 10),
          plot.margin = margin(t=0, b=0,unit="cm"),
          axis.ticks.margin = unit(0, "null"),
          text = element_text(size = 8),
          legend.position = "none") +
    scale_x_continuous(expand=c(0,0)) +
    geom_segment(data=es.df, aes_(x = ~es, xend = ~es, y = 0, yend=enrichmentScore),
                 colour = "#2479ae", linetype = "dashed", linewidth = 0.33) +
    # geom_segment(data=es.df, aes_(x = -Inf, xend = ~es, y = enrichmentScore, yend=enrichmentScore),
    #              colour = "#2479ae", linetype = "dashed", linewidth = 0.33) +
    geom_hline(yintercept = 0, color = "#a9aaaa", linetype = "dashed", linewidth = 0.33) +
    # add area
    # geom_polygon(aes_(y = ~runningScore),
    #              size=1, fill="#1A9F62", alpha = 0.3) +
    # add line
    geom_line(aes_(y = ~runningScore, color= ~Description),
              size=2/.pt, color = "#2479ae") +
    # position for cohordinates based on positive or negative enrichment score
    {if(enrichmentScore > 0)annotate("text",
                                     x=nrow(gsdata)*0.97, y=enrichmentScore*0.97,
                                     hjust = "right", vjust = "top",size=3,
                                     label = signif_label)} +
    {if(enrichmentScore < 0)annotate("text",
                                     x=nrow(gsdata)*0.03, y=enrichmentScore*0.97,
                                     hjust = "left", vjust = "bottom",size=3,
                                     label = signif_label)} +
    ylab("Enrichment Score (ES)") +
    ggtitle(geneSetID)
  
  # now banded plot ---------------------------------------
  
  p2 <- ggplot() +
    geom_vline(data = gsdata %>% filter(position == 1), aes(xintercept=x), color="black", size = 0.3, alpha = 0.5) +
    scale_x_continuous(expand = c(0,0), limits = c(min(gsdata$x), max(gsdata$x))) +
    scale_y_continuous(expand=c(0,0))  +
    coord_cartesian(clip = "off") +
    xlab(NULL) + ylab(NULL) + theme_void(8) +
    theme(legend.position = "none",
          #plot.margin = margin(t=-.1, b=0,unit="cm"),
          plot.margin = margin(t=0, b=0,unit="cm"),
          plot.background = element_blank(),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.ticks.margin = unit(0, "null"),
          panel.border = element_blank()) +
    xlab("Gene Rank")
  
  
  p1/p2 + plot_layout(heights = unit(c(2.13, 0.1), "cm"), widths = unit(3.3, "cm"))
}
