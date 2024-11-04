library(RColorBrewer)

custom_gseaplot2 <- function(x, geneSetID) {
  statistics <- x %>% as.data.frame() %>% filter(Description == geneSetID)
  
  gsdata <- enrichplot:::gsInfo(x, geneSetID)
  
  # statistics label
  signif_label <- ifelse(nrow(GSEA_YAPTAZ@result) > 1,
                         paste0("NES: ", round(statistics$NES,2),"\np-value: ", format(statistics$pvalue, digits=2, scipen = 0), "\nq-value: ",  format(statistics$qvalues, digits=2, scipen = 0)),
                         paste0("NES: ", round(statistics$NES,2),"\np-value: ", format(statistics$pvalue, digits=2, scipen = 0))
  )
  
  
  enrichmentScore <- x@result[geneSetID, "enrichmentScore"]
  
  es.df <- data.frame(es = which.min(abs(gsdata$runningScore - enrichmentScore)))
  # es.df <- data.frame(es = which.min(abs(p$data$runningScore - enrichmentScore)))
  
  # setup theme
  p1 <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(8) +
    theme(panel.border = element_rect(linewidth = 0.5, color="black", fill=NA),
          # panel.grid.major = element_line(colour = "grey92", linewidth = 0.5),
          #panel.grid.minor = element_line(colour = "grey92"),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 8, color = "black", margin = margin(t=0, r=0, b=0, l=0, unit="cm")),
          plot.margin = margin(t=0, b=0,unit="cm"),
          axis.ticks.length.x = unit(0, "null"),
          text = element_text(size = 8),
          legend.position = "none") +
    scale_x_continuous(expand=c(0,0)) +
    geom_segment(data=es.df, aes_(x = ~es, xend = ~es, y = 0, yend=enrichmentScore),
                 colour = "red", linetype = "dashed", linewidth = 0.5) +
    geom_segment(data=es.df, aes_(x = -Inf, xend = ~es, y = enrichmentScore, yend=enrichmentScore),
                 colour = "red", linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
    # add area
    # geom_polygon(aes_(y = ~runningScore),
    #              size=1, fill="#1A9F62", alpha = 0.3) +
    # add line
    geom_line(aes_(y = ~runningScore, color= ~Description),
              linewidth=0.7, color = "#1A9F62") +
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
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  
  v <- seq(1, sum(gsdata$position), length.out=9)
  inv <- findInterval(rev(cumsum(gsdata$position)), v)
  if (min(inv) == 0) inv <- inv + 1
  
  col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
  
  ymin <- min(gsdata$ymin)
  yy <- max(gsdata$ymax - gsdata$ymin) * 1
  xmin <- which(!duplicated(inv))
  xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
  d <- data.frame(ymin = ymin, ymax = yy,
                  xmin = xmin,
                  xmax = xmax,
                  col = col[unique(inv)])
  
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_rect(data=d, 
              aes_(xmin=~xmin, xmax=~xmax, ymin=~ymin, ymax=~ymax, fill=~I(col)),
              alpha=1, inherit.aes=FALSE) +
    geom_linerange(aes_(ymin=~ymin+ymax*0.25, ymax=~ymax*0.75), color="black", linewidth = 0.3, alpha = 0.5) +
    xlab(NULL) + ylab(NULL) + theme_classic(8) +
    theme(legend.position = "none",
          #plot.margin = margin(t=-.1, b=0,unit="cm"),
          plot.margin = margin(t=0, b=0,unit="cm"),
          plot.background = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 8, color = "black", margin = margin(t=0, r=0, b=0, l=0, unit="cm")),
          axis.ticks.length.y = unit(0, "null"),
          panel.border = element_rect(linewidth = 0.5, color="black", fill=NA)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))  +
    xlab("Gene Rank")
  
  
  p1/p2 + plot_layout(heights = unit(c(3.2, 0.3), "cm"), widths = unit(5, "cm"))
  
}
