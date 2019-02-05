ewce.plot.mod <- function (total_res, mtc_method = "bonferroni") 
{
  if (!mtc_method %in% c("holm", "hochberg", "hommel", "bonferroni", 
                         "BH", "BY", "fdr", "none")) {
    stop("ERROR: Invalid mtc_method argument. Please see '?p.adjust' for valid methods.")
  }
  multiList = TRUE
  if (is.null(total_res$list)) {
    multiList = FALSE
  }
  total_res$q = p.adjust(total_res$p, method = mtc_method)
  ast_q = rep("", dim(total_res)[1])
  ast_q[total_res$q < 0.05] = "*"
  total_res$ast_q = ast_q
  total_res$sd_from_mean[total_res$sd_from_mean < 0] = 0
  graph_theme = theme_bw(base_size = 12, base_family = "Helvetica") + 
    theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
          axis.line = element_line(size = 0.7, color = "black"), 
          text = element_text(size = 14), axis.title.y = element_text(vjust = 0.6))
  upperLim = max(abs(total_res$sd_from_mean))
  total_res$y_ast = total_res$sd_from_mean * 1.05
  total_res$abs_sd = abs(total_res$sd_from_mean)
  if ("Direction" %in% colnames(total_res)) {
    the_plot = ggplot(total_res) + geom_bar(aes_string(x = "CellType", 
                                                       y = "abs_sd", fill = "Direction"), position = "dodge", 
                                            stat = "identity") + graph_theme
  }
  else {
    the_plot = ggplot(total_res) + geom_bar(aes_string(x = "CellType", 
                                                       y = "abs_sd"), fill = "red", stat = "identity") + 
      graph_theme + theme(legend.position = "none")
  }
  the_plot = the_plot + theme(plot.margin = unit(c(1, 0, 0, 
                                                   0), "mm"), axis.text.x = element_text(angle = 55, hjust = 1)) + 
    theme(panel.border = element_rect(colour = "black", 
                                      fill = NA, size = 1)) + xlab("") + theme(strip.text.y = element_text(angle = 0)) + 
    coord_cartesian(ylim = c(0, 1.1 * upperLim)) + ylab("Std.Devs. from the mean") + 
    theme(plot.margin = unit(c(0, 0, 0, 1.5), "cm"))
  the_plot = the_plot + scale_y_continuous(breaks = c(0, ceiling(upperLim * 
                                                                   0.66))) + geom_text(aes_string(label = "ast_q", x = "CellType", 
                                                                                                  y = "y_ast"), size = 10)
  if (multiList) {
    the_plot = the_plot + facet_wrap("~ list", scales = "free")
  }
  return(the_plot)
}