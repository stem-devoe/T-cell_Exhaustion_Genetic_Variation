my.theme = theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), # add border to whole plot # 0.2
                 axis.line = element_line(linewidth = 0.5), # 0.2
                 axis.title = element_text(size = 6, margin = margin(0,0,0,0)),
                 axis.text = element_text(size = 6, margin = margin(0,0,0,0)), # tickmark labels
                 plot.title = element_text(size = 7, margin = margin(0,0,0,0)), # plot title
                 axis.ticks.length = unit(0.05, "cm"), # adjust length of tick (outside panel)
                 axis.ticks =  element_line(linewidth = 0.5)) # 0.2

xaxis.gap.fix = scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
yaxis.gap.fix = scale_x_continuous(expand = expansion(mult = c(0, 0.05)))

x_at_zero = scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0,0.05)))
                               
set_linewidth = 0.5