library(ggplot2)

# load data
growthglc <- read_csv("~/github/citramalate/kinetic/result/growthglc.csv")

# XCH_GLC vs GLC_feed
ggplot(growthglc, aes(x=GLC_feed, y=XCH_GLC)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks=seq(0,1,0.1)) +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    theme_bw() +
    theme(
        plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
    ) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))

# various species concentrations vs GLC_feed

library(tidyverse)
#library(dplyr)

growthglc_collapsed <- growthglc %>%
    select(GLC_feed, GLCx, GLCp, G6P) %>%
    gather(key = "variable", value = "value", -GLC_feed)

ggplot(growthglc, aes(x=GLC_feed)) +
    geom_point(aes(y=GLCx/1000), colour = "red") +
    geom_line(aes(y=GLCx/1000), colour = "red") +
    #geom_point(aes(y=GLCp), colour = "blue") +
    #geom_line(aes(y=GLCp), colour = "blue") +
    geom_point(aes(y=G6P/1000), colour = "black") +
    geom_line(aes(y=G6P/1000), colour = "black") +
    scale_x_continuous(breaks=seq(0,1,0.1)) +
    scale_y_continuous(name="Species concentration (1,000x)",breaks=seq(0,25,2)) +
    labs(colour="g") +
    scale_colour_manual(labels=c("a", "b"),values=c("red", "blue")) +
    #scale_y_log10() +
    theme_bw() +
    theme(
        plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
    ) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))