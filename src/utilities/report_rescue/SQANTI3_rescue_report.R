#!/bin/bash Rscript
#
#-------------------------------------------------------------------------------
#                       _________________________
#
#                         SQANTI3 RESCUE REPORT
#                       _________________________
#
# Author: Ángeles Arzalluz-Luque
# 
# Contact: angeles.arzalluz@gmail.com
#
# Affiliation: Institute for Integrative Systems Biology, CSIC, Valencia, Spain
#
# Last updated: 14/June/2022
#
#-------------------------------------------------------------------------------



####------------------------- FILTER REPORT ARG LIST -----------------------####

#### Define script arguments ####

option_list <- list(
  optparse::make_option(c("-d","--dir"), type = "character", 
                        help = "Output/input directory - must be same as SQ3 filter."),
  optparse::make_option(c("-o","--output"), type = "character", default = "SQANTI3", 
                        help = "Output report file prefix."),
  optparse::make_option(c("-f", "--filter_type"), type = "character", default = "ml",
                        help = "Type of SQ3 filter that was run to generate input.
                        Must be one of the following: 'ml' or 'rules'.
                        Default: 'ml'.")
)

# Parse arguments
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)


message("\n--------------------------------------------------")
message("\n \t SQANTI3 rescue report")
message("\n--------------------------------------------------")

# Import pipe operator
require(magrittr)

####------------------------ INPUTS ---------------------------####




####------------------------ PLOT THEME ---------------------------####

# Load ggplot2
require(ggplot2)

# Install RColorConesa from GitHub if not already installed
pkg <- installed.packages() %>% rownames

if(!("RColorConesa" %in% pkg)){
  suppressMessages(install.packages("RColorConesa"))
}


# Set theme parameters (from SQANTI3_report.R)
sq_theme <- theme_classic(base_family = "Helvetica") +
  theme(plot.title = element_text(lineheight=.4, size=15, hjust = 0.5)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm")) +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=12) ) +
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size=12), 
        legend.key.size = unit(0.5, "cm"))

theme_set(sq_theme)


# Set category palette
cat.palette = c(FSM="#6BAED6", ISM="#FC8D59", NIC="#78C679", 
                NNC="#EE6A50", Genic="#969696", Antisense="#66C2A4", 
                Fusion="goldenrod1", Intergenic = "darksalmon", `Genic \nintron`="#41B6C4")
