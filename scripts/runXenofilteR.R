#!/usr/bin/env Rscript
# 
# ./runXenoFilterR.R -s samples.txt -o ./ test -n test
# input sample file contains 2 columns (reference_bam, host_bam)
#
suppressPackageStartupMessages(library("XenofilteR"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
	make_option(c("-s", "--samples"), help="Sample list"),
	make_option(c("-o", "--out"), help="outputFile name"),
	make_option(c("-n", "--name"), help="outputFile name")
    )
opt <- parse_args(OptionParser(option_list=option_list))
sample.file <-opt$s
destination.folder <- opt$o
output.names <- opt$n

sample.list <- read.table(sample.file)
bp.param <- SnowParam(workers = 1, type = "SOCK")
XenofilteR(sample.list, destination.folder = destination.folder, bp.param = bp.param, output.names=output.names)