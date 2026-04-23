args <- commandArgs(TRUE)

library(HMMcopy)
library(MASS)

rfile      <- args[1]
gfile      <- args[2]
mfile      <- args[3]
resolution <- args[4]   # "50kb" or "10mb"

pattern     <- ifelse(resolution == "50kb", ".50kb.wig", ".10mb.wig")
output_file <- paste0(gsub(pattern, "", rfile), pattern, ".Normalization.txt")

raw <- wigsToRangedData(rfile, gfile, mfile)

# correctReadcount uses LOESS to fit GC bias.  It fails when a BAM has very
# few reads (e.g. fetus/mom TLEN-split BAMs at 10 mb resolution).  In that
# case write a fallback table so downstream tools receive a well-formed file.
normal <- tryCatch(
    correctReadcount(raw),
    error = function(e) {
        message(sprintf(
            "[HMMcopy] WARNING: correctReadcount failed (%s). Writing fallback uncorrected table.",
            conditionMessage(e)
        ))
        NULL   # signal that we must use the raw fallback path
    }
)

if (is.null(normal)) {
    # Build a plain data.frame from the RangedData and add the expected columns
    # so downstream parsers see the same schema as a successful run.
    raw_df <- as.data.frame(raw)
    # RangedData → data.frame gives: space, start, end, reads, gc, map
    # Add dummy correction columns; copy = log2(reads+1) as a crude proxy.
    reads_vec         <- as.numeric(raw_df[["reads"]])
    raw_df$valid      <- !is.na(reads_vec) & reads_vec > 0
    raw_df$ideal      <- raw_df$valid
    raw_df$cor.gc     <- reads_vec
    raw_df$cor.map    <- reads_vec
    raw_df$corrected.copy <- log2(pmax(reads_vec, 0) + 1)
    raw_df$copy       <- raw_df$corrected.copy

    write.table(raw_df, output_file,
                quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
} else {
    na.filter <- normal[-which(is.na(normal$copy))]
    all       <- as.data.frame(na.filter)

    for (i in c(1:22, "X", "Y")) {
        assign(paste0("chr", i), subset(all, all[, 1] == paste0("chr", i)))
    }

    write.table(normal, output_file,
                quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}
