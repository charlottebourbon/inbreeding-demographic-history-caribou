#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  FROH by ROH length class
#
#  Companion code for:
#    Bourbon et al. (2026) Inbreeding and Demographic History of Caribou
#    (Rangifer tarandus) in Western Canada Inferred From Genome-Wide SNP Data.
#    Evolutionary Applications. doi:10.1111/eva.70311
#
#  A modified version of detectRUNS::Froh_inbreedingClass() that accepts custom
#  ROH length categories. Partitions each individual's ROHs into length classes
#  and computes the FROH contributed by each class (class ROH length / genome
#  length covered by the map).
#
#  Length classes used in the paper (Mb):
#    short        0.3 - 2
#    short-medium 2   - 4
#    medium       4   - 6
#    medium-long  6   - 8
#    long         > 8
#
#  Requires: plyr, data.table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(plyr)
library(data.table)

# Total mapped length per chromosome, from a PLINK .map file.
# Returns a data frame with columns CHROMOSOME and CHR_LENGTH (bp).
chromosome_length <- function(map_file) {

  if (!file.exists(map_file)) stop(paste("file", map_file, "does not exist"))

  map <- data.table::fread(map_file, header = FALSE)
  colnames(map) <- c("CHR", "SNP_NAME", "CM", "POSITION")
  map <- map[map$POSITION != 0, ]   # drop unplaced markers (position 0)

  CHR <- POSITION <- NULL           # silence R CMD check notes

  length_genome <- ddply(map, .(CHR), summarize, CHR_LENGTH = max(POSITION))
  length_genome$CHR_LENGTH <- as.numeric(length_genome$CHR_LENGTH)

  message("Total genome length: ", sum(length_genome$CHR_LENGTH))
  length_genome
}

# FROH contributed by each ROH length class, per individual.
#   runs     : data frame of ROHs from detectRUNS::consecutiveRUNS.run()
#   map_file : path to the PLINK .map file used for ROH detection
# Returns one row per individual with a Froh_Class_* column per length class.
froh_inbreeding_class <- function(runs, map_file) {

  # Class breaks (Mb) and labels; the final open-ended class is ">8".
  range_mb   <- c(0.3, 2, 4, 6, 8, Inf)
  name_class <- c("0.3-2", "2-4", "4-6", "6-8", ">8")

  runs$MB    <- runs$lengthBps / 1e6
  runs$CLASS <- cut(runs$MB, breaks = range_mb, labels = name_class, right = FALSE)

  length_genome <- chromosome_length(map_file)
  genome_bp     <- sum(length_genome$CHR_LENGTH)

  message("Calculating FROH by length class")

  froh_class <- unique(runs[c("group", "id")])

  for (i in seq_along(name_class)) {
    class_name <- name_class[i]

    if (class_name == ">8") {
      subset_roh <- runs[runs$MB >= 8, ]
    } else {
      subset_roh <- runs[runs$MB >= range_mb[i] & runs$MB < range_mb[i + 1], ]
    }
    if (nrow(subset_roh) < 1) next

    froh_temp <- ddply(subset_roh, .(id), summarize, sum = sum(lengthBps))
    froh_temp[[paste0("Froh_Class_", class_name)]] <- froh_temp$sum / genome_bp
    colnames(froh_temp)[2] <- paste0("Sum_Class_", class_name)

    froh_class <- merge(froh_class, froh_temp, by = "id", all = TRUE)
  }

  froh_class[is.na(froh_class)] <- 0
  froh_class
}
