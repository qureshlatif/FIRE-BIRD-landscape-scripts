##################################################################################
# Generate reference table that lists the state, year, and full ID for each fire #
##################################################################################

## The breakpoint for a year is July 1.
## Thus, all fires from July 1, 1984, to July 1, 1985, create habitat for 1985.

setwd("T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/wildfire_data/tars")
library(dplyr)
library(stringr)

files <- list.files()

# Columns for reference table
tar.file <- files[which(str_detect(files, ".tar.gz"))]
State <- tar.file %>% str_sub(1, 2) %>% toupper()
ID <- tar.file %>% str_sub(1, 21)
Date <- ID %>% str_sub(21-7, 21)
Date <- str_c(str_sub(Date, 1, 4), "-",
              str_sub(Date, 5, 6), "-",
              str_sub(Date, 7, 8))
Year <- as.integer(str_sub(Date, 1, 4))
Month <- as.integer(str_sub(Date, 6, 7))
late.months <- which(Month >= 7)
Year[late.months] <- Year[late.months] + 1
  # Fires that occurred July - December contribute to habitat in the following year.

write.table(data.frame(State, Date, Year, Month, ID, stringsAsFactors = F),
            "T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/wildfire_data/tars.txt",
            sep=",", row.names = F)
