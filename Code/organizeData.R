library(readxl)
library(tidyverse)
library(stringr)

# Read first sheet and rearrange columns
d = read_excel("data/raw/Data_MLubell1_Clean2.xlsx", na = "NODATA")
d = d[, c(5:6, 1:4, 7:ncol(d))]  # First two columns must be the edgelist
# Write the edgelist to csv
structure(d[, c(1, 2, 3, 4, 8)], names = c("cat1", "cat2", "id", "source", "prey_kg")) %>% 
  write.csv("data/derived/edgelist.csv", row.names = FALSE)

# Read second sheet and clean it up
dyads = read_excel("data/raw/Data_MLubell1_Clean2.xlsx", na = "x", sheet = "All possible combos")
dyads = dyads[-1, c(1, 2, 5, 6)]  # Remove the second row and keep only meaningful columns
dyads = dyads[!apply(dyads, 1, function(x) sum(is.na(x)) == length(x)), ]  # Remove rows at the bottom that are totally empty


# Get them in the same form and row-bind them together

## First sheet:
e1 = d[, c(1, 2, 7, 10)]  # Just the four columns of interest from the first sheet
e1 = e1[!duplicated(e1), ]   # Just the not-duplicated rows
names(e1) = c("cat1", "cat2", "relatedness", "propOverlap")

## Second sheet
# Create columns for each cat:
dyads = 
  strsplit(dyads$`Cat Pair`, "_") %>%  # Split the pairs into component cats
  do.call(rbind, .) %>%  # Each row is returned as an item in a list. Bind them together as rows.
  data.frame(stringsAsFactors = FALSE) %>%  # Make it a data.frame, keeping from being factors
  structure(names = c("cat1", "cat2")) %>%   # Give it column names
  cbind(dyads, .)  # Column bind it to the dyads data.frame
# Want in the order given and in the other order:
e2a = structure(dyads[, c(5, 6, 2, 3)],  # Cat1 first and proportion of its overlap
                names = c("cat1", "cat2", "relatedness", "propOverlap"))
e2b = structure(dyads[, c(6, 5, 2, 4)],  # Cat2 first and proportion of its overlap
                names = c("cat1", "cat2", "relatedness", "propOverlap"))
allDyads = rbind(e1, e2a, e2b)  # Row-bind them together

# Which cats appear in second sheet but not the first?
# unique(c(dyads$cat1, dyads$cat2))[!unique(c(dyads$cat1, dyads$cat2)) %in% get.vertex.attribute(n, "vertex.names")]
## Mark didn't mention M62. Leave him in and ask about it.

remove = c("F59", "F96", "F97")
allDyads = filter(allDyads, !cat1 %in% remove & !cat2 %in% remove)

# We know UncF and UncF2 are going to be missing info, so remove them to test for what's missing/inconsistent
unc = c("UncF", "UncF2")
filtered = filter(allDyads, !cat1 %in% unc & !cat2 %in% unc) %>%
  arrange(cat1, cat2)

# All the pairs:
cats = unique(c(filtered$cat1, filtered$cat2))
allPairs = 
  expand.grid(cat1 = cats, cat2 = cats) %>% 
  filter(cat1 != cat2) %>%
  arrange(cat1, cat2)

# Which pairs are missing?
# Don't appear in dataset:
missing = allPairs[!paste0(allPairs$cat1, allPairs$cat2) %in% paste0(filtered$cat1, filtered$cat2), ]


if(FALSE) { # Wrote these, sent to ME, he sent back responses
  # Two rows have missing entries. Add them
  miss = 
    filtered[!complete.cases(filtered), 1:2] %>%
    rbind(missing) %>%
    cbind(relatedness = NA, propOverlap = NA)
  
  # Which pairs have conflicting data
  conflict = 
    group_by(filtered, cat1, cat2) %>%
    summarise(numRelatednesses = n_distinct(round(relatedness, 3)),
              numOverlaps = n_distinct(round(propOverlap, 3))) %>%
    filter(numRelatednesses > 1 | numOverlaps > 1) %>%
    semi_join(allDyads, .) %>%
    mutate(propOverlap = NA) %>%
    filter(!duplicated(.)) 
  
  write.csv(miss, "data/raw/missingData.csv", row.names = FALSE, na = "")
  write.csv(conflict, "RangeOverlapConflicts.csv", row.names = FALSE, na = "")
  
} else {  # Read ME's responses
  conflict = read.csv("data/raw/conflictCorrections.csv")
  miss = read.csv("data/raw/missingData.csv")
  
  # Put original data together with ME's corrections
  
  d = rbind(allDyads, miss)  # Add what was missing
  
  # Remove duplicates (including duplicated relatedness and overlap)
  d = d[!duplicated(d), ]
  # Remove three rows where one cat is an uncF
  d = d[!grepl("UncF", d$cat2), ]
  
  # Remove any pairs in the original dataset that were in conflict, then add corrected conflicts
  d = anti_join(d, conflict[, 1:2]) %>%
    rbind(conflict)
  
  # Leaves four dyads where one entry was missing or there was a tiny difference. 
  d[duplicated(d[, 1:2]) | duplicated(d[, 1:2], fromLast = T), ]
  # Average with removing NAs get us to 132 = 12 * 11 directed dyads!
  d = group_by(d, cat1, cat2) %>%
    summarise(relatedness = unique(relatedness),
              propOverlap = mean(propOverlap, na.rm = TRUE))
  
  # Create some variables
  d = 
    mutate(d,
           sex1 = str_sub(cat1, 1, 1),
           sex2 = str_sub(cat2, 1, 1),
           sexPair = ifelse(sex1 == sex2, paste0(sex1, sex2), "MF"))
  write.csv(d, "data/derived/dyadInfo.csv", row.names = FALSE)
}

