featureNames <- c(
  paste0("mi", 1:20), "N", paste0("nc", 1:9), 
  sapply(paste0(paste0("c", 1:9), "_"), function (x) paste0(x, c(
    "end3'distMean", "end3'distMedian", "end3'distVar", 
    "stopdistMean", "stopdistMedian", "stopdistVar", "bindEnMin", 
    "bindEnMean", "bindEnVar", "bindEnMedian", "5'UTRFreq", "CDSFreq", 
    "3'UTRFreq", "RISCAccEngMean", "RISCAccEngMedian", "RISCAccEngMin", 
    "RISCAccEngVar", "NearEngDiffMean", "NearEngDiffMin", 
    'NearEngDiffMedian', "NearEngDiffVar"
    ))),
  sapply(paste0(paste0("c", 6:9), "_"), function (x) paste0(x, c(
    "freqStart3", "freqStart4", "freqStart5"))),
  sapply(paste0(paste0("c", 1:9), "_"), function (x) paste0(x, c(
         paste0("freqA", c(-10:-1, 1:30)), 
         paste0("freqU", c(-10:-1, 1:30)), 
         paste0("freqG", c(-10:-1, 1:30)))
  ))
)

featureDescriptions <- c(
  paste("miRNA: Identity of nuleotide at position", 1:20),
  "Total number of sites found",
  paste("Number of sites belonging to class", c(
        "1: Perfect wc match to 2-7 plus either: adenine at position 1, match at position 8, match at positions 1 and 8",
        "2: Match to 2-7 with GU pairing plus either: adenine at position 1, match at position 8, match at positions 1 and 8",
        "3: WC match to 2-7 with one mismatch plus either: adenine at position 1, match at position 8, match at positions 1 and 8",
        "4: Match to 2-7 with GU pairing and one mismatch plus either: adenine at position 1, match at position 8, match at positions 1 and 8",
        "5: Perfect wc match to 2-7",
        "6: Perfect 11nt WC match starting from 3, 4 or 5",
        "7: 11nt match including GU pairing starting from 3, 4 or 5",
        "8: 11nt match except for one mismatch starting from 3, 4 or 5",
        "9: 11nt match including GU pairing except for one mismatch starting from 3, 4 or 5"
        )
    ),
  sapply(paste0("For sites of class ", 1:9, ":"), function(x) paste(x, c(
    "Mean distance from 3' end", "Median distance from 3' end",
    "Variance of distance from 3' end", "Mean distance from stop codon", 
    "Median distance from stop codon", 
    "Variance of distance from stop codon",
    "Minimum site binding energy value", "Mean site binding energy",
    "Variance of site binding energy", "Median site binding energy",
    "Frequency of sites in 5'UTR", "Frequency of sites in CDS",
    "Frequency of sites in 3'UTR", 
    "Mean energy needed to make sites accessible to RISC",
    "Median energy needed to make sites accessible to RISC",
    "Min value of energy needed to make a site accessible to RISC",
    "Variance of energy needed to make sites accessible to RISC",
    "Mean diffence in binding energy between sites and their nearby regions",
    "Median diffence in binding energy between sites and their nearby regions",
    "Min value of diffence in binding energy between a site and its nearby region",
    "Variance of diffence in binding energy between sites and their nearby regions"
  ))),
  sapply(paste0("For sites of class ", 6:9, ":"), function(x) paste(x, c(
    "Frequency of centered sites where the match starts at position 3",
    "Frequency of centered sites where the match starts at position 4",
    "Frequency of centered sites where the match starts at position 5"
  ))),
  sapply(paste0("For sites of class ", 1:9, ":"), function(x) paste(x, c(
      paste("Frequency of adenine at position", c(-10:-1, 1:30)), 
      paste("Frequency of uracil at position", c(-10:-1, 1:30)), 
      paste("Frequency of guanine at position", c(-10:-1, 1:30))
  )))
)

featureTable <- cbind(featureNames, featureDescriptions)
write.table(featureTable, "features.csv", sep = ",", row.names = FALSE)
