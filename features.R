featureDescription <- list(
  c("N", "total number of all sites"),
  c("mi1", "ide")
  )

featureNames <- c(
  paste0("mi", 1:20), "N", paste("nc", 1:7))

featureDescriptions <- c(
  paste("miRNA: Identity of nuleotide at position", 1:20),
  "Total number of sites found",
  paste("Number of sites belonging to class", c(
        "1: Perfect wc match to 2-7 plus either: Adenosine at position 1, match at position 8, match at positions 1 and 8",
        "2: Match to 2-7 with GU pairing plus either: Adenosine at position 1, match at position 8, match at positions 1 and 8",
        "3: WC match to 2-7 with one mismatch plus either: Adenosine at position 1, match at position 8, match at positions 1 and 8",
        "4: Match to 2-7 with GU pairing and one mismatch plus either: Adenosine at position 1, match at position 8, match at positions 1 and 8",
        "5: Perfect wc match to 2-7",
        "6: Centered. Perfect WC match "
        )
    )
  
  )
