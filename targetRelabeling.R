mir23b_matrix <- mir23b_featMat[[3]]
mir23b_matrix <- mir23b_matrix[ , !colnames(mir23b_matrix) %in% c("Target")]
colnames(accToTarget_23b) <- c("GenBank.Accession", "Target")
mir23b_matrix <- merge(mir23b_matrix, accToTarget_23b, by = "GenBank.Accession")
mir23b_matrix <- mir23b_matrix[!(mir23b_matrix$Target == 1 & mir23b_matrix$N == 0), ]
featMatmir23b <- merge(mir23b_matrix, freeEngFeats_23b, 
                       by = "GenBank.Accession")
featMatmir23b <- convertMissingValues(featMatmir23b)

mir27a_matrix <- mir27a_featMat[[3]]
mir27a_matrix <- mir27a_matrix[ , !colnames(mir27a_matrix) %in% c("Target")]
colnames(accToTarget_27a) <- c("GenBank.Accession", "Target")
mir27a_matrix <- merge(mir27a_matrix, accToTarget_27a, by = "GenBank.Accession")
mir27a_matrix <- mir27a_matrix[!(mir27a_matrix$Target == 1 & mir27a_matrix$N == 0), ]
featMatmir27a <- merge(mir27a_matrix, freeEngFeats_27a, 
                       by = "GenBank.Accession")
featMatmir27a <- convertMissingValues(featMatmir27a)

mir3118_matrix <- mir3118_featMat[[3]]
mir3118_matrix <- mir3118_matrix[ , !colnames(mir3118_matrix) %in% c("Target")]
colnames(accToTarget_3118) <- c("GenBank.Accession", "Target")
mir3118_matrix <- merge(mir3118_matrix, accToTarget_3118, by = "GenBank.Accession")
mir3118_matrix <- mir3118_matrix[!(mir3118_matrix$Target == 1 & mir3118_matrix$N == 0), ]
featMatmir3118 <- merge(mir3118_matrix, freeEngFeats_3118, 
                       by = "GenBank.Accession")
featMatmir3118 <- convertMissingValues(featMatmir3118)

mir4307_matrix <- mir4307_featMat[[3]]
mir4307_matrix <- mir4307_matrix[ , !colnames(mir4307_matrix) %in% c("Target")]
colnames(accToTarget_4307) <- c("GenBank.Accession", "Target")
mir4307_matrix <- merge(mir4307_matrix, accToTarget_4307, by = "GenBank.Accession")
mir4307_matrix <- mir4307_matrix[!(mir4307_matrix$Target == 1 & mir4307_matrix$N == 0), ]
featMatmir4307 <- merge(mir4307_matrix, freeEngFeats_4307, 
                        by = "GenBank.Accession")
featMatmir4307 <- convertMissingValues(featMatmir4307)