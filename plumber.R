# plumber.R

#* Return the result of pqsfinder() ... example: http://localhost:8000/pqs?dnaString=ATGGGGTGGGGAGGGGTGGGGATA
#* @param dnaString The dna string on which to look for quadruplexes
#* @get /pqs
function(dnaString){
  msg <- paste0("DNA sequence is: ", dnaString)
  library(pqsfinder)
  pqs <- pqsfinder(DNAString(dnaString), strand="+")
  msg2 <- paste0("the result is: ", as.character(pqs))
  list(msg, msg2)
}

