# plumber.R

#* @filter cors
cors <- function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$setHeader("Access-Control-Allow-Methods","*")
    res$setHeader("Access-Control-Allow-Headers", req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)
    res$status <- 200 
    return(list())
  } else {
    plumber::forward()
  }
}

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

#* Get options and dna strings, return job ID ... example: 
#* @param opts Options to configure algorithm
#* @param sequences The data on which to look for quadruplexes
#* @post /analyze
function(opts, sequences) {
  maxLength <- opts['maxLength']
  minScore <- opts['minScore']
  strandSense <- opts['strandSense']
  strandAnti <- opts['strandAnti']
  minLL <- opts['minLL']
  maxLL <- opts['maxLL']
  maxNB <- opts['maxNB']
  maxNM <- opts['maxNM']
  maxND <- opts['maxND']
  # TODO this should be resolved on FE
  if (strandSense == TRUE && strandAnti == TRUE) { 
    strand <- "*" 
  }
  else if (strandSense == TRUE) {
    strand <- "+"
  }
  else {
    strand <- "-" 
  }
  library(pqsfinder)
  for(seqDnaString in sequences['dnaString'][,1]) {
    pqs <- pqsfinder(DNAString(seqDnaString), strand=strand, max_len = maxLength, min_score = minScore, loop_min_len = minLL, loop_max_len = maxLL,
                     max_bulges = maxNB, max_mismatches = maxNM, max_defects = maxND);
    print(pqs)
  }
  for(seqDesc in sequences['seqDescription'][,1]) {
    print(seqDesc);
  }
  cat('\n')
  return(id)
}

