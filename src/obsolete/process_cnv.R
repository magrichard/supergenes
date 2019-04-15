if (!exists("cnv_lusc")) {
  print("Loading LUSC copy number variations (cnv) data...")
  cnv_lusc = readRDS("~/projects/tcga_studies/study_TCGA-LUSC_cnv.rds")
}



process_cnv <- function(data_cnv = cnv_lusc$data, treshold = 0.3){

processed_cnv<- epimedtools::monitored_apply(data_cnv,1,
                  function(gene){
                    sapply(gene,
                      function(patient){
                        if(patient <= abs(treshold) & patient >= -abs(treshold)){patient <- 1}
                          else{patient <- 0}
                                   
    })
  })
  return(t(processed_cnv))
}