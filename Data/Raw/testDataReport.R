# run data reports for each species, 
# and save to the current directory

library(bookdown)


# First, try for dover sole
specName <- "dover"

specList <- c("Dover","English","Rock","Petrale","Arrowtooth")


renderDataReport <- function( specName )
{

  # Create parameter list for rendering the document
  params <- list( species= specName,
                  fYear= 1954,
                  lYear= 2018,
                  lowL= 15,
                  hiL= 45,
                  rootDir= ".")
  # Make an output file name
  outFile <- paste( "dataReport_", specName, ".html", sep = "")

  # Render
  rmarkdown::render(  input = "dataReportTemplate.Rmd", 
                      output_file = outFile,
                      params = params,
                      envir = new.env(),
                      output_format = "bookdown::html_document2")

  # remove temporary files
  dataReportFiles <- paste("dataReport_",specName,"_files",sep = "")
  unlink(dataReportFiles, recursive = TRUE)


}

lapply( X = specList, FUN = renderDataReport )