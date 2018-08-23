#nacteni dat

Nacteni <- function()
{
  nazev_souboru <- file.choose()
  mojedata <- read.table(nazev_souboru,
                      header = T, sep = "", dec = ".",
                      colClasses=c("Clock" ="character"))
  return(mojedata)

}
