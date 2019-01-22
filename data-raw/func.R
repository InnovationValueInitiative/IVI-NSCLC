cpi_adj <- function(year){
  cpi <- data.table(read_excel("medical-cpi.xlsx", skip = 11))
  row <- match(year, cpi$Year)
  adj <- cpi[Year == 2018, Jan]/cpi[row, Jan]
  return(adj)
}