alldata = data.frame()
dshot = data.frame()
dscold = data.frame()
dnhot = data.frame()
dncold = data.frame()
dndshot = data.frame()
dndscold = data.frame()

min.rows = 100
take.n = 20
drops = 0
min.window.population = 20

dropn = function(df, n) {
    if (n > 0) df[order(df$numgenes, decreasing=FALSE)[-seq(n)],]
    else df
}
                  
do.file = function(filename) {
    print(paste("reading:", filename))
    data = read.delim(filename, header=TRUE)
    data = data[complete.cases(data)]
    data = data[data$numgenes > min.window.population, ]
    if (nrow(data) < min.rows) return()

    data[, c("dnz", "dsz", "dndsz")] = scale(data[, c("dn", "ds", "dnds")])
    alldata = rbind(alldata, data)
        
