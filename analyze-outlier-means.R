# initialize dataframes
alldata <- data.frame()
dshot <- data.frame()
dscold <- data.frame()
dnhot <- data.frame()
dncold <- data.frame()
dndshot <- data.frame()
dndscold <- data.frame()


dropn = function(df, n) {
    if (n > 0) {
        df[order(df$numgenes, decreasing=FALSE)[-seq(n)],]
    } else {
        df
    }
}

min.rows = 100
take.n = 20
drops = 0
min.window.population = 30


for (filename in list.files("sliding_window",  full.names = TRUE)) {
    
    print(paste0("reading: ", filename))
    data <- read.delim(filename, header = TRUE)
    data <- data[complete.cases(data), ]
    data <- data[data$numgenes > min.window.population,]
    if (nrow(data) < min.rows) next
    
    

    data[, c("dnz", "dndsz", "dsz")] <- scale(data[, c("dn", "dnds", "ds")])
    

    alldata <- rbind(alldata, data)
    
    newdshot = head(data[order(-data$dsz), ], n = take.n)
    newdshot = dropn(newdshot, drops)
    dshot <- rbind(dshot, newdshot)

    
    newdscold = head(data[order(data$dsz), ], n = take.n)
    newdscold = dropn(newdscold, drops)
    dscold <- rbind(dscold, newdscold)
    
    
    newdnhot = head(data[order(-data$dndsz), ], n = take.n)
    newdnhot = dropn(newdnhot, drops)
    dnhot <- rbind(dnhot, newdnhot)
    
    newdncold = head(data[order(data$dndsz), ], n = take.n)
    newdncold = dropn(newdncold, drops)
    dncold <- rbind(dncold, newdncold)

    
    newdndshot = head(data[order(-data$dnds), ], n = take.n)
    newdndshot = dropn(newdndshot, drops)
    dndshot <- rbind(dndshot, newdndshot)
    
    newdndscold = head(data[order(data$dnds), ], n = take.n)
    newdndscold = dropn(newdndscold, drops)
    dndscold <- rbind(dndscold, newdndscold)
    
}

do.test = function(df) {
    return(wilcox.test(df$numgenes, alldata$numgenes))
}


make.hist = function(df, title, color) {
    jpeg(paste0("window.pop.graphs/", title, ".jpg"))
    hist(df$numgenes, col=color, xlab="#genes in window", main=title)
    dev.off()
}


all.hists = function() {
    make.hist(alldata, "All.Data", "dodgerblue")
    make.hist(dshot, "dS.Hotspots", "darkred")
    make.hist(dscold, "dS.Coldspots", "dodgerblue")
    make.hist(dnhot, "dN.Hotspots", "darkred")
    make.hist(dncold, "dN.Coldspots", "dodgerblue")
    make.hist(dndshot, "dNdS.Hotspots", "darkred")
    make.hist(dndscold, "dNdS.Coldspots", "dodgerblue")
}
