



dir_path = "sliding_window"

get10 = function(df, col, hilo) {
    if (hilo == "hi") {
        df.sorted = df[order(-df[[ paste0(col, "z")      ]]),]
    } else {
        df.sorted = df[order(df[[ paste0(col, "z")      ]]),]
    }
    df.short = head(df.sorted, 10)
    return(df.short)
}






plot.hilo.z = function(species, parameter) {
    path = paste0(dir_path, "/", species, ".txt.window")
    df = read.table(path, header=TRUE, sep='\t')
    df = df[complete.cases(df),]
    if (nrow(df) < 20) {
        return("not enough data")
    }

    numgenes.mean = mean(df$numgenes)
    numgenes.sd = sd(df$numgenes)
    numgenes.plusone = numgenes.mean + numgenes.sd
    numgenes.minusone = numgenes.mean - numgenes.sd
    numgenes.max = max(df$numgenes)
    numgenes.min = min(df$numgenes)
    
    for (col in c("dnds", "dn", "ds", "numgenes")) {
        df[[paste0(col, "z")]] = scale(df[[col]])
    }

    outliers = get10(df, parameter, "hi")
    coldspots = get10(df, parameter, "lo")

    maxh = max(outliers[[paste0(parameter, "z")]])
    minc = min(coldspots[[paste0(parameter, "z")]])
    maxh = ifelse( maxh > 0, 1.05 * maxh, 0.95 * maxh)
    minc = ifelse( minc > 0, .95 * minc, 1.05 * minc)

    jpeg(paste0("window_populations/", species, ".jpg"))
    plot(outliers$numgenes, outliers[[paste0(parameter, "z")]],
         col="darkred",
         pch = 19,
         ylab = paste0(parameter, " Z-scores"),
         xlab = "Gene Count",
         xlim = c(.95 * numgenes.min, 1.05 * numgenes.max),
         ylim = c(minc, maxh),
         main = paste(parameter, "vs.", "gene population:", "Hot and Coldspots for", species))
    abline(v=numgenes.mean, col="red", lty="solid")
    abline(v=numgenes.plusone, col="red", lty=2)
    abline(v=numgenes.minusone, col="red", lty=2)
    abline(v=numgenes.max, col="blue", lty="solid")
    abline(v=numgenes.min, col="blue", lty="solid")
    points(coldspots$numgenes, coldspots[[paste0(parameter, "z")]],
           pch=16, col="dodgerblue")
    dev.off()
    
    return(outliers)
}
