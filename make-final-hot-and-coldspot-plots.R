    


species.data.dir = "pi_merge_output"


get.chromosome.lengths = function() {
    chromosome.lengths = data.frame(filename=character(), max_loc=numeric(), stringsAsFactors=FALSE)
    file.list = list.files(species.data.dir)
    for (f in file.list) {
        file.path = paste0(species.data.dir, "/", f)
        file.data = read.table(file.path, header=TRUE, sep="\t")
        maxloc = max(file.data$loc)
#        print(paste0("f=", f))
        filename = gsub("\\.txt$", "", f)
#        print(paste0("filename=", filename))
        chromosome.lengths = rbind(chromosome.lengths, data.frame(filename, maxloc))
    }
    return(chromosome.lengths)
}




plotspots = function(species) {

    highdir = "zscores.reduced.high/"
    lowdir = "zscores.reduced.low/"
    high_suffix = ".txt.window.dnds.high"
    low_suffix = ".txt.window.dnds.low"
    ec_high_name = paste0(highdir, species, high_suffix)
    ec_low_name = paste0(lowdir, species, low_suffix)

    echigh = read.table(ec_high_name, header=T, sep='\t')
    eclow = read.table(ec_low_name, header=T, sep='\t')
    eclim = c(0, subset(results_df, filename=="Escherichia_coli")$max_loc)




    par(mfrow=c(2,1), oma=c(0, 0, 0, 0), mar=c(5, 4, 4, 2))
    title = paste0(species, ": Hotspot Distribution")
    plot(echigh$`low.high`, rep(1,10), pch=21, bg="indianred2", col="indianred2", ylim=c(0,1.5), xlab="Distance from origin of replication (bases)", ylab="", main=title, xlim=eclim, yaxt="n")
    axis(2, at=NULL, labels=FALSE, tick=FALSE)
    points(eclim, c(1.3, 1.3), pch=23, bg='blue')
    segments(eclim[1], 1.3, eclim[1], 0, lwd=2, col="blue")
    segments(eclim[2], 1.3, eclim[2], 0, lwd=2, col="blue")
    for (i in 1:10) {
        segments(echigh$`low.high`[i], 1, echigh$`low.high`[i], 0, lwd=2, col="indianred2")}
    abline(h=0, col="gray")
    text(x=0, y=1.3, labels="Ori", pos=3)
    text(x=eclim[2], y=1.3, labels="Ter", pos=3)

    
    title = paste0(species, ": Coldspot Distribution")
    plot(eclow$`low.high`, rep(1,10), pch=21, bg="skyblue2", col="skyblue2", ylim=c(0,1.5), xlab="Distance from origin of replication (bases)", ylab="", main=title, xlim=eclim, yaxt="n")
    axis(2, at=NULL, labels=FALSE, tick=FALSE)
    points(eclim, c(1.3, 1.3), pch=23, bg='red')
    segments(eclim[1], 1.3, eclim[1], 0, lwd=2, col="red")
    segments(eclim[2], 1.3, eclim[2], 0, lwd=2, col="red")
    for (i in 1:10) {
        segments(eclow$`low.high`[i], 1, eclow$`low.high`[i], 0, lwd=2, col="skyblue2")}
    abline(h=0, col="gray")
    text(x=0, y=1.3, labels="Ori", pos=3)
    text(x=eclim[2], y=1.3, labels="Ter", pos=3)
}
