library(data.table)
library(dplyr)

    


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


chromosome.lengths = get.chromosome.lengths()


compile.histogram.data = function(which.one) {
    directory.path = which.one
    name.parts = strsplit(which.one, split="\\.")[[1]]
    parameter = name.parts[1]
    hilo = name.parts[2]
    filename.suffix = paste0("\\.txt\\.window\\.", parameter, "\\.", hilo , "$")
    # print(paste0("filename.suffix = ", filename.suffix))
    results.all.df = data.frame(species=character(), percentages=numeric(), stringsAsFactors=FALSE)
    for (f in list.files(directory.path)) {
        # print(paste0("f = ", f))        
        species = gsub(filename.suffix, "", f)
        # print(paste0("sp = ", species))
        file.path = paste0(directory.path, "/", f)
        # print(paste0("file.path = ", file.path))
        file.data = read.table(file.path, header=TRUE, sep="\t")
        # print(file.data)
        max.loc = subset(chromosome.lengths, filename==species)$maxloc
        # print(paste0(species, "\t", max.loc))
        percentages = file.data$`low.high` / max.loc * 100
        temp.df = data.frame(species, percentages)
        results.all.df = bind_rows(results.all.df, temp.df)
    }
    return(results.all.df)
}



make.histogram = function(which.one) {
    data = compile.histogram.data(which.one)
    name.parts = strsplit(which.one, split="\\.")[[1]]
    parameter = name.parts[1]
    hilo = name.parts[2]
    if (parameter == "dnds") {
        cap = "dN/dS"
    } else if (parameter == "dn") {
        cap = "dN"
    } else {
        cap = "dS"
    }
    hotcold = ifelse(hilo=="high", "Hot", "Cold")
    color = ifelse(hilo=="high", "darkred", "skyblue2")
    title = paste0("Distribution of ", cap, " ", hotcold, "spot Locations")
    filename = paste0("all\\.", tolower(hotcold), "spots\\.hist\\.", parameter, "\\.jpg")
    path = paste("thesis_graphs/", filename)
    jpeg(path)
    hist(data$percentages, col=color, breaks=50, main=title, xlab="Fractional distance between Ori and Ter (%)", ylim=c(0, 150))
    dev.off()
}

