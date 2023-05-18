    


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

