args <- commandArgs(trailingOnly = TRUE)

data = read.table(args[1])
d = data.frame(x=diff(diff(data$V2)), y=1:(nrow(data)-2))
cutoff=which.min(d$x^2+d$y^2)

write.table(cutoff, args[2], quote=FALSE, row.names=FALSE, col.names=FALSE)
