load("rtfbsdb-score.rdata")

write.starchbed <- function(bed, file.starch) {
        # pipe bed into starch file
        write.table(bed, file = pipe(paste("sort-bed - | starch - >", file.starch)), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
        if( !file.exists(file.starch) )
                cat("! Failed to write starch file (", file.starch, ") using the sort-bed and starch commands.\n");
}

fix.str <- function(tf) {
  str <- tf[,10]
  isna <- is.na(str);str[isna] <- tf[isna,11]
  isna <- is.na(str);str[isna] <- tf[isna,12]

  return(str)
}

str.u <- fix.str(tf.u)
str.pi<- fix.str(tf.pi)

#str.u<- sapply(1:NROW(tf.u), function(x) { y <- tf.u[x, c(10:12)]; y <- y[!is.na(y)]; return(y[1]) })
#str.pi<- sapply(1:NROW(tf.pi), function(x) { y <- tf.pi[x, c(10:12)]; y <- y[!is.na(y)]; return(y[1]) })

write.starchbed(data.frame(tf.u$chr, tf.u$start, tf.u$stop, tf.u$motif.id, tf.u$max.score, str.u, tf.u$tf.name, tf.u$score.human, tf.u$score.chimp, tf.u$score.rhesus), "tf.u.hg19.bed.starch")
write.starchbed(data.frame(tf.pi$chr, tf.pi$start, tf.pi$stop, tf.pi$motif.id, tf.pi$max.score, str.pi, tf.pi$tf.name, tf.pi$score.human, tf.pi$score.chimp, tf.pi$score.rhesus), "tf.pi.hg19.bed.starch")

