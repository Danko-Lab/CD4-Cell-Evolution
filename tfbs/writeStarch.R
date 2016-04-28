load("rtfbsdb-score.rdata")

write.starchbed <- function(bed, file.starch) {
        # pipe bed into starch file
        write.table(bed, file = pipe(paste("sort-bed - | starch - >", file.starch)), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
        if( !file.exists(file.starch) )
                cat("! Failed to write starch file (", file.starch, ") using the sort-bed and starch commands.\n");
}

str.u<- sapply(1:NROW(tf.u), function(x) { y <- tf.u[x, c(10:12)]; y <- y[!is.na(y)]; return(y[1]) })
str.pi<- sapply(1:NROW(tf.pi), function(x) { y <- tf.pi[x, c(10:12)]; y <- y[!is.na(y)]; return(y[1]) })

write.starchbed(data.frame(tf.u$V1, tf.u$V2, tf.u$V3, tf.u$motif.id, tf.u$max.score, str.u, tf.u$tf.name, tf.u$score1, tf.u$score2, tf.u$score3), "tf.u.hg19.bed.starch")
write.starchbed(data.frame(tf.pi$V1, tf.pi$V2, tf.pi$V3, tf.pi$motif.id, tf.pi$max.score, str.pi, tf.pi$tf.name, tf.pi$score1, tf.pi$score2, tf.pi$score3), "tf.pi.hg19.bed.starch")

