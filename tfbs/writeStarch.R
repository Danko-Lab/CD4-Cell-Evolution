load("rtfbsdb-score.rdata")

write.starchbed <- function(bed, file.starch) {
        # pipe bed into starch file
        write.table(bed, file = pipe(paste("sort-bed - | starch - >", file.starch)), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
        if( !file.exists(file.starch) )
                cat("! Failed to write starch file (", file.starch, ") using the sort-bed and starch commands.\n");
}

write.starchbed(data.frame(tf.u$chr, tf.u$start, tf.u$end, tf.u$motif.id, tf.u$max.score, tf.u$tf.name, tf.u$score.human, tf.u$score.chimp, tf.u$score.rhesus), "tf.u.hg19.bed.starch")
write.starchbed(data.frame(tf.pi$chr, tf.pi$start, tf.pi$end, tf.pi$motif.id, tf.pi$max.score, tf.pi$tf.name, tf.pi$score.human, tf.pi$score.chimp, tf.pi$score.rhesus), "tf.pi.hg19.bed.starch")

