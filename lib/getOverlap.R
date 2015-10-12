## OverlapSepect in R using rphast
require(rphast)
getOverlap <- function(BED1, BED2) {
  BED1_feat  <- feat(seqname= BED1[,1], start= BED1[,2], end= BED1[,3])
  BED2_feat <- feat(seqname= BED2[,1], start= BED2[,2], end= BED2[,3])
  ol <- overlap.feat(x= BED2_feat, filter= BED1_feat)
  pos_indx <- match(paste(ol$seqname, ol$start, ol$end), paste(BED2_feat$seqname, BED2_feat$start, BED2_feat$end))
  return(pos_indx)
}


