library(rtfbsdb)

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP( db );

motif_id   <- c( "M6203_1.02", "M4505_1.02", "M6151_1.02", "M4490_1.02", "M6352_1.02", "M6494_1.02", "M6180_1.02")
tf.name    <- c( "ELF1",       "GABPA",      "ARNT",       "YY1",        "MYCN",       "STAT2",      "CREB1")

tfbs.drawLogo(tfs, file.pdf="tfbs.drawLogo1.pdf", motif_id= motif_id);


#tfbs.drawLogo(tfbs, file.pdf, index=NA, tf_id=NA, motif_id=NA, tf_name=NA, family_name=NA, tf_status=NA, groupby=NA) 
