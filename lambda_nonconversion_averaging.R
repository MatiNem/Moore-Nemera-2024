##This code averages lambda values for background nonconversion rate subtraction
library(data.table)
library(dplyr)

PV_WT1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041642_summary.tsv")
PV_WT2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041644_summary.tsv")

PV_KO1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041641_summary.tsv")
PV_KO2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041643_summary.tsv")

L5_WT1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041645_summary.tsv")
L5_WT2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041648_summary.tsv")

L5_KO1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041646_summary.tsv")
L5_KO2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041647_summary.tsv")

L4_WT1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041634_summary.tsv")
L4_WT2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB042981_summary.tsv")

L4_KO1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041633_summary.tsv")
L4_KO2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB042982_summary.tsv")

SST_WT1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB042979_summary.tsv")
SST_WT2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB042980_summary.tsv")
SST_WT3 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB045574_22CGWFLT4_S263_summary.tsv")

SST_KO1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB042978_summary.tsv")
SST_KO2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB045575_22CGWFLT4_S264_summary.tsv")

#PV WT LIB041642
#PV WT LIB041644
#PV KO LIB041641
#PV KO LIB041643

#L5 WT LIB041645
#L5 WT LIB041647
#L5 KO LIB041646
#L5 KO LIB041648

#L4 WT LIB041634
#L4 WT LIB042981
#L4 KO LIB041633
#L4 KO LIB042982

#SST WT LIB042979
#SST WT LIB042980
#SST WT LIB045574
#SST KO LIB042978
#SST KO LIB045575

PV_WT_lambda <- mean(c(PV_WT1[V1=="Lambda", V2], PV_WT2[V1=="Lambda", V2]))
PV_KO_lambda <- mean(c(PV_KO1[V1=="Lambda", V2], PV_KO2[V1=="Lambda", V2]))
L5_WT_lambda <- mean(c(L5_WT1[V1=="Lambda", V2], L5_WT2[V1=="Lambda", V2]))
L5_KO_lambda <- mean(c(L5_KO1[V1=="Lambda", V2], L5_KO2[V1=="Lambda", V2]))
L4_WT_lambda <- mean(c(L4_WT1[V1=="Lambda", V2], L4_WT2[V1=="Lambda", V2]))
L4_KO_lambda <- mean(c(L4_KO1[V1=="Lambda", V2], L4_KO2[V1=="Lambda", V2]))
PV_WT_KO_lambda <- mean(c(PV_WT1[V1=="Lambda", V2], PV_WT2[V1=="Lambda", V2], PV_KO1[V1=="Lambda", V2], PV_KO2[V1=="Lambda", V2]))
L5_WT_KO_lambda <- mean(c(L5_WT1[V1=="Lambda", V2], L5_WT2[V1=="Lambda", V2], L5_KO1[V1=="Lambda", V2], L5_KO2[V1=="Lambda", V2]))
L4_WT_KO_lambda <- mean(c(L4_WT1[V1=="Lambda", V2], L4_WT2[V1=="Lambda", V2], L4_KO1[V1=="Lambda", V2], L4_KO2[V1=="Lambda", V2]))

SST_WT_lambda <- mean(c(SST_WT1[V1=="Lambda", V2], SST_WT2[V1=="Lambda", V2], SST_WT3[V1=="Lambda", V2]))
SST_KO_lambda <- mean(c(SST_KO1[V1=="Lambda", V2], SST_KO2[V1=="Lambda", V2]))
SST_WT_KO_lambda <- mean(c(SST_WT1[V1=="Lambda", V2], SST_WT2[V1=="Lambda", V2], SST_WT3[V1=="Lambda", V2], SST_KO1[V1=="Lambda", V2], SST_KO2[V1=="Lambda", V2]))

avg_nonconv <- rbind(
  cbind(label="PV_WT_lambda", nonconversion_rate = PV_WT_lambda),
  cbind(label="PV_KO_lambda", nonconversion_rate = PV_KO_lambda),
  cbind(label="L5_WT_lambda", nonconversion_rate = L5_WT_lambda),
  cbind(label="L5_KO_lambda", nonconversion_rate = L5_KO_lambda),
  cbind(label="L4_WT_lambda", nonconversion_rate = L4_WT_lambda),
  cbind(label="L4_KO_lambda", nonconversion_rate = L4_KO_lambda),
  cbind(label="SST_WT_lambda", nonconversion_rate = L4_WT_lambda),
  cbind(label="SST_KO_lambda", nonconversion_rate = L4_KO_lambda),
  cbind(label="PV_WT_KO_lambda", nonconversion_rate = PV_WT_KO_lambda),
  cbind(label="L5_WT_KO_lambda", nonconversion_rate = L5_WT_KO_lambda),
  cbind(label="L4_WT_KO_lambda", nonconversion_rate = L4_WT_KO_lambda),
  cbind(label="SST_WT_KO_lambda", nonconversion_rate = SST_WT_KO_lambda)
) %>% data.table

write.table(avg_nonconv, file="HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

label1 <- "PV_WT_lambda"

avg_nonconv[label==label1, as.numeric(nonconversion_rate)]

#table of all lambda nonconversion rates, not averaged


lambda_nonconv <- data.table(rbind(
    cbind(GTAC_ESP_ID="LIB041642", label="PV_WT_rep1", nonconversion_rate=PV_WT1[V1=="Lambda", V2]), #PV WT LIB041642
    cbind(GTAC_ESP_ID="LIB041644", label="PV_WT_rep2", nonconversion_rate=PV_WT2[V1=="Lambda", V2]), #PV WT LIB041644
    cbind(GTAC_ESP_ID="LIB041641", label="PV_KO_rep1", nonconversion_rate=PV_KO1[V1=="Lambda", V2]), #PV KO LIB041641
    cbind(GTAC_ESP_ID="LIB041643", label="PV_KO_rep2", nonconversion_rate=PV_KO2[V1=="Lambda", V2]), #PV KO LIB041643
    
    cbind(GTAC_ESP_ID="LIB041645", label="L5_WT_rep1", nonconversion_rate=L5_WT1[V1=="Lambda", V2]), #L5 WT LIB041645
    cbind(GTAC_ESP_ID="LIB041648", label="L5_WT_rep2", nonconversion_rate=L5_WT2[V1=="Lambda", V2]), #L5 WT LIB041648
    cbind(GTAC_ESP_ID="LIB041646", label="L5_KO_rep1", nonconversion_rate=L5_KO1[V1=="Lambda", V2]), #L5 KO LIB041646
    cbind(GTAC_ESP_ID="LIB041647", label="L5_KO_rep2", nonconversion_rate=L5_KO2[V1=="Lambda", V2]), #L5 KO LIB041647
    
    cbind(GTAC_ESP_ID="LIB041634", label="L4_WT_rep1", nonconversion_rate=L4_WT1[V1=="Lambda", V2]), #L4 WT LIB041634
    cbind(GTAC_ESP_ID="LIB042981", label="L4_WT_rep2", nonconversion_rate=L4_WT2[V1=="Lambda", V2]), #L4 WT LIB042981
    cbind(GTAC_ESP_ID="LIB041633", label="L4_KO_rep1", nonconversion_rate=L4_KO1[V1=="Lambda", V2]),  #L4 KO LIB041633
    cbind(GTAC_ESP_ID="LIB042982", label="L4_KO_rep2", nonconversion_rate=L4_KO2[V1=="Lambda", V2]),  #L4 KO LIB042982
    
    cbind(GTAC_ESP_ID="LIB042979", label="SST_WT_rep1", nonconversion_rate=SST_WT1[V1=="Lambda", V2]), #SST WT LIB042979
    cbind(GTAC_ESP_ID="LIB042980", label="SST_WT_rep2", nonconversion_rate=SST_WT2[V1=="Lambda", V2]), #SST WT LIB042980
    cbind(GTAC_ESP_ID="LIB045574", label="SST_WT_rep3", nonconversion_rate=SST_WT3[V1=="Lambda", V2]), #SST WT LIB045574
    cbind(GTAC_ESP_ID="LIB042978", label="SST_KO_rep1", nonconversion_rate=SST_KO1[V1=="Lambda", V2]), #SST KO LIB042978
    cbind(GTAC_ESP_ID="LIB045575", label="SST_KO_rep2", nonconversion_rate=SST_KO2[V1=="Lambda", V2]))) #SST KO LIB045575


write.csv(lambda_nonconv, file="HG_lab/Mati/GabelLab/bisulfite/summary_tsv/intact_lambda_nonconversion_table.csv", quote=FALSE, row.names=FALSE)

