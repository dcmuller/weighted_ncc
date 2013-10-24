sink("./analysis/output/l05_getdata_epic.rlog", split = TRUE)
#! p05_getdata_epic.r
#! dcmuller
#! Create EPIC dataset for application of weighted analysis
##
#####################################################
library(foreign)
library(data.table)
library(rms)
sessionInfo()

#####################################################
epic <- read.dta("./data/source/epic_ls_dmuller_sept2013.dta")
epic <- data.table(epic)
setkey(epic, country, center, idepic, sex, d_birth, d_recrui)

## drop some variables
todrop <- c("mar_stat", "desc_per", "stillmen", "try_smok", "occ_smok", 
            "i_occ_smok", "shs_pchild", "shs_tchild", "shs_hom_wrk", 
            "shs_home", "shs_home_hrs", "shs_work", "shs_work_hrs", 
            "job_empl", "kidney_st", "cigars", "pipe", "i_cigret", 
            "a_diabet", "n_per_12m", "n_livbor_c", "shs_home_pac", 
            "n_occ_smok", "t_diabet", "s_kidney_st", "a_cigars", "a_givcs", 
            "n_cigars", "t_cigars", "i_cigars", "a_pipe", "a_givpi", "n_pipe", 
            "i_pipe", "n_pipe_t", "occ_stat", "pa_occ_stat", "self_empl", 
            "empl_people", "n_sbirth_c", "hip_adj", "waist_adj", "whr_adj", 
            "pa_score", "pa_score_m", "pa_work", "pa_vig", "m_diy", 
            "m_houswrk", "m_floors", "m_vigpa", "m_walk", "m_cycl", "m_gard", 
            "m_sport", "pa_mets", "pa_recr", "pa_hhld", "pa_index", 
            "pa_score_c", "pa_score_m_c", "pa_mets_c", "pa_mets_c_sx", 
            "pa_total", "pa_total_sx", "alc_beer20", "alc_beer30", 
            "alc_beer40", "alc_beer50", "alc_beerre", "alc_wine20", 
            "alc_wine30", "alc_wine40", "alc_wine50", "alc_winere", 
            "alc_spir20", "alc_spir30", "alc_spir40", "alc_spir50", 
            "alc_spirre", "alc_fwin20", "alc_fwin30", "alc_fwin40", 
            "alc_fwin50", "alc_fwinre", "alc_20", "alc_30", "alc_40", 
            "alc_50", "alc_30", "alc_40", "alc_lifetime_beer", 
            "alc_age_start_beer", "alc_age_stop_beer", "alc_drinktime_beer", 
            "alc_lifetime_wine", "alc_age_start_wine", "alc_age_stop_wine", 
            "alc_drinktime_wine", "alc_lifetime_spir", "alc_age_start_spir", 
            "alc_age_stop_spir", "alc_drinktime_spir", "alc_lifetime_fwin", 
            "alc_age_start_fwin", "alc_age_stop_fwin", "alc_drinktime_fwin", 
            "alc_drinker", "alc_re_c", "alc_lifetime_c", "alc_pattern")
todrop <- unique(todrop)
epic[, todrop := NULL, with = FALSE]
rm(todrop)


epic_fup <- read.dta("./data/source/eligible_ctrl_2013_jul.dta")
names(epic_fup) <- tolower(names(epic_fup))
epic_fup <- data.table(epic_fup)
setkey(epic_fup, country, center, idepic, sex, d_birth, d_recrui)

epic <- merge(epic, epic_fup)
rm(epic_fup)


###################################################
## Generate derived variables

## calculated age at baseline
epic$baseline_age <- as.numeric(epic$d_recrui - epic$d_birth)/365.25

## Integer age at recruitment
age_int <- floor(epic$age_recr)

## Country
epic$country_cde <- ifelse(epic$country == "B", 10, as.numeric(epic$country))
with(epic, table(country_cde, country, useNA="ifany"))
countrylabs <- c("France", "Italy", "Spain", "United Kingdom",     
                 "The Netherlands", "Greece", "Germany", "Sweden",
                 "Denmark", "Norway")
epic$country_cde <- factor(epic$country_cde, labels=countrylabs)
rm(countrylabs)

## lung cancers
lungcodes <- c("C340", "C341", "C342", "C343", "C348", "C349")
epic$lung_ca <- 0L
epic$lung_ca[epic$siteanyc %in% lungcodes] <- 1L
table(epic$behaanyc[epic$lung_ca==1], useNA = "ifany")
epic$lung_ca <- ifelse(epic$behaanyc==3 & !is.na(epic$behaanyc), 
                       epic$lung_ca, 0L)
rm(lungcodes)

## look at where to set end of follow up by country
diag_q <- quarters(epic$d_dganyc[epic$lung_ca==1])
diag_y <- year(epic$d_dganyc[epic$lung_ca==1])
table(diag_y, diag_q, epic$country[epic$lung_ca==1])
rm(diag_q, diag_y)

## set up end of follow-up
epic$endfup <- as.Date("2009-01-01")
epic$endfup[epic$country==1] <- as.Date("2005-07-01")
epic$endfup[epic$country==3] <- as.Date("2008-01-01")
with(epic, table(country, endfup))

epic$endfup <- with(epic, 
                    pmin(d_dganyc, endfup, d_endfup_frst, d_dthlst, na.rm=TRUE))
epic$endfup_age <- as.numeric(epic$endfup - epic$d_birth)/365.25


## time from baseline to diagnosis and end of follow up
epic$ttodiag_y <- as.numeric(epic$d_dganyc - epic$d_recrui) / 365.25
epic$ttoendfup_y <- as.numeric(epic$endfup - epic$d_recrui) / 365.25

## flag incident cases up to different follow up times
epic$case_all <- ifelse(epic$lung_ca == 1 & epic$endfup == epic$d_dganyc, 
                        1L, 
                        0L)

## generate a smoking status factor
epic$smoke <- ifelse(epic$smoke_stat==4, NA, epic$smoke_stat)
with(epic, table(smoke, smoke_stat, useNA="ifany"))
epic$smoke <- factor(epic$smoke, labels=c("never", "former", "current"))

## save the dataset dataset
save(epic, file = "./data/derived/d01_epic.RData")
sink()