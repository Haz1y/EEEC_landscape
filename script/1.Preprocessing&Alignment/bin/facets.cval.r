#!/usr/bin/env Rscript

library(facets)
library(data.table)
library(readr)
options(datatable.fread.input.cmd.message=FALSE)
source(system.file("extRfns", "readSnpMatrixDT.R", package="facets"))

datafile = commandArgs(TRUE)[1] # "MESO_054_T.chr1.pileup.gz"
genome = commandArgs(TRUE)[2] # "hg19"
cur_params = as.numeric(c(commandArgs(TRUE)[3:6])) # c(1000,25,300,150)
ndepth_param = as.numeric(commandArgs(TRUE)[7]) # 20

load('arm_position.Rdata')

#we check if we have to compute multiples cvalues
if (!is.na(commandArgs(TRUE)[8]) && commandArgs(TRUE)[8]=="MCVAL") {
	m_cval = TRUE
} else {
	m_cval = FALSE
}

if (!is.na(commandArgs(TRUE)[9]) && commandArgs(TRUE)[9]=="PDF") {
	plot_pdf = TRUE
} else {
	plot_pdf = FALSE
}

# function to calculate the size of each integral copy
get_intergral_size = function(fit, method = 'em') {
  cncf = fit$cncf
  cncf$size = cncf$end - cncf$start
  if(method == 'cncf') {
    tt = table(cncf$tcn)
    }
  else {
    tt = table(cncf$tcn.em)
    }
  if(length(tt) == 1) {
    integral_size = sum(cncf$size, na.rm = T)
    names(integral_size) = names(tt)
  } else {
    tt2 = names(tt)
    integral_size=vector()
    for(y in 1:length(tt2)){
      if(method == 'cncf') {
        integral_size[y] = sum(cncf$size[which(cncf$tcn == tt2[y])])
      } else {
        integral_size[y] = sum(cncf$size[which(cncf$tcn.em == tt2[y])])
      }
    }
    names(integral_size) = tt2
  }
  names(integral_size) = as.numeric(names(integral_size))
  return(integral_size)
}

# main function to calculate FGA / GAIN / LOSS / LOH / WGD
get_FGA = function(fit, out = NULL, sampleID = NULL, method = 'em', include_loh = F) {
  cncf = fit$cncf
  cncf$size = cncf$end - cncf$start
  integral_size = get_intergral_size(fit, method)
  cut_off = as.numeric(names(integral_size)) > 2
  is_WGD = sum(integral_size[cut_off], na.rm = T)/(sum(integral_size, na.rm = T)) > .5
  major_cn = as.numeric(names(which.max(integral_size)))
  gain = sum(integral_size[as.numeric(names(integral_size)) > major_cn], na.rm = T)/sum(integral_size, na.rm = T)
  loss = sum(integral_size[as.numeric(names(integral_size)) < major_cn], na.rm = T)/sum(integral_size, na.rm = T)
  if(method == 'cncf') {
    LOH = sum(cncf$size[which(cncf$tcn == major_cn & cncf$lcn == 0)], na.rm = T)/sum(integral_size, na.rm = T)
  } else {
    LOH = sum(cncf$size[which(cncf$tcn.em == major_cn & cncf$lcn.em == 0)], na.rm = T)/sum(integral_size, na.rm = T)
  }
  if(include_loh) {
    FGA = gain + loss + LOH
  } else {
    FGA = gain + loss
  }
  if(!is.null(out)) {
    sample_name = as.character(out$IGV[1,1])
  } else {
    sample_name = sampleID
  }
  if(is.null(fit$emflags)) { fit$emflags = '' }
	if(is.null(fit$loglik)) { fit$loglik = '' }
	if(is.na(fit$purity)) { fit$purity = 'NA' }
  output = data.frame('SAMPLE' = sample_name, 'PURITY' = fit$purity, 'PLOIDY' = fit$ploidy, 'FGA' = FGA, 'GAIN' = gain, 'LOSS'= loss, 'LOH' = LOH, 'N_BREAKPOINTS' = nrow(cncf), 'IS_WGD' = is_WGD, 'INTERGAL_PLOIDY' = major_cn, 'EM_FLAGS' = fit$emflags,'dipLogR'=fit$dipLogR,'loglik'= fit$loglik)
  return(output)
}

# function to compute FGA by chromosomal arm # need arm_position.Rdata
get_ARM_FGA = function(fit, arm_position, out= NULL, sampleID = NULL, method = 'em', include_loh = F, calls_threshold = 0.9) {
  require(data.table)
  integral_size = get_intergral_size(fit, method)
  major_cn = as.numeric(names(which.max(integral_size)))

  if(!is.null(out)) {
    sample_name = as.character(out$IGV[1,1])
  } else {
    sample_name = sampleID
  }

  cncf = fit$cncf
  cncf_pos = data.table(chr=cncf$chrom, loc.start=cncf$start,
                        loc.end=cncf$end)
  setkey(cncf_pos, chr, loc.start, loc.end)
  arm_position = arm_position[ !arm_position$arm %in% c('13p', '14p', '15p', '21p', '22p','23p', '23q')] # remove acrocentric and sex chromosome

  setkey(arm_position, chr, start, end)
  fo_impact.idx <- foverlaps(arm_position, cncf_pos, nomatch=NA, which = T)
  fo_impact = foverlaps(arm_position, cncf_pos, nomatch=NA)
  fo_impact[, loc.start := ifelse(loc.start < start, start, loc.start)]
  fo_impact[, loc.end := ifelse(loc.end > end, end, loc.end)]
  fo_impact$size = fo_impact$loc.end - fo_impact$loc.start

  cncf = cncf[fo_impact.idx$yid,]
  cncf$arm = fo_impact$arm
  cncf$size = fo_impact$size

  n_breaks = cf_gain = gain = cf_loss = loss = vector()
  for(k in 1:length(arm_position$arm)) {
    tmp = cncf[which(cncf$arm == arm_position$arm[k]),]
    if( sum(is.na(tmp$size)) ==1 ) {
      n_breaks[k] = gain[k] = loss[k] = NA
    } else {
      n_breaks[k] = nrow(tmp)
      if(method == 'cncf') {
        gain[k] = sum(tmp$size[which(tmp$tcn > major_cn)], na.rm = T)/sum(tmp$size, na.rm = T)
        cf_gain_tmp =  if(length(tmp$cf[which(tmp$tcn > major_cn)]) == 0) { NA } else { tmp$cf[which(tmp$tcn > major_cn)]}
        cf_gain[k] =   if(length(cf_gain_tmp) == 1) { cf_gain_tmp} else { cf_gain_tmp[which.max(tmp$size[which(tmp$tcn > major_cn)])] }
        loss[k] = sum(tmp$size[which(tmp$tcn < major_cn)], na.rm = T)/sum(tmp$size, na.rm = T)
        cf_loss_tmp =  if(length(tmp$cf[which(tmp$tcn < major_cn)]) == 0) { NA } else { tmp$cf[which(tmp$tcn < major_cn)]}
        cf_loss[k] = if(length(cf_loss_tmp) == 1) { cf_loss_tmp } else { cf_loss_tmp[which.max(tmp$size[which(tmp$tcn < major_cn)])] }

      } else {
        gain[k] = sum(tmp$size[which(tmp$tcn.em > major_cn)], na.rm = T)/sum(tmp$size, na.rm = T)
        cf_gain_tmp =  if(length(tmp$cf.em[which(tmp$tcn.em > major_cn)]) == 0) { NA } else { tmp$cf.em[which(tmp$tcn.em > major_cn)]}
        cf_gain[k] = if(length(cf_gain_tmp) == 1) { cf_gain_tmp} else { cf_gain_tmp[which.max(tmp$size[which(tmp$tcn.em > major_cn)])] }
        loss[k] = sum(tmp$size[which(tmp$tcn.em < major_cn)], na.rm = T)/sum(tmp$size, na.rm = T)
        cf_loss_tmp =  if(length(tmp$cf.em[which(tmp$tcn.em < major_cn)]) == 0) { NA } else { tmp$cf.em[which(tmp$tcn.em < major_cn)]}
        cf_loss[k] = if(length(cf_loss_tmp) == 1) { cf_loss_tmp } else { cf_loss_tmp[which.max(tmp$size[which(tmp$tcn.em < major_cn)])] }

      }
    }
  }

  altered_gain = if(sum(gain[!is.na(gain)] > calls_threshold) == 0) {
    NA
  } else {
    paste0(arm_position$arm[which(gain[!is.na(gain)] > calls_threshold)], '_gain')
  }

  altered_loss = if(sum(loss[!is.na(loss)] > calls_threshold) == 0) {
    NA
  } else {
    paste0(arm_position$arm[which(loss[!is.na(loss)] > calls_threshold)], '_loss')
  }
  altered_arm = c(altered_gain, altered_loss)
  altered_arm = altered_arm[which(!is.na(altered_arm))]
  if(length(altered_arm) > 0) {
    output = data.frame('SAMPLE_ID' = sample_name, 'altered_arm' = altered_arm,
                        'altered_arm_cf' = c(as.numeric(cf_gain[which(gain > calls_threshold)]), as.numeric(cf_loss[which(loss > calls_threshold)])), 'purity' = as.numeric(fit$purity))
    output$altered_arm_ccf = output$altered_arm_cf/output$purity
    return(output)
  } else {
	output = data.frame('SAMPLE_ID' = sample_name, 'altered_arm' = 'NA',
                        'altered_arm_cf' = 'NA', 'purity' = as.numeric(fit$purity),'altered_arm_ccf' = 'NA')
	return(output)
	}

}

#we set a random set for reproducible analysis PreprocSamples
set.seed(1234)

sample_name = gsub("snppileup.csv.gz","",datafile)
sample_name_no_dot = gsub(".snppileup.csv.gz","",datafile)
rcmat = readSnpMatrixDT(datafile)


plot_facets = function (oo_facets, fit_facets , text_title, plot_name, pref = "cval300", pdf = F) {
  if (pdf) {
  	pdf(paste(plot_name,pref,"_CNV.pdf",sep=""), width = 12, height = 11)
  } else {
  	png(paste(plot_name,pref,"_CNV.png",sep=""),width=30,height=27.5,units="cm",res=300)
  }
  plotSample(x = oo_facets, emfit = fit_facets, sname = text_title)
  dev.off()
  pdf(paste(plot_name,pref,"_CNV_spider.pdf",sep=""),width = 6,height = 6)
  logRlogORspider(oo_facets$out, oo_facets$dipLogR)
  dev.off()
}

#def parameters snp_nbhd = 1000   cval_preproc = 35 cval_proc1 = 300 cval_proc2 = 150 min_read_count = 20
#snp_nbhd: 500-5000:500 window size, the smaller the most number of SNPs used
#cval_preproc :
xx = preProcSample(rcmat, gbuild = genome, ndepth = ndepth_param, snp.nbhd = cur_params[1], cval = cur_params[2])
fit1=procSample(xx, cval=cur_params[3], min.nhet=15, dipLogR=NULL) #fit fine
#You can specify cval that is large enough to avoid hyper-segmentation
oo_fine=procSample(xx, cval=cur_params[4], min.nhet=15, dipLogR=fit1$dipLogR)#large fit to avoid short segments
fit_fine = emcncf(oo_fine)
#we plot the results
text_title=paste(sample_name_no_dot,": Purity=",round(fit_fine$purity,3)*100,"%; Ploidy=",round(fit_fine$ploidy,2),sep="")

pref=paste("def_cval",cur_params[4],sep="")

plot_facets(oo_fine, fit_fine, text_title, sample_name,pref,plot_pdf)
#we write the output table

stat = get_FGA(fit_fine, sampleID= sample_name_no_dot)
write_tsv(stat, paste(sample_name,pref,"_stats.tsv",sep=""))

arm_FGA = get_ARM_FGA(fit_fine, arm_position, sampleID= sample_name_no_dot)
write_tsv(arm_FGA, file=paste(sample_name,pref,"_arm_fga.tsv",sep=""))

fit_fine$cncf['cnlr.median-dipLogR'] = fit_fine$cncf$cnlr.median - fit_fine$dipLogR
write.table(fit_fine$cncf, file=paste(sample_name,pref,"_CNV.txt",sep=""), quote = F, sep = "\t", row.names = F)

#we have to compute multiples cvalues 500,1000 and 1500
if(m_cval){

#You can specify cval that is large enough to avoid hyper-segmentation
oo_fine=procSample(xx, cval=500, min.nhet=15, dipLogR=fit1$dipLogR)#large fit to avoid short segments
fit_fine = emcncf(oo_fine)
#we plot the results
text_title=paste(sample_name_no_dot,": Purity=",round(fit_fine$purity,3)*100,"%; Ploidy=",round(fit_fine$ploidy,2),sep="")
plot_facets(oo_fine, fit_fine, text_title, sample_name,"cval500",plot_pdf)
#we write the ouput table
stat = get_FGA(fit_fine, sampleID= sample_name_no_dot)
write_tsv(stat, file=paste(sample_name,"cval500","_stats.tsv",sep=""))

arm_FGA = get_ARM_FGA(fit_fine, arm_position, sampleID= sample_name_no_dot)
write_tsv(arm_FGA, file=paste(sample_name,"cval500","_arm_fga.tsv",sep=""))

fit_fine$cncf['cnlr.median-dipLogR'] = fit_fine$cncf$cnlr.median - fit_fine$dipLogR
write.table(fit_fine$cncf, file=paste(sample_name,"cval500","_CNV.txt",sep=""), quote = F, sep = "\t", row.names = F)

#You can specify cval that is large enough to avoid hyper-segmentation
oo_fine=procSample(xx, cval=1000, min.nhet=15, dipLogR=fit1$dipLogR)#large fit to avoid short segments
fit_fine = emcncf(oo_fine)
#we plot the results
text_title=paste(sample_name_no_dot,": Purity=",round(fit_fine$purity,3)*100,"%; Ploidy=",round(fit_fine$ploidy,2),sep="")
plot_facets(oo_fine, fit_fine, text_title, sample_name,"cval1000",plot_pdf)
#we write the ouput table
stat = get_FGA(fit_fine, sampleID= sample_name_no_dot)
write_tsv(stat, file=paste(sample_name,"cval1000","_stats.tsv",sep=""))

arm_FGA = get_ARM_FGA(fit_fine, arm_position, sampleID= sample_name_no_dot)
write_tsv(arm_FGA, file=paste(sample_name,"cval1000","_arm_fga.tsv",sep=""))

fit_fine$cncf['cnlr.median-dipLogR'] = fit_fine$cncf$cnlr.median - fit_fine$dipLogR
write.table(fit_fine$cncf, file=paste(sample_name,"cval1000","_CNV.txt",sep=""), quote = F, sep = "\t", row.names = F)

#You can specify cval that is large enough to avoid hyper-segmentation
oo_fine=procSample(xx, cval=1500, min.nhet=15, dipLogR=fit1$dipLogR)#large fit to avoid short segments
fit_fine = emcncf(oo_fine)
#we plot the results
text_title=paste(sample_name_no_dot,": Purity=",round(fit_fine$purity,3)*100,"%; Ploidy=",round(fit_fine$ploidy,2),sep="")
plot_facets(oo_fine, fit_fine, text_title, sample_name,"cval1500",plot_pdf)
#we write the ouput table
stat = get_FGA(fit_fine, sampleID= sample_name_no_dot)
write_tsv(stat, file=paste(sample_name,"cval1500","_stats.tsv",sep=""))

arm_FGA = get_ARM_FGA(fit_fine, arm_position, sampleID= sample_name_no_dot)
write_tsv(arm_FGA, file=paste(sample_name,"cval1500","_arm_fga.tsv",sep=""))

fit_fine$cncf['cnlr.median-dipLogR'] = fit_fine$cncf$cnlr.median - fit_fine$dipLogR
write.table(fit_fine$cncf, file=paste(sample_name,"cval1500","_CNV.txt",sep=""), quote = F, sep = "\t", row.names = F)

}
