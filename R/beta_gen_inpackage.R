#' @importFrom dplyr %>%
beta_gen_inpackage<-function(outcome_name=outcome_name,
                             covariate_name=covariate_name,
                             instrument_name=instrument_name,
                             treatment_name=treatment_name,
                             P_hut=P_hut,
                             data=data,
                             se_type=se_type,
                             SS=SS){
  covariate_List<-covariate_name
  data.covariate<-as.data.frame(data[,covariate_List])
  colnames(data.covariate)<-covariate_name

  instrument_List<-instrument_name
  data.instrument<-as.data.frame(data[,instrument_List])
  colnames(data.instrument)<-instrument_name

  outcome_List<-outcome_name
  data.outcome<-as.data.frame(data[,outcome_List])
  colnames(data.outcome)<-outcome_name

  treatment_List<-treatment_name
  data.treatment<-as.data.frame(data[,treatment_List])
  colnames(data.treatment)<-treatment_name

  outcome<-data.outcome
  covariate<-data.covariate
  instrument<-data.instrument
  treatment<-data.treatment

  Covariate_NP<-matrix(nrow=SS,ncol=length(covariate_List))
  for (i in 1:length(covariate_List)){
    colnames(data)[i+1]<-"X_cov"
    npexpectedExpl_lm_our<-locpol::locpol(X_cov~P_hut, data= data,kernel = locpol::gaussK,xeval=data[, "P_hut"])
    npexpectedExpl_mat_our<-npexpectedExpl_lm_our$lpFit
    npexpectedExpl_our<-numeric(SS)
    for(ll in 1:SS){
      npexpectedExpl_our[ll]<-mean(npexpectedExpl_mat_our[npexpectedExpl_mat_our$P_hut==P_hut[ll],2])
    }
    Covariate_NP[,i]<-npexpectedExpl_our
    colnames(data)[i+1]<-covariate_List[i]
  }
  Covariate_NP<-as.data.frame(Covariate_NP)

  npexpectedOutcome0_lm_our<-locpol::locpol(D0Y~P_hut,data, kernel = locpol::gaussK,xeval=data[, "P_hut"])
  npexpectedOutcome0_mat_our<-npexpectedOutcome0_lm_our$lpFit
  npexpectedOutcome0_our<-numeric(SS)
  for(i in 1:SS){
    npexpectedOutcome0_our[i]<-mean(npexpectedOutcome0_mat_our[npexpectedOutcome0_mat_our$P_hut==P_hut[i],2])
  }

  npexpectedOutcome1_lm_our<-locpol::locpol(DY~P_hut,data, kernel = locpol::gaussK,xeval=data[, "P_hut"])
  npexpectedOutcome1_mat_our<-npexpectedOutcome1_lm_our$lpFit
  npexpectedOutcome1_our<-numeric(SS)
  for(i in 1:SS){
    npexpectedOutcome1_our[i]<-mean(npexpectedOutcome1_mat_our[npexpectedOutcome1_mat_our$P_hut==P_hut[i],2])
  }

  Semi_Cov_0<-matrix(nrow=SS,ncol=length(covariate_name))
  for (i in 1:length(covariate_name)){
    Semi_Cov_0[,i]<-(1-data$P_hut)*(data[,covariate_List]-Covariate_NP)[,i]
  }
  Semi_Cov_1<-matrix(nrow=SS,ncol=length(covariate_name))
  for (i in 1:length(covariate_name)){
    Semi_Cov_1[,i]<-data$P_hut*(data[,covariate_List]-Covariate_NP)[,i]
  }

  nppartialOutcome0_our<-data$D0Y-npexpectedOutcome0_our
  nppartialOutcome1_our<-data$DY-npexpectedOutcome1_our


  est_result0<-estimatr::lm_robust(nppartialOutcome0_our~Semi_Cov_0+0,se_type=se_type)
  est_result1<-estimatr::lm_robust(nppartialOutcome1_our~Semi_Cov_1+0,se_type=se_type)

  beta0_est<-as.matrix(est_result0$coef)
  colnames(beta0_est)<-"beta0_est"
  rownames(beta0_est)<-covariate_List

  beta1_est<-as.matrix(est_result1$coef)
  colnames(beta1_est)<-"beta1_est"
  rownames(beta1_est)<-covariate_List

  beta0_se<-as.matrix(est_result0$std)
  colnames(beta0_se)<-"beta0_se"
  rownames(beta0_se)<-covariate_List

  beta1_se<-as.matrix(est_result1$std)
  colnames(beta1_se)<-"beta1_se"
  rownames(beta1_se)<-covariate_List

  return(list(Beta=cbind(beta0_est, beta1_est), SE=cbind(beta0_se, beta1_se)))
}

