#' Coeffeicients of Parametric Part of Marginal Treatment Effect (MTE) Function .
#'
#' `bet_gen` estimates a parametric part of a MTE function through the estimation method invented by Carneiro and Lee (2009). For details, see Carneiro and Lee (2009).
#'
#' @param outcome_name The name of the outcome
#' @param covariate_name The name of covariate
#' @param instrument_name The name of instrument variables
#' @param treatment_name The name of a treatment variable
#' @param data data.frame to be used for the estimation of MTE.
#' @param family
#' The distribution used for the calculation of the propensity score.
#' You can choose "probit" or "logit". The default choice is "probit".
#' @param trim,
#' To mitigate the effect of ill behavior of estimtaed propensity scores,
#' we set those less than the value of trim or larger than (1-the value of trim) as the value of the value of trim  or (1-the value of trim).
#' @param intercept 'TRUE' or 'FLASE' option to decide whether you include an intercept
#' for estimating the propensity score or not. The default choice is 'TRUE'.
#' @param se_type The sort of standard error sought.
#' The options are "HC0", "HC1" (or "stata", the equivalent) , "HC2" (default), "HC3", or "classical".
#' @return A list that contains the following elements:
#' \item{Beta}{The estimated value of parametric parts of MTE}
#' \item{SE}{A standard error of each coeffecient of parametric parts of MTE}
#'
#' @export
#' @importFrom dplyr %>%
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' library(MASS)
#' SS<-2000
#' X1<-rnorm(n=SS,mean=0,sd=1)
#' X2<-rnorm(n=SS,mean=0,sd=1)
#' X3<-rnorm(n=SS,mean=0,sd=1)
#' X4<-rnorm(n=SS,mean=0,sd=1)
#' Z1<-rnorm(n=SS,mean=0,sd=1)
#' Z2<-rnorm(n=SS,mean=0,sd=1)
#' Z3<-rnorm(n=SS,mean=0,sd=1)
#' Z4<-rnorm(n=SS,mean=0,sd=1)
#' E0<-c(1,0.5,0.3)
#' E1<-c(0.5,1,0.8)
#' E2<-c(0.3,0.8,1)
#' E<-cbind(E0,E1,E2)
#' E4<-mvrnorm(n=SS,mu=c(0,0,0),Sigma=E)
#' U0<-E4[,1]
#' U1<-E4[,2]
#' V<-E4[,3]
#' Expl<-cbind(X1,X2,X3,X4)
#' Instr<-cbind(X1,X2,X3,X4,Z1,Z2,Z3,Z4)
#' PInstr<-cbind(Z1,Z2,Z3,Z4)
#' beta0<-c(0.5,0.1,-0.1,-0.5)
#' beta1<-c(0.8,0.4,-0.4,-0.8)
#' ganma0<-c(0.4,0.4,0.4,0.4,0.3,0.3,0.3,0.3)
#' D<-numeric(SS)
#' for(i in 1:SS){
#'   if ((ganma0 %*% (Instr[i,]))>V[i]){
#'     D[i]=1}
#' }
#' Y<-numeric(SS)
#' for(i in 1:SS){
#'   Y[i]<-D[i]*((beta1%*%Expl[i,])+U1[i])+(1-D[i])*((beta0%*%Expl[i,])+U0[i])
#' }
#' demo=as.data.frame(cbind(Y,Expl,PInstr,D))
#'
#'
#' RESULT_Beta<-beta_gen(outcome_name="Y",
#'                       covariate_name=c("X1","X2","X3","X4"),
#'                       instrument_name=c("Z1","Z2","Z3","Z4"),
#'                       treatment_name= "D",
#'                       data=demo,
#'                       family="probit",
#'                       trim=0.01,
#'                       intercept=TRUE,
#'                       se_type="HC2"
#')
#'}
#'
#'
#'
#' @references
#' Carneiro and Lee (2009)
#' Estimating distributions of potential outcomes using local instrumental variables with an application to changes in college enrollment and wage inequality.
#' Journal of Econometrics, 2009, vol. 149, issue 2, 191-208

beta_gen<-function(outcome_name=outcome_name,
                   covariate_name=covariate_name,
                   instrument_name=instrument_name,
                   treatment_name=treatment_name,
                   data=data,
                   family="probit",
                   trim=0.01,
                   intercept=TRUE,
                   se_type="HC2"){

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

  PHUT_data<-as.data.frame(cbind(treatment,covariate,instrument))
  if(intercept==TRUE){
    formula=as.matrix(PHUT_data[,treatment_name])~as.matrix(PHUT_data[,covariate_name])+as.matrix(PHUT_data[,instrument_name])
  }
  if(intercept==FALSE){
    formula=as.matrix(PHUT_data[,treatment_name])~-1+as.matrix(PHUT_data[,covariate_name])+as.matrix(PHUT_data[,instrument_name])
  }

  if(family=="probit"){
    PHUT_glm<-stats::glm(formula, family=stats::binomial(link="probit"))
    P_hut<-stats::pnorm((stats::model.matrix(formula, data =PHUT_data))%*%PHUT_glm$coefficients)
  }
  if(family=="logit"){
    PHUT_glm<-stats::glm(formula, family=stats::binomial(link="logit"))
    P_hut<-stats::plogis(stats::model.matrix(formula, data =PHUT_data)%*%PHUT_glm$coefficients, location = 0, scale = 1)
  }
  for(i in 1:nrow(PHUT_data)){
    if(P_hut[i]<trim){
      P_hut[i]<-trim
    }
    if(P_hut[i]>(1-trim)){
      P_hut[i]<-(1-trim)
    }
  }

  use_data <-as.data.frame(cbind(outcome,
                                 covariate,
                                 instrument,
                                 treatment,P_hut))

  colnames(use_data)<-c(outcome_name, covariate_name, instrument_name,treatment_name, "P_hut")
  Dataframe<-use_data
  preY1frame<-Dataframe[Dataframe[,treatment_List]==1,]
  preY0frame<-Dataframe[Dataframe[,treatment_List]==0,]
  minP_0<-min(preY1frame$P_hut)
  maxP_0<-(1-trim)
  minP_1<-trim
  maxP_1<-max(preY0frame$P_hut)
  Y1frame<-preY1frame[preY1frame$P_hut>minP_1 & preY1frame$P_hut<maxP_1,]
  Y0frame<-preY0frame[preY0frame$P_hut>minP_0 & preY0frame$P_hut<maxP_0,]
  Y1<-as.matrix(Y1frame)
  Y0<-as.matrix(Y0frame)
  Data<-rbind(Y1,Y0)
  SS<-nrow(Data)
  Dataframe<-as.data.frame(Data)


  P_hut_covariate<- data.frame(matrix(nrow=SS,ncol=length(covariate_List)))
  P_hut_covariate<-Dataframe[,covariate_List]*Dataframe$P_hut




  DY<-Dataframe[,treatment_List]*Dataframe[,outcome_List]
  D0Y<-(1-Dataframe[,treatment_List])*Dataframe[,outcome_List]
  Dataframe<-cbind(Dataframe,P_hut_covariate,DY,D0Y)
  rownames(Dataframe)<-seq(1,SS,1)
  colnames(Dataframe)<-c(outcome_name, covariate_name, instrument_name,treatment_name, "P_hut",  paste("P_",covariate_name,sep=""),"DY","D0Y")

  SS<-nrow(Dataframe)
  Covariate_NP<-matrix(nrow=SS,ncol=length(covariate_List))
  for (i in 1:length(covariate_List)){
    colnames(Dataframe)[i+1]<-"X_cov"
    npexpectedExpl_lm_our<-locpol::locpol(X_cov~P_hut, data= Dataframe,kernel=locpol::gaussK,xeval=Dataframe[, "P_hut"])
    npexpectedExpl_mat_our<-npexpectedExpl_lm_our$lpFit
    npexpectedExpl_our<-numeric(SS)
    for(ll in 1:SS){
      npexpectedExpl_our[ll]<-mean(npexpectedExpl_mat_our[npexpectedExpl_mat_our$P_hut==P_hut[ll,1],2])
    }
    Covariate_NP[,i]<-npexpectedExpl_our
    colnames(Dataframe)[i+1]<-covariate_List[i]
  }
  Covariate_NP<-as.data.frame(Covariate_NP)

  npexpectedOutcome0_lm_our<-locpol::locpol(D0Y~P_hut,Dataframe, kernel=locpol::gaussK,xeval=Dataframe[, "P_hut"])
  npexpectedOutcome0_mat_our<-npexpectedOutcome0_lm_our$lpFit
  npexpectedOutcome0_our<-numeric(SS)
  for(i in 1:SS){
    npexpectedOutcome0_our[i]<-mean(npexpectedOutcome0_mat_our[npexpectedOutcome0_mat_our$P_hut==P_hut[i,1],2])
  }

  npexpectedOutcome1_lm_our<-locpol::locpol(DY~P_hut,Dataframe, kernel =locpol::gaussK,xeval=Dataframe[, "P_hut"])
  npexpectedOutcome1_mat_our<-npexpectedOutcome1_lm_our$lpFit
  npexpectedOutcome1_our<-numeric(SS)
  for(i in 1:SS){
    npexpectedOutcome1_our[i]<-mean(npexpectedOutcome1_mat_our[npexpectedOutcome1_mat_our$P_hut==P_hut[i,1],2])
  }

  Semi_Cov_0<-matrix(nrow=SS,ncol=length(covariate_name))
  for (i in 1:length(covariate_name)){
    Semi_Cov_0[,i]<-(1-Dataframe$P_hut)*(Dataframe[,covariate_List]-Covariate_NP)[,i]
  }
  Semi_Cov_1<-matrix(nrow=SS,ncol=length(covariate_name))
  for (i in 1:length(covariate_name)){
    Semi_Cov_1[,i]<-Dataframe$P_hut*(Dataframe[,covariate_List]-Covariate_NP)[,i]
  }

  nppartialOutcome0_our<-Dataframe$D0Y-npexpectedOutcome0_our
  nppartialOutcome1_our<-Dataframe$DY-npexpectedOutcome1_our


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


