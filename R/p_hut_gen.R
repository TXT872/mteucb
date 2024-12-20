#' Propensity Score Generation
#'
#' `p_hut_gen` generates the estimated propensity scores
#'
#' @param outcome_name The name of the outcome
#' @param covariate_name The name of covariate
#' @param instrument_name The name of instrument variables
#' @param treatment_name The name of a treatment variable
#' @param data data.frame to be used for the estimation of MTE
#' @param family The distribution used for the calculation of the propensity score.
#' You can choose "probit" or "logit". The default choice is "probit".
#' @param trim, To mitigate the effect of ill behavior of estimated propensity scores,
#' we set values below the trim level to the trim level and values above (1 - trim) level to (1 - trim) level, respectively. The default choice is “0.01”.
#' @param intercept 'TRUE' or 'FLASE' option to decide whether you include an intercept
#' for estimating the propensity score or not. The default choice is 'TRUE'.
#'
#'
#'
#' @return A list that contains the following elements:
#' \item{Data}{A data frame that contains the following elements: \cr outcome, covariate, instrument variables, treatment and propensity score}
#' \item{Sample Size}{Sample Size after the trimming}
#' \item{P_hut}{The estimated propensity scores}
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom MASS mvrnorm
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
#' RESULT_P_hut<-p_hut_gen(outcome_name="Y",
#'                         covariate_name=c("X1","X2","X3","X4"),
#'                         instrument_name=c("Z1","Z2","Z3","Z4"),
#'                         treatment_name= "D",
#'                         data=demo,
#'                         family="probit",
#'                         trim=0.01,
#'                         intercept=TRUE
#' )
#' }
#'
#'
#'
#'
p_hut_gen<-function(outcome_name=outcome_name,
                   covariate_name=covariate_name,
                   instrument_name=instrument_name,
                   treatment_name=treatment_name,
                   data=data,
                   family="probit",
                   trim=0.01,
                   intercept=TRUE){
  if(trim>1 || trim<0){
    stop(paste("Error: The value of trimming is incorrect"))
  }
  covariate_List<-covariate_name
  data.covariate<-as.data.frame(data[,covariate_List])
  colnames(data.covariate) <-covariate_name

  instrument_List<-instrument_name
  data.instrument<-as.data.frame(data[,instrument_List])
  colnames(data.instrument) <-instrument_name

  outcome_List<-outcome_name
  data.outcome<-as.data.frame(data[,outcome_List])
  colnames(data.outcome) <-outcome_name

  treatment_List<-treatment_name
  data.treatment<-as.data.frame(data[,treatment_List])
  colnames(data.treatment) <-treatment_name


  outcome<-data.outcome
  covariate<-data.covariate
  instrument<-data.instrument
  treatment<-data.treatment

  use_data <-as.data.frame(cbind(outcome,
                                 covariate,
                                 instrument,
                                 treatment))

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
  for(i in 1:nrow(use_data)){
    if(P_hut[i]<trim){
      P_hut[i]<-trim
    }
    if(P_hut[i]>(1-trim)){
      P_hut[i]<-(1-trim)
    }
  }
  Dataframe<-data.frame(cbind(use_data, P_hut))
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
  rownames(Dataframe)<-seq(1,SS,1)
  colnames(Dataframe)<-c(outcome_name, covariate_name, instrument_name,treatment_name, "P_hut")

  return(list(Data=Dataframe, SampleSize=SS, P_hut=Dataframe[,"P_hut"]))
}
