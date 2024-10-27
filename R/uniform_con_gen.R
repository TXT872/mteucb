#' Uniform Confidence Bands for a Marginal Treatment Effect (MTE) Function .
#'
#' `uniform_con_gen` computes uniform confidence bands for a MTE function that is estimated by a semiparametric estimation method. For details, see Tsuda, Jin and Okui (2024+).
#'
#' @param outcome_name The name of the outcome
#' @param covariate_name The name of covariate
#' @param instrument_name The name of instrument variables
#' @param treatment_name The name of a treatment variable
#' @param data data.frame to be used for the estimation of MTE
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
#' @param bw The bandwidth used in estimating the nonparametric part of MTE
#' @param l_eval The minimum value of evaluation points for the estimated MTE
#' @param u_eval The maximum value of evaluation points for the estimated MTE
#' @param kernel The kernel function used in estimating the nonparametric part of MTE.
#' You can choose 'gaussK' (gaussian kernel), and 'EapK' (Epanechinikov kernel). The default choice is 'gaussK'.
#' @param covariate_value A covariate value to plot the estimate MTE.
#' @param significance_level The significance level for the construction of uniform confidence band.
#'
#' @return A list that contains the following elements:
#' \item{Data}{A data frame that contains the following elements: \cr outcome, covariate, instrument variables, treatment and propensity score}
#' \item{Estimation}{A data frame of the following results: \cr outcome, covariate, instrument variables, treatment and propensity score}
#' \item{Beta}{The estimated value of parametric parts of MTE}
#' \item{SE}{A standard error of each coeffecient of parametric parts of MTE}
#' \item{Plot}{A list that contains the ggplot elements for unifrom confidence band}
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom MASS
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
#' Supp_P_hut_U<-0.85
#' Supp_P_hut_L<-0.15
#' O_SS<-1000
#' B2<-seq(Supp_P_hut_L,Supp_P_hut_U,length=O_SS)
#' covariate_value<-c(0,0,0,0)
#' Final_Result_Ex<-uniform_con_gen(outcome_name="Y",
#'                                  covariate_name=c("X1","X2","X3","X4"),
#'                                  instrument_name=c("Z1","Z2","Z3","Z4"),
#'                                  treatment_name= "D",
#'                                  data=demo,
#'                                  family="probit",
#'                                  trim=0.01,
#'                                  intercept=TRUE,
#'                                  se_type="HC2",
#'                                  bw=0.1,
#'                                  l_eval=0.15,
#'                                  u_eval=0.85,
#'                                  covariate_value=covariate_value,
#'                                  significance_level=0.05)
#' }
#'
#'
#'
#'
#' @references
#' Tsuda, Jin and Okui (2024+).
#' Uniform Confidence Band for Marginal Treatment Effect Function.
#' will be available on arXiv
uniform_con_gen<-function(outcome_name=outcome_name,
                          covariate_name=covariate_name,
                          instrument_name=instrument_name,
                          treatment_name=treatment_name,
                          data=data,
                          family="probit",
                          trim=0.01,
                          intercept=TRUE,
                          se_type="HC2",
                          bw=bw,
                          l_eval=l_eval,
                          u_eval=u_eval,
                          kernel=locpol::gaussK,
                          covariate_value=covariate_value,
                          significance_level=0.05){
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
  colnames(Dataframe)<-c(outcome_name, covariate_name, instrument_name,treatment_name, "P_hut")

  covariate_List<-covariate_name
  data.covariate<-as.data.frame(Dataframe[,covariate_List])
  colnames(data.covariate) <-covariate_name


  P_hut_covariate<- data.frame(matrix(nrow=SS,ncol=length(covariate_name)))
  for (i in  1:length(covariate_name)){
    P_hut_covariate[,i]<-(data.covariate[[i]])*Dataframe$P_hut
  }
  colnames(P_hut_covariate) <- paste("P_",covariate_name,sep="")
  DY<-Dataframe[,treatment_List]*Dataframe[,outcome_List]
  D0Y<-(1-Dataframe[,treatment_List])*Dataframe[,outcome_List]
  Dataframe<-cbind(Dataframe,P_hut_covariate,DY,D0Y)
  rownames(Dataframe)<-seq(1,SS,1)
  colnames(Dataframe)<-c(outcome_name, covariate_name, instrument_name,treatment_name, "P_hut",  paste("P_",covariate_name,sep=""),"DY","D0Y")


  Beta<-beta_gen_inpackage(outcome_name=outcome_List,
                           covariate_name=covariate_List,
                           instrument_name=instrument_List,
                           treatment_name=treatment_List,
                           P_hut=Dataframe$P_hut,
                           data=Dataframe,
                           se_type=se_type,
                           SS=SS)

  Estimation_Result<-unif_gen_inpackage(outcome_name=outcome_List,
                           covariate_name=covariate_List,
                           instrument_name=instrument_List,
                           treatment_name=treatment_List,
                           P_hut_name="P_hut",
                           beta0=Beta$Beta[,1],
                           beta1=Beta$Beta[,2],
                           data=Dataframe,
                           bw=bw,
                           l_eval=l_eval,
                           u_eval=u_eval,
                           kernel=kernel,
                           covariate_value=covariate_value,
                           significance_level=significance_level)

  Dataframe<-Dataframe[,c(outcome_List,covariate_List,instrument_List,treatment_List,"P_hut")]

  Estimation<-Estimation_Result$Estimation
  Plot<-Estimation_Result$Plot
  return(list(Data=Dataframe, Estimation=Estimation, Beta=Beta$Beta, SE=Beta$SE, Plot=Plot))
}
