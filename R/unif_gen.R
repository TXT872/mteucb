#' Uniform Confidence Bands for a Marginal Treatment Effect (MTE) Function .
#'
#' `unif_gen` computes unifrom confidence bands of the MTE function given the parametric part of MTE is known. For details, see Tsuda, Jin and Okui (2024+).
#'
#' @param outcome_name The name of the outcome
#' @param covariate_name The name of covariate
#' @param instrument_name The name of instrument variables
#' @param treatment_name The name of a treatment variable
#' @param P_hut_name The name of estimated propensity scores.
#' @param beta0 A vector of parameters of a linear model for the outcome given an individual is not treated
#' @param beta1 A vector of parameters of a linear model for the outcome given an individual is treated
#' @param data data.frame including propensity scores to be used for the estimation of MTE
#' @param bw The bandwidth used in estimating the nonparametric part of MTE
#' @param l_eval The minimum value of evaluation points for the estimated MTE
#' @param u_eval The maximum value of evaluation points for the estimated MTE
#' @param kernel The kernel function used in estimating the nonparametric part of MTE.
#' You can choose 'gaussK' (gaussian kernel), and 'EapK' (Epanechinikov kernel). The default choice is 'gaussK'.
#' @param covariate_value A covariate value to plot the estimate MTE.
#' @param significance_level The significance level for the construction of uniform confidence band.
#'
#'
#' @return A list that contains the following elements:
#' \item{Estimation}{A data frame of the following results: \cr outcome, covariate, instrument variables, treatment and propensity score}
#' \item{Plot}{A list that contains the ggplot elements for unifrom confidence band}
#'
#' @export
#'
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
#'
#' RESULT_P_hut<-p_hut_gen(outcome_name="Y",
#'                                  covariate_name=c("X1","X2","X3","X4"),
#'                                  instrument_name=c("Z1","Z2","Z3","Z4"),
#'                                  treatment_name= "D",
#'                                  data=demo,
#'                                  family="probit",
#'                                  trim=0.01,
#'                                  intercept=TRUE
#'                                 )
#'
#' RESULT_unif<-unif_gen(outcome_name="Y",
#' covariate_name=c("X1","X2","X3","X4"),
#' instrument_name=c("Z1","Z2","Z3","Z4"),
#' treatment_name= "D",
#' P_hut_name= "P_hut",
#' bw=0.1,
#' beta0=beta0,
#' beta1=beta1,
#' data=RESULT_P_hut$Data,
#' l_eval=0.15,
#' u_eval=0.85,
#' kernel=locpol::gaussK,
#' covariate_value=covariate_value,
#' significance_level=0.05)
#'}
#'
#'
#' @references
#' Tsuda, Jin and Okui (2024+).
#' Uniform Confidence Band for Marginal Treatment Effect Function.
#' will be available on arXiv

unif_gen<-function(outcome_name=outcome_name,
                   covariate_name=covariate_name,
                   instrument_name=instrument_name,
                   treatment_name=treatment_name,
                   P_hut_name=P_hut_name,
                   beta0=beta0,
                   beta1=beta1,
                   data=data,
                   bw=bw,
                   l_eval=l_eval,
                   u_eval=u_eval,
                   kernel=locpol::gaussK,
                   covariate_value=covariate_value,
                   significance_level=0.05){
  if((l_eval+0.05)> u_eval){
    stop(paste("Error: The values of evaluation points are incorrect"))
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

  P_hut_List<-P_hut_name

  SS<-nrow(data)

  P_hut_covariate<- data.frame(matrix(nrow=SS,ncol=length(covariate_name)))
  for (i in  1:length(covariate_name)){
    P_hut_covariate[,i]<-(data.covariate[[i]])*data[,P_hut_List]
  }
  colnames(P_hut_covariate) <- paste("P_",covariate_name,sep="")

  Dataframe<-cbind(data,P_hut_covariate)
  rownames(Dataframe)<-seq(1,SS,1)
  colnames(Dataframe)<-c(outcome_name, covariate_name, instrument_name,treatment_name, "P_hut",  paste("P_",covariate_name,sep=""))

  beta_hutc<-cbind(t(as.matrix(beta0)),t(as.matrix(beta1-beta0)))
  KExpl<-as.matrix(Dataframe[,c(covariate_name, paste("P_",covariate_name,sep=""))])
  Ytildac<-Dataframe[,outcome_name]-(KExpl %*% t(beta_hutc))

  #Estimation of first derivative of K(p)
  LAST_Dc<-cbind(Ytildac,data[,P_hut_List])
  colnames(LAST_Dc)<-c("Ytildac", "P_hut")
  LAST_Dcf<-as.data.frame(LAST_Dc)

  dimX<-length(covariate_name)
  dimX_L<-dimX+1
  dimX_R<-2*dimX

  ac<-t(as.matrix(beta_hutc[1,dimX_L:dimX_R])) %*% covariate_value

  Supp_P_hut_U<-u_eval
  Supp_P_hut_L<-l_eval
  O_SS<-1000
  eval<-seq(Supp_P_hut_L,Supp_P_hut_U,length=O_SS)

  Density_c<-locpol::PRDenEstC(LAST_Dcf$P_hut,bw=bw, kernel=kernel, xeval=eval)
  DDc<-locpol::locpol(Ytildac~P_hut,data=LAST_Dcf,bw=bw,kernel =kernel, xeval=eval,deg=2)
  BHC<-DDc$lpFit[,3]+rep(ac,length(eval))

  sigma_hut<-numeric(length(eval))
  for(l in 1:length(eval)){
    b<-(1/(4*sqrt(pi)))*((1/Density_c[l,2]))*((1/locpol::mu2K(kernel))^2)*DDc$lpFit[l,6]
    sigma_hut[l]<-b
  }


  length_supp<-max(eval)-min(eval)

  if(-2*log((bw*2*pi)/((length_supp)*sqrt(1.5)))<=0){
    stop(paste("Error: The critical value cannot be calculated"))
  }

  critical_value<-sqrt(-2*log((bw*2*pi)/((length_supp)*sqrt(1.5)))-2*log(-(1/2)*log(1-significance_level)))

  BHC_L<-numeric(length(eval))
  BHC_L<-BHC-critical_value*sqrt(sigma_hut/(SS*((bw)^3)))
  BHC_U<-numeric(length(eval))
  BHC_U<-BHC+critical_value*sqrt(sigma_hut/(SS*((bw)^3)))

  Plot_R<-cbind(BHC_L,BHC, BHC_U)

  Result_plot<-as.data.frame(cbind(Plot_R,eval))

  plot<-Result_plot%>%ggplot2::ggplot()
  plot<-plot+ggplot2::geom_ribbon(fill="red",ggplot2::aes(x =eval,
                                                 y = BHC,ymin =BHC, ymax =BHC_U), alpha = 0.3)
  plot<-plot+ggplot2::geom_ribbon(fill="red",ggplot2::aes(x = eval,
                                                 y = BHC,ymin =BHC_L, ymax =BHC), alpha = 0.3)
  plot<-plot+ggplot2::geom_point(ggplot2::aes(x = eval,
                                     y = BHC,color="Estimated MTE"),size = 0.8)
  plot<-plot+ggplot2::scale_color_manual(values = c("Estimated MTE" = "red"))
  plot<-plot+ggplot2::labs(x = "P(Z)",y=outcome_name)+ggplot2::theme(axis.title.x =ggplot2::element_text(size=15,hjust = 0.5))+ggplot2::theme(legend.position = "none")

  Eval_U<-0.05*floor(u_eval/0.05)
  Eval_L<-0.05*ceiling(l_eval/0.05)
  B2<-seq(Eval_L,Eval_U,by=0.05)


  Density_B2<-locpol::PRDenEstC(LAST_Dcf$P_hut,bw=bw, kernel=kernel, xeval=B2)
  DDB2<-locpol::locpol(Ytildac~P_hut,data=LAST_Dcf,bw=bw,kernel =kernel, xeval=B2,deg=2)
  B2_MTE<-DDB2$lpFit[,3]+rep(ac,length(B2))

  sigma_hut_B2<-numeric(length(B2))
  for(l in 1:length(B2)){
    b<-(1/(4*sqrt(pi)))*((1/Density_B2[l,2]))*((1/locpol::mu2K(kernel))^2)*DDB2$lpFit[l,6]
    sigma_hut_B2[l]<-b
  }



  BHC_L_B2<-numeric(length(B2))
  BHC_L_B2<-B2_MTE-critical_value*sqrt(sigma_hut_B2/(SS*((bw)^3)))
  BHC_U_B2<-numeric(length(eval))
  BHC_U_B2<-B2_MTE+critical_value*sqrt(sigma_hut_B2/(SS*((bw)^3)))


  Estimation<-as.data.frame(cbind(B2,B2_MTE, sqrt(sigma_hut_B2), BHC_L_B2,BHC_U_B2))
  colnames(Estimation)<-c("P_hut","MTE", "SE", "UB of UCB","LB of UCB")
  return(list(Estimation=Estimation, Plot=plot))
}
