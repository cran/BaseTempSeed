
#' @title Preparation of Data for Germination Percentages
#'
#' @param Data Dataframe with Temperature in Row and Time in Column
#' @param TSeed Total Number of Seeds
#' @param GPercentage Vector of Germination Fractions
#' @import stats
#' @return
#' \itemize{
#'   \item InterData: Final Data for Estimation
#' }
#' @export
#'
#' @examples
#' library(BaseTempSeed)
#' Interdata<-TempInter(Data=Data,GPercentage=c(0.1,0.2, 0.3,0.4, 0.5,0.6, 0.7,.8,0.9))
#'
#' @references
#' \itemize{
#'\item Ellis, R. H., Simon, G., & Covell, S. (1987). The Influence of Temperature on Seed Germination Rate in Grain LegumesIII. A Comparison of Five Faba Bean Genotypes at Constant Temperatures Using a New Screening Method. Journal of Experimental Botany, 38(6), 1033–1043.
#'\item Garcia-Huidobro, J., Monteith, J. L., & Squire, G. R. (1982). Time, temperature and germination of pearl millet (Pennisetum typhoides S. & H.) I. Constant temperature. Journal of experimental botany, 33(2), 288-296.
#'}

TempInter<-function(Data,TSeed=50,GPercentage=c(0.1, 0.25, 0.5, 0.75, 0.9)){
  data <- as.data.frame(Data)
  tot_seed <- TSeed# Total number of seeds taken in petridish
  g <- GPercentage # Germination percentages to be interpolated


  final <- NULL
  for (n in 1:nrow(data)) {
    .in <- as.vector(data[n, -1])
    times <- as.numeric(names(.in))
    cum_sum <- cumsum(.in)
    p_germ <- cum_sum/tot_seed

    g_inter <- subset(g, g<=max(p_germ))

    values <- approx(p_germ, times, g_inter, ties = "max")
    inter_time <- values$y
    na_data <- rep(NA, length(g)-length(inter_time))
    inter_final <- append(inter_time, na_data)
    final <- rbind(final, inter_final)
  }
  colnames(final) <- paste0(g*100,"%")
  rownames(final) <- data[, 1]
  InterData<-final
  return(InterData)
}

#' @title Estimation of Base Temperature
#'
#' @param Data Output TempInter Function
#' @import stats NlcOptim
#' @return
#' \itemize{
#'   \item Coefficients: Estimate of Coefficients
#' }
#' @export
#'
#' @examples
#' library("BaseTempSeed")
#' Inter_data <- TempInter(Data=Data,GPercentage=c(0.1,0.2, 0.3,0.4, 0.5,0.6, 0.7,.8,0.9))
#' Est<-EstCoeff(Inter_data)
#' @references
#' \itemize{
#'\item Ellis, R. H., Simon, G., & Covell, S. (1987). The Influence of Temperature on Seed Germination Rate in Grain LegumesIII. A Comparison of Five Faba Bean Genotypes at Constant Temperatures Using a New Screening Method. Journal of Experimental Botany, 38(6), 1033–1043.
#'\item Garcia-Huidobro, J., Monteith, J. L., & Squire, G. R. (1982). Time, temperature and germination of pearl millet (Pennisetum typhoides S. & H.) I. Constant temperature. Journal of experimental botany, 33(2), 288-296.
#'}

EstCoeff<-function(Data){
  Inter_data<-Data
  int<-t(na.omit(t(Inter_data)))
  if(ncol(int)!=ncol(Inter_data)){
    message(paste0("Interpolation at fraction ",names(attr(int,"na.action"))," has missing value at ","temperature"," ",as.vector(attr(int,"na.action"))," and this fraction is/are ommited"))
    Inter_data<-Inter_data[,!(colnames(Inter_data) %in% (names(attr(int,"na.action"))))]
  }
  reg_data <- 1/Inter_data
  Coeff <- NULL
  for (c in 1:ncol(reg_data)) {
    # c=1
    reg<-lm(reg_data[,c]~as.numeric(row.names(reg_data)))
    Coef<-as.vector(reg$coefficients)
    names(Coef)<-c("Alpha", "Beta")
    Coeff<-rbind(Coeff,Coef)
  }

  # x_write
  coeff_vector<-as.vector(t(Coeff))
  name<-paste0("x","[",rep(seq(1,2*ncol(reg_data))),"]")
  x_name<-NULL
  for (i in seq(1,length(name),2)) {
    R<-rep(name[i:(i+1)],nrow(reg_data))
    x_name<-append(x_name,R)
  }



  y_data<-as.vector(reg_data)
  temperature<-rep(as.numeric(row.names(reg_data)),ncol(reg_data))
  e=1
  a=2*(e-1)+1
  b=a+1
  equation<-(paste0("(",y_data[e],"-",noquote(x_name[a]),"-",temperature[e],"*",noquote(x_name[b]),")^2"))
  for (e in 2:length(y_data)) {
    a=2*(e-1)+1
    b=a+1
    eq<-(paste0("(",y_data[e],"-",noquote(x_name[a]),"-",temperature[e],"*",noquote(x_name[b]),")^2"))
    equation<-paste0(equation,"+",eq)
  }
  objfunc <- function(x) {
    return(eval(parse(text=equation)
    ))
  }


  f <- NULL
  for (p in 2:ncol(reg_data)) {
    f <- rbind(f, paste0(paste0("x","[1]"),"*",paste0("x","[",2*p,"]")," -", paste0("x","[",2,"]"),"*",paste0("x","[",(2*p-1),"]")
    ))
  }
  constr <- function(x) {
    f
    return(list(ceq=eval((parse(text=f))),c=NULL))
  }

  sol <- solnl(coeff_vector, objfun = objfunc, confun = constr)
  estimate<-as.vector(sol$par)
  alpha_est_seq<- (estimate[seq_len(length(estimate)) %% 2==1])
  beta_est_seq<-(estimate[seq_len(length(estimate)) %% 2==0])
  tb<- mean(-alpha_est_seq/beta_est_seq)
  beta_est_seq1<- -alpha_est_seq/tb
  Coefficients<-data.frame(Aplha=alpha_est_seq, Beta=beta_est_seq1)
  rownames(Coefficients)<-colnames(reg_data)
  return(Coefficients)
}

#' This is data to be included in my package
#' @name Data
#' @docType data
#' @keywords datasets
#' @usage data(Data)
#' @format A data frame with 5 rows and 19 column
NULL



