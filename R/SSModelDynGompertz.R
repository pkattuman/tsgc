# Created by: Craig Thamotheram
# Created on: 27/07/2022

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

setOldClass("xts")
setOldClass("KFS")
#'
#' @title Class for Dynamic Gompertz Curve State-Space Model
#'
#' @description Class for Dynamic Gompertz Curve State-Space Model Object, whose 
#' fields contain the model settings. Contains methods to obtain FilterResults
#' object and plot the time series.
#'
#' @field Y The cumulated variable.
#' @field q The signal-to-noise ratio (ratio of slope to irregular
#'   variance). Defaults to \code{'NULL'}, in which case no
#'   signal-to-noise ratio will be imposed. Instead, it will be estimated.
#' @field sea.type Seasonal type. Options are \code{'trigonometric'}
#'   and \code{'none'}. \code{'trigonometric'} will yield a model with a
#'   trigonometric seasonal component and \code{'none'} will yield a model
#'   with no seasonal component.
#'@field sea.period The period of seasonality. For a day-of-the-week
#'   effect with daily data, this would be 7. Not required if
#'   \code{sea.type = 'none'}.
#'@field reinit.date (Only needed for reinitialization.) The reinitialisation date \eqn{r}. Should be
#' specified as an object of class \code{"Date"}. Defaults to NA, which 
#' represents the non-reinitialized version.
#' @field original.results (Only needed for reinitialization.) Rather than re-estimating the model up
#' to the \code{reinit.date}, a \code{FilterResults} class object can be
#' specified here and the parameters for the reinitialisation will be taken
#' from this object. Default is \code{NULL}. This parameter is optional.
#' @field use.presample.info (Only needed for reinitialization.) Logical value denoting whether or
#' not to use information from before the reinitialisation date in the
#' reinitialisation procedure. Default is \code{TRUE}. If \code{FALSE}, the
#' model is estimated from scratch from the reinitialisation date and no
#' attempt to use information from before the reinitialisation date is made.
#'
#' 
#' @importFrom xts periodicity last
#' @importFrom methods new
#' @importFrom xts xts
#' @importFrom zoo index
#' @importFrom KFAS SSModel fitSSM KFS
#' @importFrom magrittr %>%
#' @importFrom KFAS SSMtrend SSMseasonal
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-06")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' @export SSModelDynamicGompertz
#' @exportClass SSModelDynamicGompertz
SSModelDynamicGompertz <- setRefClass(
  "SSModelDynamicGompertz",
  fields = list(
    Y = "xts",
    q = "ANY",
    sea.type="ANY", 
    sea.period="ANY",
    reinit.date = "ANY",
    original.results = "ANY",
    use.presample.info = "ANY"
  ),
  methods = list(initialize = function(Y, q = NULL, sea.type = 'trigonometric',
                                       sea.period = 7,reinit.date=NULL, original.results=NULL,
                                       use.presample.info=TRUE)
  {
    "Create an instance of the \\code{SSModelDynamicGompertz} class. Parameters 
    are defined in `fields` section. 
      \\subsection{Usage}{\\code{SSModelDynGompertzReinit$new(y, q = 0.005,
      reinit.date = as.Date(\"2021-05-12\",format = date.format))}}"
    Y <<- Y
    q <<- q
    sea.type <<- sea.type
    sea.period <<- sea.period
    reinit.date <<- reinit.date
    original.results <<- original.results
    use.presample.info <<- use.presample.info
  }, 
  update = function(pars, model, q, sea.type) {
    "Update method for Kalman filter to implement the dynamic Gompertz curve
       model.
       A maximum of 3 parameters are used to set the observation noise
       (1 parameter), the transition equation slope and seasonal noise. If q (signal
        to noise ratio) is not null then the slope noise is set using this
        ratio.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{pars} Vector of parameters.}
        \\item{\\code{model} \\code{KFS} model object.}
        \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
         variance).}
        \\item{\\code{sea.type} Seasonal type. Options are
        \\code{'trigonometric'} and \\code{'none'}.}
      }}
      \\subsection{Return Value}{\\code{KFS} model object.}"
    estH <- any(is.na(model$H))
    estQ <- any(is.na(model$Q))
    if ((!estH) & (!estQ)) {
      # If nothing to update then return model
      return(model)
    } else {
      nparQ <- if (sea.type == 'trigonometric') { 1 } else { 0 }
      # 1. Set seasonal noise
      if (estQ) {
        Q <- as.matrix(model$Q[, , 1])
        # Update diagonal elements
        naQd <- which(is.na(diag(Q)))
        Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
        diag(Q)[naQd] <- exp(0.5 * pars[1])
        # Check for off-diagonal elements and raise error if found.
        naQnd <- which(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd]))
        if (length(naQnd) > 0) {
          stop("NotImplmentedError: Unexpected off-diaganol element updating")
        }
      }
      
      # 2. Set observation noise
      H <- as.matrix(model$H[, , 1])
      if (estH) {
        naHd <- which(is.na(diag(H)))
        H[naHd, naHd][lower.tri(H[naHd, naHd])] <- 0
        diag(H)[naHd] <- exp(0.5 * pars[(nparQ + 1)])
        model$H[naHd, naHd, 1] <- crossprod(H[naHd, naHd])
      }
      
      # 3. Set slope noise
      # Get index of slope, 1 before the seasonal component.
      model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
      i.slope <- 2
      # Estimate slope if no signal to noise ratio specified.
      if (is.null(q)) {
        Q.slope <- exp(0.5 * pars[(nparQ + 2)])
      } else {
        Q.slope <- crossprod(H[naHd, naHd]) * q
      }
      model$Q[i.slope, i.slope, 1] <- Q.slope
    }
    return(model)
  },
  estimate = function() {
    "Estimates the dynamic Gompertz curve model when applied to an object of
      class \\code{SSModelDynamicGompertz}.
      \\subsection{Return Value}{An object of class \\code{FilterResults}
      containing the result output for the estimated dynamic Gompertz curve
      model.}
      "
    # 1. Get LDL of cumulative series Y.
    y <- tsgc::df2ldl(Y)
    
    # 2. Add update / model methods
    updatefn <- purrr::partial(
      .self$update, ... =, q = q, sea.type = sea.type
    )
    model <- .self$get_model(y, q = q, sea.type, sea.period)
    # 2. Estimate via MLE unknown params
    model_fit <- fitSSM(model$model, inits = model$inits, updatefn = updatefn,
                        method = 'BFGS')
    
    # 3. Run smoother/filter
    model_output <- KFS(model_fit$model)
    
    # 4. Get truncated index from model if using a reinitialisation in
    # self$get_model
    date.index <- if (!is.null(model$index)) { model$index } else { index(y) }
    
    results <- FilterResults$new(
      data_xts = Y,
      index = date.index,
      reinit.date=reinit.date,
      output = model_output
    )
    return(results)
  },
  get_dynamic_gompertz_model = function(
    y,
    q = NULL,
    sea.type = 'trigonometric',
    sea.period = 7,
    a1 = NULL,
    P1 = NULL,
    Q = NULL,
    H = NULL)
  {
    "Returns the Dynamic Gompertz Model object, assuming no reinitialization. 
    Used in the \\code{get_model} function.
    \\subsection{Parameters}{\\itemize{
      \\item{\\code{y} The cumulated variable}
      \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
      variance). Defaults to \\code{'NULL'}, in which case no signal-to-noise
      ratio will be imposed. Instead, it will be estimated.}
      \\item{\\code{sea.type} Seasonal type. Options are \\code{'trigonometric'}
       and \\code{'none'}. \\code{'trigonometric'} will yield a model with a
       trigonometric seasonal component and \\code{'none'} will yield a model
       with no seasonal component.}
      \\item{\\code{sea.period} The period of seasonality. For a day-of-the-week
       effect with daily data, this would be 7. Not required if
       \\code{sea.type = 'none'}.}
      \\item{\\code{a1} Optional parameter specifying the prior mean of the
      states. Defaults to \\code{'NULL'}. Leave as \\code{'NULL'} for a diffuse
      prior (no prior information). If a proper prior is to be specified, both
      \\code{a1} and \\code{P1} must be given.}
      \\item{\\code{P1} Optional parameter specifying the prior mean of the
      states. Defaults to \\code{'NULL'}. Leave as \\code{'NULL'} for a diffuse
       prior (no prior information). If a proper prior is to be specified,
       both \\code{a1} and \\code{P1} must be given.}
      \\item{\\code{Q} Optional parameter specifying the state error variances
      where these are to be imposed rather than estimated. Defaults to
      \\code{'NULL'} which will see the variances estimated.}
      \\item{\\code{H} Optional parameter specifying the irregular variance
      where this is to be imposed rather than estimated. Defaults to
      \\code{'NULL'} which will see the variance estimated.}
    }}
    "
    Qt.slope <- if (is.null(Q)) { NA } else { Q[2, 2] }
    Qt.seas <- if (is.null(Q)) { NA } else { Q[3, 3] }
    Ht <- if (is.null(H)) { NA } else { H }
    
    # 1. Set prior on state as ~ N(a1, P1) if a1 supplied.
    use.prior <- if (!is.null(a1)) { TRUE } else { FALSE }
    
    if (use.prior) {
      if (sea.type == 'trigonometric') {
        ss_model <- SSModel(
          y ~
            SSMtrend(
              degree = 2,
              Q = list(matrix(0), matrix(Qt.slope)),
              a1 = a1[1:2],
              P1 = P1[1:2, 1:2]
            ) +
            SSMseasonal(
              period = sea.period,
              Q = Qt.seas,
              sea.type = sea.type,
              a1 = a1[3:dim(a1)[1]],
              P1 = P1[3:dim(a1)[1], 3:dim(a1)[1]]
            ),
          H = Ht
        )
        n.pars <- 0
      } else if (sea.type == 'none') {
        ss_model <- SSModel(
          y ~
            SSMtrend(
              degree = 2,
              Q = list(matrix(0), matrix(Qt.slope)),
              a1 = a1[1:2],
              P1 = P1[1:2, 1:2]
            ),
          H = Ht
        )
        n.pars <- 0
      } else {
        stop(sprintf("sea.type= '%s' not implemented", sea.type))
      }
    } else {
      if (sea.type == 'trigonometric') {
        ss_model <- SSModel(
          y ~
            SSMtrend(
              degree = 2,
              Q = list(matrix(0), matrix(Qt.slope))
            ) +
            SSMseasonal(
              period = sea.period,
              Q = Qt.seas,
              sea.type = sea.type),
          H = matrix(Ht)
        )
        n.pars <- if (is.null(q)) { 3 } else { 2 }
      } else if (sea.type == 'none') {
        ss_model <- SSModel(
          y ~
            SSMtrend(
              degree = 2,
              Q = list(matrix(0), matrix(Qt.slope))
            ),
          H = matrix(Ht)
        )
        n.pars <- if (is.null(q)) { 2 } else { 1 }
      } else {
        stop(sprintf("sea.type= '%s' not implemented", sea.type))
      }
    }
    out <- list(model = ss_model, inits = rep(0, n.pars))
    return(out)
  },
    get_model = function(
    y,
    q = NULL,
    sea.type = 'trigonometric',
    sea.period = 7
    )
    {
        "Returns dynamic Gompertz curve model. 
    \\subsection{Parameters}{\\itemize{
      \\item{\\code{y} The log-growth rate time series.}
      \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
      variance). Defaults to \\code{'NULL'}, in which case no signal-to-noise
      ratio will be imposed. Instead, it will be estimated.}
      \\item{\\code{sea.type} Seasonal type. Options are \\code{'trigonometric'}
       and \\code{'none'}. \\code{'trigonometric'} will yield a model with a
       trigonometric seasonal component and \\code{'none'} will yield a model
       with no seasonal component.}
      \\item{\\code{sea.period} The period of seasonality. For a day-of-the-week
       effect with daily data, this would be 7. Not required if
       \\code{sea.type = 'none'}.}
    }}
    \\subsection{Description}{
    The dynamic Gompertz with an integrated random walk (IRW) trend is
    \\deqn{\\ln g_{t}=\\delta_{t}+\\varepsilon_{t},  \\;\\;\\;\\;
    \\varepsilon_{t}\\sim NID(0,\\sigma_{\\varepsilon }^{2}), \\;\\;\\;\\;
    t=2,...,T, }
    where \\eqn{Y_t} is the cumulated variable, \\eqn{y_t = \\Delta Y_t},
    \\eqn{\\ln g_{t}=\\ln y_{t}-\\ln Y_{t-1}} and
    \\deqn{\\delta_{t} =\\delta_{t-1}+\\gamma_{t-1},}
    \\deqn{\\gamma_{t} =\\gamma_{t-1}+\\zeta_{t}, \\;\\;\\;\\;
    \\zeta_{t}\\sim NID(0,\\sigma_{\\zeta }^{2}),}
    where the observation disturbances \\eqn{\\varepsilon_{t}}  and slope
    disturbances \\eqn{\\zeta_{t}}, are iid Normal and mutually independent.
    Note that, the larger the signal-to-noise ratio,
    \\eqn{q_{\\zeta }=\\sigma_{\\zeta }^{2}/\\sigma_{\\varepsilon }^{2}},
    the faster the slope changes in response to new observations. Conversely,
    a lower signal-to-noise ratio induces smoothness.

    For the model without seasonal terms (\\code{sea.type = 'none'}) the are
    priors are
    \\deqn{\\begin{pmatrix} \\delta_1 \\ \\gamma_1 \\end{pmatrix}
    \\sim N(a_1,P_1)}.
    The diffuse prior has \\eqn{P_1 = \\kappa I_{2\\times 2}} with
    \\eqn{\\kappa \\to \\infty}. Implementation of the diffuse prior is handled
     by the package \\code{KFAS} (Helske, 2017). Where the model has a seasonal
      component (\\code{sea.type = 'trigonometric'}), the vector of prior means
       \\eqn{a_1} and the prior covariance matrix \\eqn{P_1} are extended
       accordingly.

    See the vignette for details of the variance matrix \\eqn{Q}.
    \\eqn{H = \\sigma^2_{\\varepsilon}}.
    }
      \\subsection{Return Value}{\\code{KFS} model object.}"
      if (is.null(.self$reinit.date)){
        model <- .self$get_dynamic_gompertz_model(
          y, q = q, sea.type = sea.type, sea.period = sea.period
        )
        return(model)
      }
      else{
        # 4.1. Index for reinitialisation, t_0
        stopifnot(length(Y[reinit.date]) == 1)
        Y.t.r_0 <- as.numeric(Y[reinit.date - 1])
        
        # 4.2 Reinitialisation:
        #   ln g_t^r = ln g_t + ln (Y_{t-1}/Y_{t-1}^r), where Y_t^r=Y_t-Y_{r_0}.
        idx.dates <- (index(y) > reinit.date)
        lag.Y <- stats::lag(Y)[idx.dates]
        y.reinit <- y[index(y) > reinit.date] + log(lag.Y / (lag.Y - Y.t.r_0))
        
        # 4.3 Run Kalman filter/smoother on new series with non-diffuse prior
        if (use.presample.info) {
          # Either estimate full model here or take results from previous model.
          if (is.null(original.results)) {
            # NB. Restrict sample to t<=r - date of reinitialisation.
            idx.est <- zoo::index(Y) <= reinit.date
            model <- SSModelDynamicGompertz$new(Y = Y[idx.est], q = q)
            res.original <- model$estimate()
            model_output <- output(res.original)
          } else {
            model_output <- output(original.results)
          }
          model_seasonal <- seasonalComp(model_output)
          season.type <- if (is.null(model_seasonal)) {'none'} else {
            'trigonometric'}
          season.period <- if (!is.null(model_seasonal)) {
            ncol(att(model_output)) - 1}
          
          # 4.3 Reset slope to 0 and add constant to initial value for level.
          # where reinit.date is t=r
          idx <- which(reinit.date == index(y))
          stopifnot(length(idx) == 1)
          att <- att(model_output)[idx,]
          Ptt <- Ptt(model_output)[, , idx]
          Tt <- drop(matrixKFS(model_output,"T"))
          Rt <- drop(matrixKFS(model_output,"R"))
          Qt <- drop(matrixKFS(model_output,"Q"))
          Ht <- drop(matrixKFS(model_output,"H"))
          
          # a. Take a_{r|r} and P_{r|r} through prediction step to get a_{r+1}
          # and P_{r+1}
          a1 <- Tt %*% att
          P1 <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
          
          # b. Set slope to 0 and add correction (\ln(Y_r/y_r) to level.
          a1["slope",] <- 0
          a1["level",] <- a1["level",] + log(Y[idx] / (Y[idx] - Y.t.r_0))
        } else {
          # Don't use presample info
          a1 <- NULL; P1 <- NULL; Qt <- NULL; Ht <- NULL
        }
        out <- .self$get_dynamic_gompertz_model(
          y = y.reinit, q = q, sea.type = season.type, sea.period = season.period,
          a1 = a1, P1 = P1, Q = Qt, H = Ht
        )
        out[['index']] <- index(y.reinit)
        return(out)
      }
    },
    summary = function() {
      "Supplies details of the SSModelDynamicGompertz object, such as estimated 
      parameter values, start and end dates of estimation."
      out <- suppressWarnings(output(.self$estimate()))
      if(is.null(q)){
        qest <- matrixKFS(out,"Q")[2, 2, 1]/matrixKFS(out,"H")[, , 1]
      }
      reinit<-!is.null(reinit.date)
      dates<-index(Y)
      
      cat("Summary of SSModelDynamicGompertz Model")
      if (reinit) {
        cat(" (Reinitialized)")
      }
      cat("\n")
      cat("--------------------------------------\n")
      cat("Cumulated Variable:\n")
      base::print(head(.self$Y))
      cat("Signal-to-Noise Ratio (q):", 
          ifelse(is.null(q), paste(signif(qest,3), "(estimated)"), 
                 paste(q, ("(user specified)"))), "\n")
      cat("Model Details:\n")
      cat("  - Model Type: Dynamic Gompertz Curve")
      if (reinit) {
        cat(" (Reinitialized)")
      }
      cat("\n")
      cat("  - Seasonal Component: ", ifelse(sea.type == 'none', "None", "Trigonometric"), "\n")
      cat("  - Period of Seasonality: ", ifelse(sea.type == 'none', "N/A", sea.period), "\n")
      cat("  - Dataset start date:", format(as.Date(dates[1], origin = "1970-01-01")))
      cat("\n")
      cat("  - Dataset end date:", format(as.Date(tail(dates,1), origin = "1970-01-01")))
      cat("\n")
      if (reinit){
        cat("  - Reinitialization date:",format(as.Date(reinit.date, origin = "1970-01-01")))
        cat("\n")
        cat("  - Use presample info:", use.presample.info)
        cat("\n")
      }
      cat("  - Model States and Standard Errors\n")
      base::print(out)
    },
    print = function() {
      "Provides a quick description of SSModelDynamicGompertz object, providing 
      model states and standard errors."
      reinit<-!is.null(reinit.date)
      out <- output(.self$estimate()) #KFS object
      if(is.null(q)){
        qest <- matrixKFS(out,"Q")[2, 2, 1]/matrixKFS(out,"H")[, , 1]
      }
      cat("SSModelDynamicGompertz Model")
      if (reinit) {
        cat(" (Reinitialized)")
      }
      cat("\n")
      cat("\n")
      cat("Cumulated Variable:\n")
      base::print(head(Y))
      cat("Number of observations:", length(.self$Y))
      cat("\n")
      cat("Signal-to-Noise Ratio (q):", 
          ifelse(is.null(q), paste(signif(qest,5), "(estimated)"), 
                 paste(q, ("(user specified)"))), "\n")
      cat("Seasonal components?",
          ifelse(is.null(seasonalComp(out)),
                 "No","Yes"),"\n")
      if (!is.null(reinit.date)){
        cat("Reinit date:",format(as.Date(reinit.date, origin = "1970-01-01")))
        cat("\n")
        cat("Use presample info:", use.presample.info)
      }
    },
    plot =function(title=NULL, series.name="target variable", MA=TRUE){
      "Plots the lagged differences of the cumulated dataset \\code{Y} in this 
      \\code{SSModelDynamicGompertz} object against time, which could represent 
      daily cases.
      \\subsection{Parameters}{\\itemize{
        \\item{\\code{MA} A logical value indicating whether 7-day centered moving 
        average should be plotted. Defaults to \\code{TRUE}.
        }
        \\item{\\code{title} Title for forecast plot. Enter as text string. 
        \\code{NULL} (i.e. no title) by default.}
        \\item{\\code{series.name} The name of the series the growth rate is being computed
        for. E.g. \\code{'cases'}. Default is `target variable`.}
}
      }
      \\subsection{Return Value}{A plot of the lagged differences of the 
      cumulated dataset \\code{Y} against time.}
      "
      cumulative_cases <- Y  

      # Calculate a centred 7-day moving average of daily differences.
      ma.cent.new.cases <- zoo::rollmean(diff(cumulative_cases), 7, align = "center")
      
      # Identify the date with maximum new cases.
      ma.cent.wave.3.idx.max <- tsgc::argmax(ma.cent.new.cases) %>% zoo::index()
      
      # Prepare data for plotting by combining actual new cases and the moving average.
      d <- cbind(diff(cumulative_cases), ma.cent.new.cases)
      colnames(d) <- c('New Cases', 'Centered 7-day MA')
      d.df <- data.frame(
        Date = index(d),
        New.Cases = coredata(d[, 1]),
        Centered.7.day.MA = coredata(d[, 2])
      )
      
      # Create base plot
      data_plot <- ggplot(data = d.df, aes(x = Date)) +
        geom_line(aes(y = New.Cases, color = "New Cases"), linewidth = 0.1) +
        scale_y_continuous(n.breaks = 10) +
        labs(x = "Date", y = paste("New", series.name), title = title)+
        scale_x_date(date_breaks = "60 days") +
        theme_light(base_size = 12) +
        theme(
          legend.position = "inside",
          legend.position.inside = c(0.2, 0.85),
          legend.title = element_text(size = 2),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(face = "bold")
        )
      
      # Conditionally add Centered 7-day MA line
      if (MA) {
        data_plot <- data_plot + 
          geom_line(aes(y = Centered.7.day.MA, color = "Centered 7 day MA"), linewidth = 1)+ 
          scale_color_manual(
            name = '',
            values = c('Centered 7 day MA' = 'red')
          )
      } 
      data_plot
    }
  )
)
