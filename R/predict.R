##' Posterior simulations of bayesian binomial logistic model as defined in Kogan 2022
##' @details
##' Flattens the 3d simulation array from stanfit to a 2d array for analysis. Stores the array and removes the model object.
##' @return
##' Object of class chbmod containing
##' \item{model_pars}{unique names of model parameters}
##' \item{rhsform}{character model formula for RHS}
##' \item{sflat}{2d array of posterior simulations nsims x npars}
#' @export
chbmod <- function(rhsform, fteed, model, vty, fit, shrink = T, comments = "", ...) {
  require(shedder)
  require(digest)
  fteed_vty <- fteed %>% filter(variety == vty)
  object <- fit(rhsform, fteed_vty, model, ...)
  if(shrink) {
    object <- object %>%
      shredder::stan_axe(what = 'fit_instance') %>%
      shredder::stan_axe(what = 'stanmodel')
  }
  sflat <- flatten_stan_array(as.array(object))
  pars <- dimnames(sflat)[[2]]
  par_grps <- list(
    beta = pars %>% str_subset("beta"),
    sd = pars %>% str_subset("S_"),
    seasonXfield = pars %>% str_subset("r_seasonXfield"),
    fieldXdate_fi = pars %>% str_subset("r_fieldXdate_fi")
    )
  par_grps[["fieldXdate"]] <- setdiff(pars, par_grps[["fieldXdate_fi"]]) %>% str_subset("r_fieldXdate")

  structure(list(
    comments = comments,
    datahash = digest::digest(fteed_vty),
    date = list(pretty=format(Sys.time(), "%a %b %d %X %Y"),formatted=str_replace_all(format(Sys.time(), "%Y-%m-%d-%X"), ":", "") %>% str_replace_all(" ", "")),
    model_pars = object@model_pars,
    model_pars_long = pars,
    par_grps = par_grps,
    rhsform = rhsform,
    sflat = sflat,
    stanfit = object,
    vty = vty),
    class = "chbmod")
}

#' @export
saveChbmod <- function(obj, ...) UseMethod("saveChbmod")

#' @export
saveChbmod.cherry_s1 <- function(obj, folder = here("dataProcessed")) {
  saveRDS(obj, file = paste0(folder, "/stage1_model_bayes_", obj$vty, "_", obj$date$formatted,".RDS"))
}

#' @export
saveChbmod.cherry_s2 <- function(obj, folder = here("dataProcessed")) {
  saveRDS(obj, file = paste0(folder, "/model_bayes_", obj$vty, "_", obj$date$formatted,".RDS"))
}

#' @export
saveChbmod.blueberry_s1 <- function(obj, folder = here("dataProcessed")) {
  saveRDS(obj, file = paste0(folder, "/stage1_model_bayes_", obj$vty, "_", obj$date$formatted,".RDS"))
}

#' @export
saveChbmod.blueberry_s2 <- function(obj, folder = here("dataProcessed")) {
  saveRDS(obj, file = paste0(folder, "/model_bayes_", obj$vty, "_", obj$date$formatted,".RDS"))
}
#' @export
print.chbmod <- function(obj, ...) {
  print(obj$summary)
  invisible(obj)
}

#' @export
cherry_s1_chbmod <- function(...) {
  rhsform <- "~std_ftemp"
  object <- chbmod(..., rhsform = rhsform, fit = fitting_cherry_stan)
  class(object) <- c("cherry_s1", class(object))
  object
}

#' @export
cherry_s2_chbmod <- function(...) {
  ns_knots <- c(-0.8, -0.45, -0.32, 0.25, 1.5)
  ns_bknots <- c(-1.18, 3.71)
  fe_formula <- as.formula(substitute(cbind(NoFlowersLive, NoFlowersDead)~std_ftemp + std_AIR_TEMP_F48 + ns(std_gddchill, knots = knots, Boundary.knots = bknots), list(knots = ns_knots, bknots = ns_bknots)))
  rhsform <- as.character(RHSForm(fe_formula, as.form = T))
  object <- chbmod(..., rhsform = rhsform, fit = fitting_cherry_stan)
  class(object) <- c("cherry_s2", class(object))
  object
}

#' @export
blueberry_s1_chbmod <-  function(...) {
  object <- chbmod(..., fit = fitting_blueberry_stan)
  class(object) <- c("blueberry_s1", class(object))
  object
}

#' @export
blueberry_s2_chbmod <-  function(...) {
  object <- chbmod(..., fit = fitting_blueberry_stan)
  class(object) <- c("blueberry_s2", class(object))
  object
}

#' @export
dataset <- function(obj, ...) UseMethod("dataset")

#' @export
dataset.chbmod <- function(obj, repo = here("reuse")) {
  name <- paste0("md5_", obj$datahash, ".RDS")
  files <- list.files(repo)
  fmatch <- match(name, files)
  readRDS(paste0(repo, "/", files[fmatch]))
}

#' @export
getformula <- function(obj, ...) {
  as.formula(obj$rhsform)
}

##' Fit model
##' @param object a \code{fe_formula_rhs} formula for model X matrix
##' @param object a \code{fteed_vty} data for X matrix
##' @param object a \code{model} compiled code for stan model
fitting_cherry_stan <- function(fe_formula_rhs, fteed_vty, model, ...) {

  fteed_vty_nc <- fteed_vty %>%
    filter(ftemp < 4) %>%
    mutate(
      can_ID = makeID(Field, date, ftemp), #
      Spur_ID = makeID(Field, date, ftemp, Spur),
      bud_ID = makeID(Field, date, ftemp, Spur, bud)
    )

  fteed_vty_c <- fteed_vty %>%
    filter(ftemp == 4) %>%
    mutate(
      Spur_ID = makeID(Field, date, ftemp, Spur),
      bud_ID = makeID(Field, date, ftemp, Spur, bud)
    )

  X <- model.matrix(formula(paste(fe_formula_rhs, collapse=" ")) , data = fteed_vty_nc)
  xcol <- ncol(X)
  colnames(X) <- paste("V", 1:ncol(X), sep = "")


  data_list <- list(
    P = xcol,
    N = nrow(fteed_vty_nc),
    NC = nrow(fteed_vty_c),
    trials = fteed_vty_nc$NoFlowers,
    trials_control = fteed_vty_c$NoFlowers,
    X = X,
    N_seasonXfield = max(fteed_vty_nc$seasonXfield_ID),
    seasonXfield_ID = fteed_vty_nc$seasonXfield_ID,
    N_fieldXdate = max(fteed_vty_nc$fieldXdate_ID),
    fieldXdate_ID = fteed_vty_nc$fieldXdate_ID,
    fieldXdate_control_ID = fteed_vty_c$fieldXdate_ID,
    N_can = max(fteed_vty_nc$can_ID),
    can_ID = fteed_vty_nc$can_ID,
    N_spur = max(fteed_vty_nc$Spur_ID),
    N_spur_control = max(fteed_vty_c$Spur_ID),
    spur_ID = fteed_vty_nc$Spur_ID,
    spur_control_ID = fteed_vty_c$Spur_ID,
    N_bud = max(fteed_vty_nc$bud_ID),
    N_bud_control = max(fteed_vty_c$bud_ID),
    bud_ID = fteed_vty_nc$bud_ID,
    bud_control_ID = fteed_vty_c$bud_ID,
    Y = fteed_vty_nc$NoFlowersLive,
    Y_control = fteed_vty_c$NoFlowersLive,
    x_cov = diag(10, xcol),# 100*as.matrix(vcov(glmerfit)),
    beta_mu = rep(0, xcol),#as.vector(coef(summary(glmerfit))[,1]),
    prior_only = F
  )
  # f_init <- function() {
  #   list(beta = MASS::mvrnorm(1, as.vector(coef(summary(glmerfit))[,1]),as.matrix(vcov(glmerfit))))
  # }
  fitted <- sampling(model, data = data_list, pars = c("S_seasonXfield","S_can", "S_Spur", "S_bud", "S_fieldXdate", "S_fieldXdate_fi", "beta", "r_seasonXfield","r_fieldXdate","r_fieldXdate_fi"), ...) #,  control = list( max_treedepth = 13),sample_file = "samples_100.csv", , control = list(adapt_delta = 0.999, max_treedepth = 13)
  fitted
}

##' posterior predictive distribution (or characteristics of) for chbmod
##' @param object a \code{chbmod} object
##' @param newdata new data for prediction
##' @param re.form  \code{NA} to specify population-level predictions (i.e., setting all random effects to zero) or a character vector with values as the names of chbmod parameters and names as the corresponding identifier in newdata
##' @param type \describe{
##' \item{"link"}{conditional mean on the scale of the link function,
##' or equivalently the linear predictor of the conditional model}
##' \item{"response"}{expected value}
##' }
##' @param stat \code{NULL} to specify posterior predictive distribution or function to specify a summary of this distribution
#' @export
predict.chbmod <- function(object, newdata = NA, re.form=NA, type=c("link", "response"), stat = NULL) {
  ## Compute linear predictor
  X <- model.matrix(getformula(object), newdata)
  lp <- X %*% beta_pars(object)

  ## Add RE to linear predictor from re.form
  if(!is.na(re.form)) {
    u_samp <- u_pars(object, newdata, re.form)
    zu_samp <- lapply(u_samp, rowSums)
    zu <- do.call(rbind, zu_samp)
    lp <- lp + zu
  }

  ## Compute inverse link if desired
  if(type == "link") {
    pred <- lp
  } else if(type == "response") {
    pred <- binomial()$linkinv(lp)
  }
  if(!is.null(stat))
    pred <- apply(pred, 1, stat)
  pred
}

##' Extract posterior simulations for random effects from chbmod
##' @param object object of class chbmod
##' @param newdata new data specifying random effects design
##' @param re.form \code{NA} to specify population-level predictions (i.e., setting all random effects to zero) or a character vector with values as the names of chbmod parameters and names as the corresponding identifier in newdata
##' @param type \describe{
##' \item{"link"}{conditional mean on the scale of the link function,
##' or equivalently the linear predictor of the conditional model}
##' \item{"response"}{expected value}
##' }
##' @return array of dims nrow(newdata) x nsim x nre
#' @export
u_pars <- function(object, newdata, re.form) {

  ## Get random effects ids from newdata
  ## FIXME: use unique() for speed savings
  re_mat <- newdata[,re.form, drop = F]

  ## Create string names to extract appropriate columns from sflat
  unames <- lapply(1:nrow(re_mat), function(i) paste0(names(re.form), "[", re_mat[i,], "]"))

  ## Extract RE posterior from sflat
  lapply(unames, function(p) object$sflat[,dimnames(object$sflat)[[2]] %in% p])
}
##' extract posterior distribution of fixed effects from chbmod
##' @param obj object of class chbmod
##' @return matrix (nsim x npar) of fixed effects from chbmod
##' @export
beta_pars <- function(obj) {
  t(samp(obj$sflat, "beta", partial = T))
}

##' Augment lethal temperatures from chbmod to newdata
##' @description First simplifies newdata by zeroing out std_ftemp and taking unique values.
#' @export
augment_lt <- function(object, newdata = weather_daily_default_fieldseason(), re.form=NULL, f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_names = "LT50", center = median, se_fit = NULL, xcol = 2) {

  newdata[["std_ftemp"]] <- 0
  newdata <- unique(newdata)
  lp <- predict(object, newdata = newdata, re.form = re.form, type = "link", stat = NULL)


  for(i in 1:length(p)) {
    lt_samples <- sweep(log(p[i]/(1-p[i])) - lp, 2, beta_pars(object)[xcol,],`/`)
    lt_samples <- f(lt_samples)
    lt <- apply(lt_samples, 1, center)
    newdata[[lt_names[i]]] <- lt
    if(!is.null(se_fit)) {
      newdata[[paste0(lt_names[i], "_se")]] <- apply(lt_samples, 1, se_fit)
    }
  }
  newdata
}

## working on
# predict_lt <- function(object, newdata, re.form, p, xvar) {
#   remove.var <- as.formula(paste0("~ . - ", xvar))
#   mf <- model.frame(update.formula(object$rhsform, remove.var), newdata)
#   udata <- unique(mf)
#   ## b0 (nrow(udata) x nsim)
#   b0 <- predict(object, newdata = udata, re.form = re.form, type = "link", stat = NULL)
# }

#' @export
predict_marginal_sim <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL, vc.form.marginal = NULL, type=c("link", "response"), stat = NULL, parallel = T, B = 100, method = c("within", "side")) {
  if(parallel)
    `%fdo%` <- `%dopar%`
  else
    `%fdo%` <- `%do%`
  #lt <- vector("numeric", nrow(newdata))
  samples <- foreach(i = 1:nrow(newdata), .export = c("object","newdata", "re.form.conditional", "re.form.marginal"), .packages = c("coldhardiness")) %fdo% {
    cat(paste0(i, "/", nrow(newdata), "\n"))

    lp <- predict(object, newdata = newdata, re.form = re.form.conditional, type = "link", stat = NULL)
    if(!is.na(re.form.marginal))
      u_samp <- u_pars(object, newdata[i,], re.form.marginal)[[1]] #because only doing for one row
    pred_conditional <- array(NA, dim = c(dim(lp), B))
    for(j in 1:ncol(lp)) {
      vc <- object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal]^2
      re <- rnorm(B, 0, sqrt(rowSums(vc)))
      if(!is.na(re.form.marginal)) {
        re2 <- apply(u_samp, 2, function(x) sample(x, B, replace = T))
        re2 <- rowSums(re2)
        re <- re + re2
      }
      #lp_j <- sweep(re, 1, lp[,j], `+`)
      lp_j <- t(outer(re, lp[,j], FUN = `+`))
      if(type == "link") {
        pred <- lp_j
      } else if(type == "response") {
        pred <- 1 / (1 + exp(-(lp_j)))
      }
      pred_conditional[,j,] <- pred
    }
    apply(pred_conditional, c(1,2), stat)
  }

  do.call(partial(abind, along = 3), samples)
}

#' @export
predict_marginal_sim2 <- function(object, newdata = NULL, re.form.conditional = NULL, vc.form.marginal = NULL, type=c("link", "response")) {
    lp <- predict(object, newdata = newdata, re.form = re.form.conditional, type = "link", stat = NULL)
    vc <- object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal]^2
    re <- rnorm(ncol(lp), 0, sqrt(rowSums(vc)))
    lp <- sweep(lp, 2, re, `+`)
    if(type == "link") {
      pred <- lp
    } else if(type == "response") {
      pred <- 1 / (1 + exp(-(lp)))
    }
    newdata$pred <- rowMeans(pred)
    newdata
}

#' @export
augment_marginal_lt <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL, vc.form.marginal = NULL,  f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_names = "50", center = median, se_fit = NULL, std_ftemp_seq = seq(-2, 2, by = 0.005), parallel = T, B = 100) {
  # beta = NULL,
  # if(is.null(beta))
  if(parallel)
    `%fdo%` <- `%dopar%`
  else
    `%fdo%` <- `%do%`
  newdata[["std_ftemp"]] <- 0
  newdata <- unique(newdata)
  #lt <- vector("numeric", nrow(newdata))
  lt_samples <- foreach(i = 1:nrow(newdata), .export = c("object","newdata", "re.form.conditional", "re.form.marginal", "p", "std_ftemp_seq"), .packages = c("coldhardiness")) %fdo% {
    cat(paste0(i, "/", nrow(newdata), "\n"))
    data <- newdata[rep(i, length(std_ftemp_seq)),]
    data$std_ftemp <- std_ftemp_seq

    lp <- predict(object, newdata = data, re.form = re.form.conditional, type = "link", stat = NULL)
    if(!is.na(re.form.marginal))
      u_samp <- u_pars(object, newdata[i,], re.form.marginal)[[1]] #because only doing for one row
    p_conditional <- array(NA, dim = c(dim(lp), B))
    for(j in 1:ncol(lp)) {
      #re <- matrix(rnorm(nrow(lp)*B, 0, rowSums(object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal])), ncol = B)
      re_sd <- sqrt(sum(object$sflat[i,dimnames(object$sflat)[[2]] %in% vc.form.marginal]^2))
      re <- rnorm(B, 0, re_sd)
      if(!is.na(re.form.marginal)) {
        re2 <- apply(u_samp, 2, function(x) sample(x, B, replace = T))
        # print(apply(u_samp, 2, sd))
        # re2 <- sum(colMeans(u_samp))
        # re2 <- do.call(cbind, re2)
        # re2 <- apply(object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal], 2, function(x) sample(x, B, replace = T))
        re2 <- rowSums(re2)
        re <- re + re2
      }
      #lp_j <- sweep(re, 1, lp[,j], `+`)
      lp_j <- t(outer(re, lp[,j], FUN = `+`))
      p_conditional[,j,] <- 1 / (1 + exp(-(lp_j)))
    }
    p_marginal <- apply(p_conditional, c(1,2), mean)
    sapply(p, function(p_i) apply(p_marginal, 2, function(x) std_ftemp_seq[which(x > p_i)[1]])) # S x p
  }

  lt_samples <- do.call(partial(abind, along = 3), lt_samples)
  lt_samples <- f(lt_samples)
  lt <- apply(lt_samples, c(2,3), center)
  if(!is.null(se_fit)) {
    lt_se <- apply(lt_samples, c(2,3), se_fit)
  }

  for(i in 1:length(p)) {
    newdata[[paste0("LT", lt_names[i])]] <- lt[i,]
    if(!is.null(se_fit)) {
      newdata[[paste0("SE",lt_names[i])]] <- lt_se[i,]
    }
  }
  newdata
}


#' @export
augment_marginal_lt_simple <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL, vc.form.marginal = NULL,  f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_names = "LT50", std_ftemp_seq = seq(-2, 2, by = 0.005), parallel = T) {
  # beta = NULL,
  # if(is.null(beta))
  if(parallel)
    `%fdo%` <- `%dopar%`
  else
    `%fdo%` <- `%do%`
  newdata[["std_ftemp"]] <- 0
  newdata <- unique(newdata)
  re <- rnorm(dim(object$sflat)[1], 0, sqrt(rowSums(object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal]^2)))

  #lt <- vector("numeric", nrow(newdata))
  lt_samples <- foreach(i = 1:nrow(newdata), .packages = c("coldhardiness", "splines")) %fdo% {
    cat(paste0(i, "/", nrow(newdata), "\n"))
    data <- newdata[rep(i, length(std_ftemp_seq)),]
    data$std_ftemp <- std_ftemp_seq

    lp <- predict(object, newdata = data, re.form = re.form.conditional, type = "link", stat = NULL)
    if(!is.na(re.form.marginal))
      u_samp <- u_pars(object, newdata[i,], re.form.marginal)[[1]] #because only doing for one row

    lp <- sweep(lp, 2, re, `+`)
    p_conditional <- 1 / (1 + exp(-(lp)))

    p_marginal <- rowMeans(p_conditional)#apply(p_conditional, c(1,2), mean)
    sapply(p, function(p_i) std_ftemp_seq[which(p_marginal > p_i)[1]]) # S x p
  }

  lt_samples <- do.call(rbind, lt_samples)
  lt <- f(lt_samples)

  for(i in 1:length(p)) {
    newdata[[paste0("LT", lt_names[i])]] <- lt[,i]
  }
  newdata
}



#
# augment_marginal_lt <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL, vc.form.marginal = NULL,  f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_names = "LT50", center = median, se_fit = NULL, std_ftemp_seq = seq(-2, 2, by = 0.005), parallel = T, B = 100) {
#   # newdata[["std_ftemp"]] <- 0
#   newdata <- unique(newdata)
#   lt <- lt_se <- list()
#   for(i in 1:nrow(newdata)) {
#     cat(paste0(i, "/", nrow(newdata)))
#     data <- newdata[rep(i, length(std_ftemp_seq)),]
#     data$std_ftemp <- std_ftemp_seq
#     # data <- newdata[rep(1:nrow(newdata), each = length(std_ftemp_seq)),]
#     # data$std_ftemp <- rep(std_ftemp_seq, nrow(newdata))
#     samples <- predict_marginal_sim(object, newdata = data, re.form.conditional = re.form.conditional, re.form.marginal = re.form.marginal, vc.form.marginal = vc.form.marginal, type="response")
#     lt_samples <- lapply(samples, function(p_marginal) sapply(p, function(p_i) apply(p_marginal, 2, function(x) std_ftemp_seq[which(x > p_i)[1]])))
#     lt_samples <- do.call(partial(abind, along = 3), lt_samples)
#     lt_samples <- f(lt_samples)
#     lt[[i]] <- apply(lt_samples, c(2,3), center)
#     if(!is.null(se_fit)) {
#       lt_se[[i]] <- apply(lt_samples, c(2,3), se_fit)
#     }
#   }
#   browser()
#
#
#   for(i in 1:length(p)) {
#     newdata[[paste0("LT", lt_names[i])]] <- lt[i,]
#     if(!is.null(se_fit)) {
#       newdata[[paste0("SE",lt_names[i])]] <- lt_se[i,]
#     }
#   }
#   newdata
# }

#' @export
flatten_stan_array <- function(x) {
  lst <- lapply(1:dim(x)[2], function(i) x[,i,, drop = T])
  do.call(rbind, lst)
}

#' @export
flatten_stanfit <- function(obj) {
  if(class(obj) == "stanfit") {
    flatten_stan_array(as.array(obj))
  } else if(length(dim(obj)) == 3) {
    flatten_stan_array(obj)
  } else {
    obj
  }
}

