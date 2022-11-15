#' @export
plot.cherry_s1 <- function(obj, which = 1:6) {
  show <- rep(FALSE, 6)
  show[which] <- TRUE
  if (show[1L]) {  # Plot the probability estimates for each date vs. ftemp
    props <- dataset(obj) %>%
      group_by(seasonXfield_ID, fieldXdate_ID, Field, date, ftemp, std_ftemp) %>%
      summarise(prop = sum(NoFlowersLive)/sum(NoFlowers))
    props$pred <- predict(model_bayes, newdata = props, re.form = c(r_seasonXfield = "seasonXfield_ID", r_fieldXdate = "fieldXdate_ID"), type = "response", stat = median)
    lts <- augment_lt(obj, newdata = props %>% ungroup() %>% select(seasonXfield_ID, fieldXdate_ID, Field, date), re.form = c(r_seasonXfield = "seasonXfield_ID", r_fieldXdate = "fieldXdate_ID"),
                      f = unstd_ftemp, p = 0.5, lt_names = "LT50", center = median)
    plt <- props %>% ggplot(aes(ftemp, prop)) + geom_point() + geom_line(aes(y = pred)) + geom_vline(data = lts, aes(xintercept = LT50))+ facet_wrap(Field~date)
    print(plt)
  }
  if (show[2L]) {
    print(traceplot(obj$stanfit, obj$par_grps$beta))
  }
  if (show[3L]) {
    print(traceplot(obj$stanfit, obj$par_grps$sd))
  }
  if (show[4L]) {
    print(traceplot(obj$stanfit, obj$par_grps$seasonXfield))
  }
  if (show[5L]) {
    print(traceplot(obj$stanfit, obj$par_grps$fieldXdate))
  }
  if (show[6L]) {
    print(traceplot(obj$stanfit, obj$par_grps$fieldXdate_fi))
  }

}
