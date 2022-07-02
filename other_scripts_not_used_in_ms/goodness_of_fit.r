lipid_gsets = get_gset_by_category(lipid_keys, all)
unique_gsets = get_unique_gsets(lipid_gsets, expressed)
gsva_out = get_gsva_per_celltype(expressed, unique_gsets, av_expression)

# show for cholesterol biosynthesis and ATF

temp = (out_lipid$sig_hits)
temp = temp[temp$celltype == 'Oli',]

fit <- lmFit(t(gsva_out[['Oli']]), design=mod)
fit <- eBayes(fit)
coefficients = fit$coefficients
predictors = mod

r_squared = list()
err = list()
sigma2 = list()
yhats = list()
actual = list()
observed = list()
for(i in temp$path){
    yhat = predictors %*% coefficients[i,]
    y = gsva_out$Oli[,i]

    ess = sum((yhat - mean(y))**2)
    rss = sum((y - yhat)**2)
    tss = ess + rss
    r2 = as.data.frame(ess/tss)
    r2$path = i

    r_squared[[i]] = r2


    error = as.data.frame(y-yhat)
    error$path = i
    err[[i]] = error

    sigma_2 = (1/length(error[,1]))*(t(error[,1])%*%error[,1])
    sigma2[[i]] = sigma_2

    yhat = as.data.frame(yhat)
    yhat$path = i
    #se = sqrt(sigma_2*((predictors%*%(solve(t(predictors)%*%predictors)))%*%t(predictors)))
    CI = qt(p=.975, df=nrow(yhat)-ncol(mod)-1, lower.tail=FALSE)*sd(yhat[,1])/sqrt(length(yhat[,1]))

    yhats[[i]] = yhat



    actuals = as.data.frame(unlist(lapply(seq(0, 1, 0.01), function(i) qnorm(i,mean=0,sd=sqrt(sigma_2)))))
    actuals$path = i

    actual[[i]] = actuals
    obser = as.data.frame(quantile(yhat[,1], probs = seq(0, 1, 0.01)))
    obser$path = i
    obser$CI_upper = obser[,1] + CI
    obser$CI_lower = obser[,1] - CI
    observed[[i]] = obser
}

# plot predictions vs residuals
yh = do.call('rbind',yhats)
er = (do.call('rbind',err))

cd = cbind(yh,er)
colnames(cd) = c('yhat', 'path','error', 'path2')
ncol(mod)

options(repr.plot.width = 20, repr.plot.height =3.5, repr.plot.res = 100)

ggplot(cd, aes(x=yhat, y=error, shape=path, color=path)) +
  geom_point() + theme_classic() +  facet_grid(. ~ path) + geom_abline(intercept = 0, slope = 0, linetype = 2)

# show QQ plots
ac = do.call('rbind',actual)
obs = (do.call('rbind',observed))
dd = cbind(ac, obs)
colnames(dd) = c('actual', 'path','observed', 'path2','CI_upper','CI_lower')
rownames(dd) = NULL
dd$path2 = NULL
base = cbind(ac[,1], 'baseline',ac[,1], NA, NA)
colnames(base) = c('actual', 'path', 'observed', 'CI_upper','CI_lower')
options(repr.plot.width = 20, repr.plot.height =3.5, repr.plot.res = 100)

ggplot(dd, aes(x=actual, y=observed, shape=path, color=path)) +
  geom_point()+ geom_abline(intercept = 0, slope = 1, linetype = 2) +#+ ylim(-.4,.4)+ xlim(-.4,.4)
    geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper, color = path), linetype=1, alpha=0.1) +  facet_grid(. ~ path) + theme_classic()



# not sure about these confidence intervals.. see if computing variance correctly
