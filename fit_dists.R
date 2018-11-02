library(rstan)
par (mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)



load("sumscores_moba.Rdata")
stan_normal = readRDS("normal.rds")
stan_binom = readRDS("binom.rds")
stan_betabinom = readRDS("betabinom.rds")


md = md[,c(3,1,4)]
par(mfrow = c(1,3))
for (v in names(md)) hist(md[,v], main = v, xlab = "")

pd = vector(length = 3, mode = "list")
names(pd) = names(md)

par(mfrow = c(1,3))
for (v in names(md)) {
  y = md[,v]
  y = y[!is.na(y)]
  range_y = range(y)
  max_y = max(y)
  y = y[sample(length(y),25000)]
  
  sf_normal = sampling(stan_normal, data = list(N = length(y), y = y), chains = 1, iter = 500)
  y_rep_normal = extract(sf_normal,"y_rep")$y_rep
  
  sf_binom = sampling(stan_binom, data = list(N = length(y), J = max_y, y = y), chains = 1, iter = 500)
  y_rep_binom = extract(sf_binom,"y_rep")$y_rep
  
  sf_betabinom = sampling(stan_betabinom, data = list(N = length(y), J = max_y, y = y), chains = 1, iter = 500)
  y_rep_betabinom = extract(sf_betabinom,"y_rep")$y_rep
  
  tmp = rbind(y_rep_betabinom,y_rep_binom,y_rep_normal, y)
  breaks = seq(round(min(tmp))-.5,round(max(tmp))+.5)
  
  counts_normal = c()
  for (i in 1:250) counts_normal = rbind(counts_normal,
                                        c(0,hist(y_rep_normal[i,],
                                                 breaks = breaks,
                                                 plot = F)$counts))
  
  counts_binom = c()
  for (i in 1:250) counts_binom = rbind(counts_binom,
                                       c(0,hist(y_rep_binom[i,],
                                                breaks = breaks,
                                                plot = F)$counts))

  counts_betabinom = c()
  for (i in 1:250) counts_betabinom = rbind(counts_betabinom,
                                       c(0,hist(y_rep_betabinom[i,],
                                                breaks = breaks,
                                                plot = F)$counts))
  
  counts_data = c(0,hist(y,
                     breaks = breaks,
                     plot = F)$counts)
  
  
  pd[[v]] = list(counts_normal = counts_normal,
                 counts_binom = counts_binom,
                 counts_betabinom = counts_betabinom,
                 counts_data = counts_data,
                 resid_normal = apply(y_rep_normal, 1, function(x) x -y),
                 resid_binom = apply(y_rep_binom, 1, function(x) x -y),
                 resid_betabinom = apply(y_rep_betabinom, 1, function(x) x -y),
                 max_y = max_y,
                 range_y = range_y,
                 breaks = breaks)
  
  # ylim = range(rbind(pd[[v]]$counts_binom,
  #                    pd[[v]]$counts_betabinom,
  #                    pd[[v]]$counts_normal,
  #                    pd[[v]]$counts_data))
  # xlim = quantile(pd[[v]]$y_rep_normal,c(.0001,1))
  # xlim[2] = ifelse(xlim[2]<pd[[v]]$max_y,pd[[v]]$max_y,xlim[2])
  # par(mar = c(2,.1,2,.1))
  # plot(0,type = "n", xlim = xlim, ylim = ylim, 
  #      main = v, yaxt = "n", 
  #      ylab  = "", bty = "n")
  # abline(v = range_y, lty = 2, col = "grey")
  # lines(breaks,counts_data, type='S', lwd = 3)
  # for(i in 1:100) {
  #   lines(breaks,counts_normal[i,], type='S', col = adjustcolor("red", alpha = .075))
  #   lines(breaks,counts_binom[i,], type='S', col = adjustcolor("blue", alpha = .075))
  #   lines(breaks,counts_betabinom[i,], type='S', col = adjustcolor("green", alpha = .075))
  # }
  # legend("topleft",
  #        col = c("red","blue","green"),
  #        lty = 1, lwd = 2,
  #        legend = c("normal","binomial","beta-binomial"),
  #        bty = "n")
}

save(pd, file = "pd_dists.Rdata")
