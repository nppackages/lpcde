#Loading dataset
data("iris")
attach(iris)

#loading conditional density estimation packages
library("np")
library("hdrcde")
library(latex2exp)

#selecting variables
y = as.matrix(iris[, 1], ncol=1)
x = as.matrix(iris[, 3], ncol=1)
x_evals = quantile(x)

#scatter plot of data with chosen conditioning values
plot(x,y, xlab = "Petal length", ylab = "Sepal length")
abline(v=x_evals[2], col="green")
abline(v=x_evals[3], col="green")
abline(v=x_evals[4], col="green")

#estimation at first conditional value
lpcde_est_q1 = lpcde::lpcde(x_data=x, y_data=y, x=x_evals[2], bw=0.9, ng=30, nonneg = TRUE, normalize = TRUE)
hdr_est_q1 = cde(x, y, x.margin=x_evals[2], y.margin = lpcde_est_q1$Estimate[,1])
npbws_q1 = npcdensbw(x, y)
np_est_q1 = fitted(npcdens(exdat=rep(x_evals[2], 30), eydat=lpcde_est_q1$Estimate[,1], bws=npbws_q1))

#estimation at second conditional value
lpcde_est_q2 = lpcde::lpcde(x_data=x, y_data=y, x=x_evals[3], ng=30, bw=0.9, nonneg = TRUE, normalize = TRUE)
hdr_est_q2 = cde(x, y, x.margin=x_evals[3], y.margin = lpcde_est_q2$Estimate[,1])
npbws_q2 = npcdensbw(x, y)
np_est_q2 = fitted(npcdens(exdat=rep(x_evals[3], 30), eydat=lpcde_est_q2$Estimate[,1], bws=npbws_q2))

#estimation at third conditional value
lpcde_est_q3 = lpcde::lpcde(x_data=x, y_data=y, x=x_evals[4], ng=30, bw=1, nonneg = TRUE, normalize = TRUE)
hdr_est_q3 = cde(x, y, x.margin=x_evals[4], y.margin = lpcde_est_q3$Estimate[,1])
npbws_q3 = npcdensbw(x, y)
np_est_q3 = fitted(npcdens(exdat=rep(x_evals[4], 30), eydat=lpcde_est_q3$Estimate[,1], bws=npbws_q3))

#plotting all CDEs
plot(1, main=TeX(r'($f(y|x=1.6)$)'),  xlab="Sepal length", ylab="density", ylim=c(-0.8, 1.6), xlim=c(4.8, 6.9), cex=1.5)
lines(lpcde_est_q1$Estimate[, 1], lpcde_est_q1$Estimate[,3])
lines(hdr_est_q1$y, hdr_est_q1$z, col=2)
lines(lpcde_est_q1$Estimate[, 1], np_est_q1, col=3)
legend('topright',lwd=1, legend=c('lpcde', 'hdrcde', 'np'), col=c(1,2,3))

plot(1, main=TeX(r'($f(y|x=4.35)$)'),  xlab="Sepal length", ylab="density", ylim=c(-0.8, 1.6), xlim=c(4.8, 6.9))
lines(lpcde_est_q2$Estimate[, 1], lpcde_est_q2$Estimate[,3])
lines(hdr_est_q2$y, hdr_est_q2$z, col=2)
lines(lpcde_est_q2$Estimate[, 1], np_est_q2, col=3)
legend('topright',lwd=1, legend=c('lpcde', 'hdrcde', 'np'), col=c(1,2,3))

plot(1, main=TeX(r'($f(y|x=5.1)$)'),  xlab="Sepal length", ylab="density", ylim=c(-0.8, 1.6), xlim=c(4.8, 6.9))
lines(lpcde_est_q3$Estimate[, 1], lpcde_est_q3$Estimate[,3])
lines(hdr_est_q3$y, hdr_est_q3$z, col=2)
lines(lpcde_est_q3$Estimate[, 1], np_est_q3, col=3)
legend('topright',lwd=1, legend=c('lpcde', 'hdrcde', 'np'), col=c(1,2,3))

# lpcde estimates with CI
plot(lpcde_est_q1, xlabel="Sepal length", ylabel="density", title=TeX(r'($f(y|x=1.6)$ with confidence bands)')) + ggplot2::ylim(-0.8, 1.6) + ggplot2::xlim(4.8, 6.9) + ggplot2::theme(legend.position='none', panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), text = ggplot2::element_text(size = 7))
plot(lpcde_est_q2,xlabel="Sepal length", ylabel="density", title=TeX(r'($f(y|x=4.35)$ with confidence bands)'))+ ggplot2::ylim(-0.8, 1.6) + ggplot2::xlim(4.8, 6.9) + ggplot2::theme(legend.position='none', panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                                                                                                                                                                                     panel.background = ggplot2::element_blank(), text = ggplot2::element_text(size = 7))
plot(lpcde_est_q3, xlabel="Sepal length", ylabel="density", title=TeX(r'($f(y|x=5.1)$ with confidence bands)'))+ ggplot2::ylim(-0.8, 1.6) + ggplot2::xlim(4.8, 6.9) + ggplot2::theme(legend.position='none', panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                                                                                                                                                                                     panel.background = ggplot2::element_blank(), text = ggplot2::element_text(size = 7))
