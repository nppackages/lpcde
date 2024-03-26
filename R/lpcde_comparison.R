data("iris")
attach(iris)
library("np")
library("hdrcde")
library(latex2exp)
y = as.matrix(iris[, 1], ncol=1)
x = as.matrix(iris[, 3], ncol=1)
x_evals = quantile(x)

plot(x,y)
abline(v=x_evals[2], col="green")
abline(v=x_evals[3], col="green")
abline(v=x_evals[4], col="green")


lpcde_est_q1 = lpcde::lpcde(x_data=x, y_data=y, x=x_evals[2], bw=0.9, ng=30, nonneg = TRUE, normalize = TRUE)
hdr_est_q1 = cde(x, y, x.margin=x_evals[2], y.margin = lpcde_est_q1$Estimate[,1])
npbws_q1 = npcdensbw(x, y)
np_est_q1 = fitted(npcdens(exdat=rep(x_evals[2], 30), eydat=lpcde_est_q1$Estimate[,1], bws=npbws_q1))

lpcde_est_q2 = lpcde::lpcde(x_data=x, y_data=y, x=x_evals[3], ng=30, bw=0.9, nonneg = TRUE, normalize = TRUE)
hdr_est_q2 = cde(x, y, x.margin=x_evals[3], y.margin = lpcde_est_q2$Estimate[,1])
npbws_q2 = npcdensbw(x, y)
np_est_q2 = fitted(npcdens(exdat=rep(x_evals[3], 30), eydat=lpcde_est_q2$Estimate[,1], bws=npbws_q2))

lpcde_est_q3 = lpcde::lpcde(x_data=x, y_data=y, x=x_evals[4], ng=30, bw=1, nonneg = TRUE, normalize = TRUE)
hdr_est_q3 = cde(x, y, x.margin=x_evals[4], y.margin = lpcde_est_q3$Estimate[,1])
npbws_q3 = npcdensbw(x, y)
np_est_q3 = fitted(npcdens(exdat=rep(x_evals[4], 30), eydat=lpcde_est_q3$Estimate[,1], bws=npbws_q3))

plot(1, main=TeX(r'($f(y|x=1.6)$)'),  xlab="Sepal length", ylab="density", ylim=c(-0.8, 1.6), xlim=c(4.8, 6.9))
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

# lpcde estimate with CI
plot(lpcde_est_q1, xlabel="Sepal length", ylabel="density", title=TeX(r'($f(y|x=1.6)$ with confidence bands)')) + ggplot2::ylim(-0.8, 1.6) + ggplot2::xlim(4.8, 6.9) + ggplot2::theme(legend.position='none')
plot(lpcde_est_q2,xlabel="Sepal length", ylabel="density", title=TeX(r'($f(y|x=4.35)$ with confidence bands)'))+ ggplot2::ylim(-0.8, 1.6) + ggplot2::xlim(4.8, 6.9) + ggplot2::theme(legend.position='none')
plot(lpcde_est_q3, xlabel="Sepal length", ylabel="density", title=TeX(r'($f(y|x=5.1)$ with confidence bands)'))+ ggplot2::ylim(-0.8, 1.6) + ggplot2::xlim(4.8, 6.9) + ggplot2::theme(legend.position='none')





#open3d()
#
#lines3d(x_evals[2], lpcde_est_q1$Estimate[,1], lpcde_est_q1$Estimate[,3])
#lines3d(x_evals[3], lpcde_est_q2$Estimate[,1], lpcde_est_q2$Estimate[,3])
#lines3d(x_evals[4], lpcde_est_q3$Estimate[,1], lpcde_est_q3$Estimate[,3])
#
#lines3d(x_evals[2], hdr_est_q1$y, hdr_est_q1$z, color='red')
#lines3d(x_evals[3], hdr_est_q2$y, hdr_est_q2$z, color='red')
#lines3d(x_evals[4], hdr_est_q3$y, hdr_est_q3$z, color='red')
#
#lines3d(x_evals[2], lpcde_est_q1$Estimate[,1], np_est_q1, color='blue')
#lines3d(x_evals[3], lpcde_est_q2$Estimate[,1], np_est_q2, color='blue')
#lines3d(x_evals[4], lpcde_est_q3$Estimate[,1], np_est_q3, color='blue')
#
#axes3d()
#title3d('', '', 'Petal Length', 'Sepal Length', 'density')
#legend3d("topright", c("lpcde", "hdrcde", "np"), pch = c(1,25, 16), col = c('black', 'red', 'blue'), cex=1)
