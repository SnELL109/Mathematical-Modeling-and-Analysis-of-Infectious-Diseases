library(deSolve)
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdaL = betaL * (theta1 * I1 + theta2 * I2 + theta3 * I3)/(SL+SH+PO+PR+I1+I2+I3)
    lambdaH = betaH * (theta1 * I1 + theta2 * I2 + theta3 * I3)/(SL+SH+PO+PR+I1+I2+I3)
    lambdaR = betaR * (theta1 * I1 + theta2 * I2 + theta3 * I3)/(SL+SH+PO+PR+I1+I2+I3)
    dSL <- Lambda * (SL+SH+PO+PR+I1+I2+I3) + WHL * SH + phiL * PO - (lambdaL + WLH + mu) * SL
    dSH <- WHL * SL + phiH * PO + phiR * PR - (lambdaH + WHL + tauR + mu) * SH
    dPO <- lambdaL * PL * SL + lambdaH * PH * SH - (sigmaO1 + phiL + phiH + WOR + mu) * PO
    dPR <- tauR * SH + WOR * PO - (lambdaR + phiR + mu) * PR
    dI1 <- lambdaL * (1 - PL) * SL + lambdaH * (1 - PH) * SH + sigmaO1 * PO + lambdaR *PR - (sigma12 + mu) * I1 
    dI2 <- sigma12 * I1 - (sigma23 + mu) * I2
    dI3 <- sigma23 * I2 - (mud + mu) * I3
    return(list(c(dSL, dSH, dPO, dPR, dI1, dI2, dI3)))
  })
}

init <- c(SL = 0.999, SH = 0.0009, PO = 0, PR = 0, I1 = 0.000001, I2 = 0.000059, I3 = 0.00004)
parameters <- c(Lambda = 0.00817, mu = 0.00747, mud = 0.0909, WLH = 0.06, WHL = 0.04, WOR = 0.001, 
                sigmaO1 = 0.023, sigma12 = 8.67, sigma23 = 0.05, PL = 0.5, PH = 0.5, phiL = 0.05, phiH = 1/24, tauR = 0.01, phiR = 0.005,
                betaL = 0.25, betaH = 1.25, betaR = 0.125, theta1 = 1, theta2 = 0.43, theta3 = 1.5)
times <- seq(0, 70, by = 0.1)
out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
out$time <- NULL
head(out, 10)

matplot(times, out, type = "l", xlab = "Time (years)", ylab = "Compartments", main = "HIV Transmission Model (PEP + PrEP)", lwd = 1, lty = 1, bty = "l", col = 2:8)
legend(40, 0.7, c(expression(S[L]), expression(S[H]), expression(P[O]), expression(P[R]), expression(I[1]), expression(I[2]), expression(I[3])), pch = 1, col = 2:8)
