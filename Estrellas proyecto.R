# FUNCIONES GENERALES--------------------------------------------

# Funcion para graficar los contornos de verosilitud relativa de (mu, sigma)
plotRelative<-function(l, aG, bG, xL, xR, yL, yR, levels, n, xlab = expression(paste(beta[0])),
                       ylab = expression(paste(beta[1])), main="Contornos de verosimilitud relativa"){
  x_vec = seq(from = xL, to = xR, length.out = n)
  y_vec = seq(from = yL, to = yR, length.out = n)
  R<-function(a, b){
    return(exp(l(a, b) - l(aG, bG)))
  }
  Rmat = matrix(nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      Rmat[i, j] = R(x_vec[i], y_vec[j])
    }
  }
  contour(x_vec,y_vec,Rmat,level=levels,xlab=xlab,ylab=ylab,main=main)
}

# Funcion para grafica PP
pp_plot = function(X, ag, bg, confidence){
  Y = pnorm(sort(X), ag, bg)
  X = seq(1/(n+1), 1, length.out = n)
  # puntos de la muestra
  plot(X, Y,
       main = "Gráfica PP",
       xlab = "Probabilidades empíricas",
       ylab = "Probabilidades teóricas",
       pch = 19,
       cex = 0.5)
  # identidad
  abline(a = 0, b = 1, col = "red", lwd = 2)
  # bandas de confianza
  points(X, qbeta((1 - confidence)/2, 1:n, n + 1 - 1:n),
         type = "l",
         lty = 2)
  points(X, qbeta((1 + confidence)/2, 1:n, n + 1 - 1:n),
         type = "l",
         lty = 2)
}


Y <- c(5778,3042,5810,5260,3134,2800,3526,3240,3058,5073,3700,2966,4450,4120,3395,4620,2840,5380,3105,3380,2904,3570,3599,3310,5942,5181)
X <- c(1,0.0017,1.529,0.6,0.00346,0.00002, 0.02, 0.0042, 0.00011, 0.31, 0.012, 0.00029, 0.153, 0.095, 0.0013, 0.22, 0.00065, 0.75, 0.00085, 0.091, 0.0073, 0.004, 0.06, 0.0034,1.08, 0.68)

n <- length(X)
n1 <-length(Y)


# Vector de 1's
Aux <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

# Matriz de diseno
XX <- cbind(Aux, X)

# Matriz transpuesta
XXt <- rbind(Aux, X)

# Estimadores de maxima verosimilitud
b1g <- (n * sum(X*Y) - sum(X)*sum(Y)) / det(XXt %*% XX)
b0g <- (sum(Y) / n) - b1g * (sum(X) / n)
bg <- rbind(b0g, b1g)
sigmag <- (t(Y - XX %*% bg) %*% (Y - XX %*% bg)) / n

# Residuos observados
R <- Y - b0g - (b1g * X)

# Parametros asociados a los residuos con distribucion Normal
t1 <- sum(R)
t2 <- sum(R^2)
mg <- t1 / n
H <- t2 - (t1^2 / n)
sg <- sqrt(H / n)


t2aux = sum((R/sqrt(sg))^2)
t1aux = sum(R)

# FUNCIONES EMVR's, VEROSIMILITUD, LOGVEROSIMILITUD (+RELATIVAS) -------------------------------

# Emvr de beta1
emvrb1 <- function(beta0){
  return(((sum(X * Y) - beta0 * sum(X)) / (sum(X^2))))
}

# Emvr de beta0
emvrb0 <- function(beta1){
  return(((sum(Y) - beta1 * sum(X)) / n))
}

# Inversa del emvr de beta0
Iemvrb0 <- function(beta1){
  return((sum(Y) - n * beta1) / sum(X))
}

# Logverosimilitud perfil de beta_0
lpb0 <- function(beta0){
  return((-n/2)*log((sum((Y - beta0 - emvrb1(beta0)*X)^2)) / n) - n/2)
}

# Logverosimilitud perfil de beta_1
lpb1 <- function(beta1){
  return((-n/2)*log((sum((Y - emvrb0(beta1) - beta1*X)^2)) / n) - n/2)
}

# Variables auxiliares para la logverosimilitud perfil de sigma
a1 <- t(Y - (XX %*% bg))
a2 <- Y - (XX %*% bg)

# Logverosimilitud perfil de sigma
lps <- function(sigma){
  return((-n/2)*log(sigma) - ((1/(2*sigma))*(a1 %*% a2)))
}

# Verosimilitud perfil relativa de beta_0
Rpb0 <- function(beta0){
  return(exp(lpb0(beta0) - lpb0(b0g)))
}

# Verosimilitud perfil relativa de beta_1
Rpb1 <- function(beta1){
  return(exp(lpb1(beta1) - lpb1(b1g)))
}

# Verosimilitud perfil relativa de sigma
Rps <- function(sigma){
  return(exp(lps(sigma) - lps(sigmag)))
}

# Logverosimilitud de (beta0, beta1)
Blogver <- function(beta0, beta1){
  return((-n/2) * log((t(Y - XX %*% c(beta0, beta1)) %*% (Y - XX %*% c(beta0, beta1))) / n) - n/2)
}

# Logverosimilitud relativa de (beta0, beta1)
rlogver <- function(beta0, beta1){
  return(Blogver(beta0, beta1) - Blogver(b0g, b1g))
}

# FUNCIONES EMVR's, VEROSIMILITUD, LOGVEROSIMILITUD (+RELATIVAS) -------------------------------

# Emvr de beta1
emvrb1 <- function(beta0){
  return(((sum(X * Y) - beta0 * sum(X)) / (sum(X^2))))
}

# Emvr de beta0
emvrb0 <- function(beta1){
  return(((sum(Y) - beta1 * sum(X)) / n))
}

# Inversa del emvr de beta0
Iemvrb0 <- function(beta1){
  return((sum(Y) - n * beta1) / sum(X))
}

# Logverosimilitud perfil de beta_0
lpb0 <- function(beta0){
  return((-n/2)*log((sum((Y - beta0 - emvrb1(beta0)*X)^2)) / n) - n/2)
}

# Logverosimilitud perfil de beta_1
lpb1 <- function(beta1){
  return((-n/2)*log((sum((Y - emvrb0(beta1) - beta1*X)^2)) / n) - n/2)
}

# Variables auxiliares para la logverosimilitud perfil de sigma
a1 <- t(Y - (XX %*% bg))
a2 <- Y - (XX %*% bg)

# Logverosimilitud perfil de sigma
lps <- function(sigma){
  return((-n/2)*log(sigma) - ((1/(2*sigma))*(a1 %*% a2)))
}

# Verosimilitud perfil relativa de beta_0
Rpb0 <- function(beta0){
  return(exp(lpb0(beta0) - lpb0(b0g)))
}

# Verosimilitud perfil relativa de beta_1
Rpb1 <- function(beta1){
  return(exp(lpb1(beta1) - lpb1(b1g)))
}

# Verosimilitud perfil relativa de sigma
Rps <- function(sigma){
  return(exp(lps(sigma) - lps(sigmag)))
}

# Logverosimilitud de (beta0, beta1)
Blogver <- function(beta0, beta1){
  return((-n/2) * log((t(Y - XX %*% c(beta0, beta1)) %*% (Y - XX %*% c(beta0, beta1))) / n) - n/2)
}

# Logverosimilitud relativa de (beta0, beta1)
rlogver <- function(beta0, beta1){
  return(Blogver(beta0, beta1) - Blogver(b0g, b1g))
}


# INTERVALOS DE CONFIANZA SIGMA CUADRADA (VARIANZA)------------------------------------------------------------------------------

# Funcion para obtener el la probabilidad de un nivel 
ints <- function(x){
  return(log(Rps(x)) - log(0.0165))
}

i1 = uniroot(ints, c(2, 5))$root
i2 = uniroot(ints, c(10, 20))$root

p <- pchisq((n * sigmag)/ i1, 18) - pchisq((n * sigmag) / i2, 18)

# Tras probar con distintos valores de de c y ver la probabilidad correspondiente
# en a función anterior, se llegó que las c's son 0.1870, 0.1500, 0.0940 y 0.016.

# Extremos del intervalo de verosimilitud de nivel c = 0.1870
ints90 <- function(sigma){
  return(log(Rps(sigma)) - log(0.187))
}
uno90 = uniroot(ints90, c(0.0001, sigmag))$root
dos90 = uniroot(ints90, c(sigmag, 15*sigmag))$root

# Extremos del intervalo de verosimilitud de nivel c = 0.1500
ints92 <- function(sigma){
  return(log(Rps(sigma)) - log(0.150))
}
uno92 = uniroot(ints92, c(0.0001, sigmag))$root
dos92 = uniroot(ints92, c(sigmag, 15*sigmag))$root

# Extremos del intervalo de verosimilitud de nivel c = 0.0940
ints95 <- function(sigma){
  return(log(Rps(sigma)) - log(0.094))
}
uno95 = uniroot(ints95, c(0.0001, sigmag))$root
dos95 = uniroot(ints95, c(sigmag, 15*sigmag))$root

# Extremos del intervalo de verosimilitud de nivel c = 0.0165
ints99 <- function(sigma){
  return(log(Rps(sigma)) - log(0.0165))
}
uno99 = uniroot(ints99, c(0.0001, sigmag))$root
dos99 = uniroot(ints99, c(sigmag, 15*sigmag))$root

# INTERVALOS DE CONFIANZA BETA 0 -------------------------------------

# Extremos del intervalo de verosimilitud de nivel c
extgb0 <- function(c){
  ints <- function(beta){
    return(log(Rpb0(beta)) - log(c))
  }
  
  uno = uniroot(ints, c(b0g / 10, b0g))$root
  dos = uniroot(ints, c(b0g, b0g * 10))$root
  
  int <- c(uno, dos)
  
  return(int)
}

# Funcion para obtener c para 90%
gb090 <- function(c){
  ints <- function(beta){
    return(log(Rpb0(beta)) - log(c))
  }
  
  uno = uniroot(ints, c(b0g / 10, b0g))$root
  dos = uniroot(ints, c(b0g, b0g * 10))$root
  
  piv1 <- sqrt((n-2)/n) * (b0g - uno) * sqrt(det(XXt %*% XX) / (sigmag * sum(X^2)))
  piv2 <- sqrt((n-2)/n) * (b0g - dos) * sqrt(det(XXt %*% XX) / (sigmag * sum(X^2)))
  
  return(pt(piv1, n-2, lower.tail = TRUE) - pt(piv2, n-2, lower.tail = TRUE) - 0.90)
}

# Funcion para obtener c para 92%
gb092 <- function(c){
  ints <- function(beta){
    return(log(Rpb0(beta)) - log(c))
  }
  
  uno = uniroot(ints, c(b0g / 10, b0g))$root
  dos = uniroot(ints, c(b0g, b0g * 10))$root
  
  piv1 <- sqrt((n-2)/n) * (b0g - uno) * sqrt(det(XXt %*% XX) / (sigmag * sum(X^2)))
  piv2 <- sqrt((n-2)/n) * (b0g - dos) * sqrt(det(XXt %*% XX) / (sigmag * sum(X^2)))
  
  return(pt(piv1, n-2, lower.tail = TRUE) - pt(piv2, n-2, lower.tail = TRUE) - 0.92)
}

# Funcion para obtener c para 95%
gb095 <- function(c){
  ints <- function(beta){
    return(log(Rpb0(beta)) - log(c))
  }
  
  uno = uniroot(ints, c(b0g / 10, b0g))$root
  dos = uniroot(ints, c(b0g, b0g * 10))$root
  
  piv1 <- sqrt((n-2)/n) * (b0g - uno) * sqrt(det(XXt %*% XX) / (sigmag * sum(X^2)))
  piv2 <- sqrt((n-2)/n) * (b0g - dos) * sqrt(det(XXt %*% XX) / (sigmag * sum(X^2)))
  
  return(pt(piv1, n-2, lower.tail = TRUE) - pt(piv2, n-2, lower.tail = TRUE) - 0.95)
}

# Funcion para obtener c para 99%
gb099 <- function(c){
  ints <- function(beta){
    return(log(Rpb0(beta)) - log(c))
  }
  
  uno = uniroot(ints, c(b0g / 10, b0g))$root
  dos = uniroot(ints, c(b0g, b0g * 10))$root
  
  piv1 <- sqrt((n-2)/n) * (b0g - uno) * sqrt(det(XXt %*% XX) / (sigmag * sum(X^2)))
  piv2 <- sqrt((n-2)/n) * (b0g - dos) * sqrt(det(XXt %*% XX) / (sigmag * sum(X^2)))
  
  return(pt(piv1, n-2, lower.tail = TRUE) - pt(piv2, n-2, lower.tail = TRUE) - 0.99)
}

# Obtenemos las c's
c090 = uniroot(gb090, c(0.0001, 0.8))$root
c092 = uniroot(gb092, c(0.0001, 0.8))$root
c095 = uniroot(gb095, c(0.0001, 0.8))$root
c099 = uniroot(gb099, c(0.0001, 0.8))$root

# Obtenemos los intervalos correspondientes
ext090 = extgb0(c090)
ext092 = extgb0(c092)
ext095 = extgb0(c095)
ext099 = extgb0(c099)

# INTERVALOS DE CONFIANZA BETA 1 ------------------------------------

# Extremos del intervalo de verosimilitud de nivel c
extgb1 <- function(c){
  ints <- function(beta){
    return(log(Rpb1(beta)) - log(c))
  }
  
  uno = uniroot(ints, c(b1g / 10, b1g))$root
  dos = uniroot(ints, c(b1g, b1g * 10))$root
  
  int <- c(uno, dos)
  
  return(int)
}

# Funcion para obtener c para 90%
gb190 <- function(c){
  ints <- function(beta){
    return(log(Rpb1(beta)) - log(c))
  }
  
  uno = uniroot(ints, c(b1g / 10, b1g))$root
  dos = uniroot(ints, c(b1g, b1g * 10))$root
  
  piv1 <- sqrt((n-2)/n) * (b1g - uno) * sqrt(det(XXt %*% XX) / (n * sigmag))
  piv2 <- sqrt((n-2)/n) * (b1g - dos) * sqrt(det(XXt %*% XX) / (n * sigmag))
  
  return(pt(piv1, n-2, lower.tail = TRUE) - pt(piv2, n-2, lower.tail = TRUE) - 0.90)
}

# Funcion para obtener c para 92%
gb192 <- function(c){
  ints <- function(beta){
    return(log(Rpb1(beta)) - log(c))
  }
  
  uno = uniroot(ints, c(b1g / 10, b1g))$root
  dos = uniroot(ints, c(b1g, b1g * 10))$root
  
  piv1 <- sqrt((n-2)/n) * (b1g - uno) * sqrt(det(XXt %*% XX) / (n * sigmag))
  piv2 <- sqrt((n-2)/n) * (b1g - dos) * sqrt(det(XXt %*% XX) / (n * sigmag))
  
  return(pt(piv1, n-2, lower.tail = TRUE) - pt(piv2, n-2, lower.tail = TRUE) - 0.92)
}

# Funcion para obtener c para 95%
gb195 <- function(c){
  ints <- function(beta){
    return(log(Rpb1(beta)) - log(c))
  }
  
  uno = uniroot(ints, c(b1g / 10, b1g))$root
  dos = uniroot(ints, c(b1g, b1g * 10))$root
  
  piv1 <- sqrt((n-2)/n) * (b1g - uno) * sqrt(det(XXt %*% XX) / (n * sigmag))
  piv2 <- sqrt((n-2)/n) * (b1g - dos) * sqrt(det(XXt %*% XX) / (n * sigmag))
  
  return(pt(piv1, n-2, lower.tail = TRUE) - pt(piv2, n-2, lower.tail = TRUE) - 0.95)
}

# Funcion para obtener c para 99%
gb199 <- function(c){
  ints <- function(beta){
    return(log(Rpb1(beta)) - log(c))
  }
  
  uno = uniroot(ints, c(b1g / 10, b1g))$root
  dos = uniroot(ints, c(b1g, b1g * 10))$root
  
  piv1 <- sqrt((n-2)/n) * (b1g - uno) * sqrt(det(XXt %*% XX) / (n * sigmag))
  piv2 <- sqrt((n-2)/n) * (b1g - dos) * sqrt(det(XXt %*% XX) / (n * sigmag))
  
  return(pt(piv1, n-2, T) - pt(piv2, n-2, T) - 0.99)
}

# Obtenemos las c's
c190 = uniroot(gb190, c(0.0001, 0.8))$root
c192 = uniroot(gb192, c(0.0001, 0.8))$root
c195 = uniroot(gb195, c(0.0001, 0.8))$root
c199 = uniroot(gb199, c(0.0001, 0.8))$root

# Obtenemos los intervalos correspondientes
ext190 = extgb1(c190)
ext192 = extgb1(c192)
ext195 = extgb1(c195)
ext199 = extgb1(c199)

# GRÁFICAS -----------------------------------------------------------

# Grafica de puntos y linea de regresion

linea <-function(x){
  return(b0g + b1g*x)
}

plot(X,linea(X),
     main = "Regresion Lineal",
     xlab = "Intensidad Luminica L",
     ylab = "Temperatura K", type = 'l', col = 'blue')
points(X,Y, pch = 20)

# Grafica PP de residuos
pp_plot(R, mg, sg, 0.95)



# Grafica de verosimilitud perfil de sigma
plot.function(x = function(t) Rps(t),
              from = 100000, #Rango
              to = 550000, lwd = 1.9,
              col = "black",
              main = expression(paste("Verosimilitud perfil de ", sigma^2)),
              ylab = expression(paste(Rp(sigma^2))),
              xlab = expression(paste(sigma^2)))
segments(x0 = uno90, y0 = 0.187, x1 = dos90, y1 = 0.187, col = "red", lwd = 2.5)
segments(x0 = uno92, y0 = 0.150, x1 = dos92, y1 = 0.150, col = "steelblue2", lwd = 2.5)
segments(x0 = uno95, y0 = 0.094, x1 = dos95, y1 = 0.094, col = "steelblue2", lwd = 2.5)
segments(x0 = uno99, y0 = 0.0165, x1 = dos99, y1 = 0.0165, col = "steelblue2", lwd = 2.5)
text(uno90-1, 0.1870, expression(c == 0.1870), cex = 0.75)
text(uno92-1, 0.1500, expression(c == 0.1500), cex = 0.75)
text(uno95-1, 0.0940, expression(c == 0.0940), cex = 0.75)
text(uno99-1, 0.0265, expression(c == 0.0165), cex = 0.75)
abline(v = sigmag, lty = 2, col = "red")
legend("topright", legend = expression(paste(hat(sigma)^2)),
       lwd = 2, lty = 2, col = "red")

#Contornos de la funci?n de verosimilitud de (mu, sigma)
plotRelative(rlogver, b0g, b1g, 3000, 4000, 1000 , 3500, c(0.01, 0.05, 0.10), 1000)
segments(x0 = 3167.568, y0 = 1000, x1 = 3167.568, y1 = 3500, col = "red",
         lwd = 2.5, lty=2)
segments(x0 = 3596.309, y0 = 1000, x1 = 3596.309, y1 = 3500, col = "red",
         lwd = 2.5, lty=2)
segments(x0 = 3000, y0 = 1780.285, x1 = 4000, y1 =1780.285 , col = "red",
         lwd = 2.5, lty=2)
segments(x0 = 3000, y0 = 2720.735, x1 = 4000, y1 = 2720.735, col = "red",
         lwd = 2.5, lty=2)
plot.function(x = function(t) emvrb1(t),
                            from = 0, #Rango
                            to = 10000, lwd = 1.5,
                           col = "blue",
                            add = T)
              plot.function(x = function(t) Iemvrb0(t),
                            from = 0, #Rango
                            to = 10000, lwd = 1.5,
                            col = "blue",
                            add = T)
points(b0g, b1g, type = "p", pch = 8, col = "red")
legend("topright", legend=expression(paste("(", hat(beta)[0],"," ,hat(beta)[1],")")),
       pch = 8, col = "red")

