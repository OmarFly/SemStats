# ==============================================================================
# Script para Análisis Bootstrap y Monte Carlo
# Autor: Gemini AI
# Fecha: 21 de mayo de 2025
#
# Descripción:
# Este script implementa los análisis solicitados en las imágenes adjuntas:
# a) Simulación Monte Carlo para estimar E(tau_hat) y V(tau_hat) con theta conocido.
# b) Bootstrap no paramétrico para una muestra observada para estimar tau,
#    V(tau_hat) y un intervalo de confianza.
# c) Estudio de simulación para evaluar la cobertura del intervalo de confianza
#    bootstrap.
# ==============================================================================

# Fijar la semilla para reproducibilidad
set.seed(1234)

# --- Funciones Auxiliares ---

# Estimador para tau(theta) = P(X=0) usando la UMVUE proporcionada
# tau_hat = ((n-1)/n)^(sum(Xi))
# Argumentos:
#   datos_muestra: un vector numérico de observaciones (Xi).
# Retorna:
#   La estimación de tau.
calcular_tau_hat <- function(datos_muestra) {
  n_val <- length(datos_muestra)
  if (n_val == 0) {
    warning("La muestra está vacía.")
    return(NA)
  }
  if (n_val == 1) { # Caso especial n=1
    if (datos_muestra[1] == 0) return(1) else return(0)
  }
  suma_xi <- sum(datos_muestra)
  return(((n_val - 1) / n_val)^suma_xi)
}

# ==============================================================================
# PARTE A: Simulación Monte Carlo (Distribución Conocida)
# ==============================================================================
# Objetivo:
#   Aproximar E(tau_hat) y V(tau_hat) mediante el método Monte Carlo (MC)
#   cuando theta es conocido.
#
# Parámetros para la simulación MC:
#   theta_verdadero_a: Parámetro lambda de la distribución Poisson.
#   n_a: Tamaño de cada muestra generada.
#   B_mc_a: Número de réplicas Monte Carlo.
# ------------------------------------------------------------------------------

cat("Iniciando Parte a: Simulación Monte Carlo...\n")

# Parámetros
theta_verdadero_a <- 1.3
n_a <- 20
B_mc_a <- 10000

# Almacenar las estimaciones de tau_hat
tau_hat_mc_valores <- numeric(B_mc_a)

# Bucle de simulación Monte Carlo
for (b in 1:B_mc_a) {
  # Generar muestra de la distribución Poisson(theta_verdadero_a)
  muestra_mc <- rpois(n = n_a, lambda = theta_verdadero_a)
  # Calcular tau_hat para la muestra generada
  tau_hat_mc_valores[b] <- calcular_tau_hat(muestra_mc)
}

# Calcular aproximaciones de E(tau_hat) y V(tau_hat)
E_tau_hat_mc <- mean(tau_hat_mc_valores)
V_tau_hat_mc <- var(tau_hat_mc_valores)

# Valor verdadero de tau(theta)
tau_verdadero_a <- exp(-theta_verdadero_a)

# Imprimir resultados de la Parte a
cat("--- Resultados Parte a ---\n")
cat("Valor verdadero de tau(theta) para theta =", theta_verdadero_a, ":", tau_verdadero_a, "\n")
cat("Aproximación Monte Carlo de E(tau_hat):", E_tau_hat_mc, "\n")
cat("Aproximación Monte Carlo de V(tau_hat):", V_tau_hat_mc, "\n")

# Generar y guardar histograma de las estimaciones de tau_hat
pdf("histograma_mc_tau_parte_a.pdf")
hist(tau_hat_mc_valores,
     main = expression(paste("Histograma de Estimaciones MC de ", hat(tau))),
     xlab = expression(hat(tau)),
     ylab = "Frecuencia",
     breaks = 50,
     col = "lightblue",
     border = "black")
abline(v = tau_verdadero_a, col = "red", lwd = 2, lty = 2)
abline(v = E_tau_hat_mc, col = "blue", lwd = 2, lty = 3)
legend("topright",
       legend = c(paste("tau(", theta_verdadero_a, ") =", round(tau_verdadero_a, 4)),
                  paste("E_MC(hat(tau)) =", round(E_tau_hat_mc, 4))),
       col = c("red", "blue"),
       lwd = 2, lty = c(2, 3))
dev.off()
cat("Histograma de la Parte a guardado como 'histograma_mc_tau_parte_a.pdf'\n")
cat("--- Fin Parte a ---\n\n")

# ==============================================================================
# PARTE B: Bootstrap No Paramétrico (Muestra Observada)
# ==============================================================================
# Objetivo:
#   Estimar P(X=0) (tau), su varianza, y un intervalo de confianza utilizando
#   el método bootstrap no paramétrico a partir de una muestra observada.
#
# Datos y Parámetros:
#   muestra_obs_b: La muestra de datos observada.
#   B_boot_b: Número de réplicas bootstrap.
#   alpha_b: Nivel de significancia para el intervalo de confianza (1 - confianza).
# ------------------------------------------------------------------------------

cat("Iniciando Parte b: Bootstrap No Paramétrico...\n")

# Muestra observada
muestra_obs_b <- c(1, 2, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0)
n_b <- length(muestra_obs_b)
B_boot_b <- 10000
alpha_b <- 0.05 # Para un IC del 95%

# Calcular tau_hat para la muestra observada original
tau_hat_obs_b <- calcular_tau_hat(muestra_obs_b)

# Almacenar las estimaciones bootstrap de tau_hat
tau_hat_boot_valores_b <- numeric(B_boot_b)

# Bucle de bootstrap no paramétrico
for (b_iter in 1:B_boot_b) {
  # Generar muestra bootstrap (muestreo con reemplazo de la muestra observada)
  muestra_boot_b <- sample(muestra_obs_b, size = n_b, replace = TRUE)
  # Calcular tau_hat para la muestra bootstrap
  tau_hat_boot_valores_b[b_iter] <- calcular_tau_hat(muestra_boot_b)
}

# Estimar la varianza de tau_hat usando bootstrap
V_tau_hat_boot_b <- var(tau_hat_boot_valores_b)

# Calcular el intervalo de confianza bootstrap (método de los percentiles)
CI_tau_boot_b <- quantile(tau_hat_boot_valores_b, probs = c(alpha_b / 2, 1 - alpha_b / 2))

# Imprimir resultados de la Parte b
cat("--- Resultados Parte b ---\n")
cat("Estimación de tau (P(X=0)) para la muestra observada (tau_hat_obs):", tau_hat_obs_b, "\n")
cat("Estimación Bootstrap de V(tau_hat):", V_tau_hat_boot_b, "\n")
cat(paste0((1 - alpha_b) * 100, "% Intervalo de Confianza Bootstrap para tau (percentiles): ["),
    CI_tau_boot_b[1], ", ", CI_tau_boot_b[2], "]\n")

# Generar y guardar histograma de las estimaciones bootstrap de tau_hat
pdf("histograma_bootstrap_tau_parte_b.pdf")
hist(tau_hat_boot_valores_b,
     main = expression(paste("Histograma de Estimaciones Bootstrap de ", hat(tau))),
     xlab = expression(paste(hat(tau),"*")),
     ylab = "Frecuencia",
     breaks = 50,
     col = "lightgreen",
     border = "black")
abline(v = tau_hat_obs_b, col = "blue", lwd = 2, lty = 2)
abline(v = CI_tau_boot_b[1], col = "purple", lwd = 2, lty = 3)
abline(v = CI_tau_boot_b[2], col = "purple", lwd = 2, lty = 3)
legend("topright",
       legend = c(paste("hat(tau)_obs =", round(tau_hat_obs_b, 4)), "Límites IC 95%"),
       col = c("blue", "purple"),
       lwd = 2, lty = c(2, 3))
dev.off()
cat("Histograma de la Parte b guardado como 'histograma_bootstrap_tau_parte_b.pdf'\n")

# Comentario sobre los resultados si la muestra proviniera de Poisson(theta=1.3)
# (valor verdadero de tau(1.3) = exp(-1.3) ~ 0.2725)
tau_verdadero_comparacion <- exp(-1.3)
cat("Comparación con tau(1.3) =", tau_verdadero_comparacion, ":\n")
cat("  La estimación observada tau_hat_obs =", tau_hat_obs_b, "\n")
cat("  El IC Bootstrap [", CI_tau_boot_b[1], ", ", CI_tau_boot_b[2], "] ",
    ifelse(tau_verdadero_comparacion >= CI_tau_boot_b[1] && tau_verdadero_comparacion <= CI_tau_boot_b[2],
           "contiene", "no contiene"), " el valor tau(1.3).\n")
cat("--- Fin Parte b ---\n\n")


# ==============================================================================
# PARTE C: Estudio de Simulación para Cobertura del IC Bootstrap
# ==============================================================================
# Objetivo:
#   Evaluar el desempeño del intervalo de confianza bootstrap no paramétrico
#   mediante un estudio de simulación para verificar su probabilidad de cobertura.
#
# Parámetros:
#   M_c: Número de repeticiones del procedimiento completo (simulaciones).
#   theta_verdadero_c: Parámetro lambda de la Poisson para generar datos.
#   n_c: Tamaño de cada muestra generada en las simulaciones.
#   B_boot_c: Número de réplicas bootstrap para cada IC.
#   alpha_c: Nivel de significancia para los ICs.
# ------------------------------------------------------------------------------

cat("Iniciando Parte c: Estudio de Simulación para Cobertura del IC...\n")
cat("ADVERTENCIA: Esta parte puede tardar varios minutos en completarse.\n")

# Parámetros
M_c <- 1000 # Número de simulaciones
theta_verdadero_c <- 1.3
n_c <- 20
B_boot_c <- 10000 # Número de réplicas bootstrap por simulación
alpha_c <- 0.05 # Para IC del 95%

# Valor verdadero de tau para la simulación
tau_verdadero_c <- exp(-theta_verdadero_c)

# Vector para almacenar los indicadores de cobertura (Z_i)
indicadores_cobertura_Z <- logical(M_c) # TRUE si cubre, FALSE si no

# Bucle principal de simulación (M repeticiones)
for (m_iter in 1:M_c) {
  if (m_iter %% (M_c/10) == 0) { # Imprimir progreso
    cat("  Simulación ", m_iter, " de ", M_c, "...\n")
  }
  
  # i. Generar n_c números aleatorios de Poisson(theta_verdadero_c)
  muestra_sim_c <- rpois(n = n_c, lambda = theta_verdadero_c)
  
  # Almacenar las estimaciones bootstrap para esta simulación
  tau_hat_boot_valores_c_actual <- numeric(B_boot_c)
  
  # ii. Bucle interno de bootstrap no paramétrico
  for (b_iter_c in 1:B_boot_c) {
    muestra_boot_c_actual <- sample(muestra_sim_c, size = n_c, replace = TRUE)
    tau_hat_boot_valores_c_actual[b_iter_c] <- calcular_tau_hat(muestra_boot_c_actual)
  }
  
  # Calcular el IC del (1-alpha_c)*100% (método de percentiles)
  CI_c_actual <- quantile(tau_hat_boot_valores_c_actual,
                          probs = c(alpha_c / 2, 1 - alpha_c / 2),
                          na.rm = TRUE) # na.rm por si alguna muestra bootstrap da problemas (poco probable aquí)
  
  # iii. Definir la variable Z
  if (!any(is.na(CI_c_actual))) { # Asegurarse que el IC se calculó correctamente
    indicadores_cobertura_Z[m_iter] <- (tau_verdadero_c >= CI_c_actual[1] && tau_verdadero_c <= CI_c_actual[2])
  } else {
    indicadores_cobertura_Z[m_iter] <- FALSE # Considerar no cobertura si el IC es NA
  }
}

# Calcular la probabilidad de cobertura empírica
cobertura_observada_c <- mean(indicadores_cobertura_Z)
nivel_confianza_nominal_c <- 1 - alpha_c

# Imprimir resultados de la Parte c
cat("--- Resultados Parte c ---\n")
cat("Valor verdadero de tau(theta) para theta =", theta_verdadero_c, ":", tau_verdadero_c, "\n")
cat("Nivel de confianza nominal:", nivel_confianza_nominal_c, "\n")
cat("Probabilidad de cobertura observada (de", M_c, "simulaciones):", cobertura_observada_c, "\n")

# Prueba de hipótesis para la cobertura (H0: p_cobertura = nivel_nominal)
# (Opcional, pero útil para el reporte)
num_coberturas <- sum(indicadores_cobertura_Z)
prueba_binomial_cobertura <- binom.test(x = num_coberturas,
                                        n = M_c,
                                        p = nivel_confianza_nominal_c)

cat("Prueba binomial para la cobertura (H0: cobertura =", nivel_confianza_nominal_c, "):\n")
cat("  Número de intervalos que cubrieron el valor verdadero:", num_coberturas, "de", M_c, "\n")
cat("  p-valor:", prueba_binomial_cobertura$p.value, "\n")

if (prueba_binomial_cobertura$p.value < alpha_c) {
  cat("  Conclusión: Se rechaza H0. La cobertura observada difiere significativamente de la nominal (alpha = ", alpha_c, ").\n")
} else {
  cat("  Conclusión: No se rechaza H0. La cobertura observada no difiere significativamente de la nominal (alpha = ", alpha_c, ").\n")
}

cat("--- Fin Parte c ---\n")
cat("Análisis completado.\n")

# Los resultados numéricos y las rutas a los archivos PDF generados
# deben ser copiados manualmente al archivo LaTeX del reporte.
# Valores clave para el reporte:
# Parte a: tau_verdadero_a, E_tau_hat_mc, V_tau_hat_mc
# Parte b: tau_hat_obs_b, V_tau_hat_boot_b, CI_tau_boot_b[1], CI_tau_boot_b[2]
# Parte c: cobertura_observada_c, prueba_binomial_cobertura$p.value

