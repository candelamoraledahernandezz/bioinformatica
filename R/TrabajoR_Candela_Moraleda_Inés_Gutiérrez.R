##########################################################################################
#####################                                        \t\t\t##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       \t\t\t##################
##########################################################################################

# Colores usados en TODO el trabajo
col_control   <- "steelblue"
col_condicion <- "tomato"


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente.
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

# Si el archivo no está en el directorio de trabajo, se abre una ventana para elegirlo
archivo <- "archivo-entrada.txt"
if (!file.exists(archivo)) archivo <- file.choose()

# check.names = FALSE conserva los nombres originales de las columnas ("Gene Name", "C2-TPM"...)
datos <- read.table(archivo, header = TRUE, sep = "\t",
                    check.names = FALSE, stringsAsFactors = FALSE)

head(datos)
tail(datos)
dim(datos)   # 619 genes x 7 columnas (nombre del gen + 6 muestras)

# Columnas 2:4 = Control (C2, C1, C3); columnas 5:7 = Condición (1W5-TPM, 1W6-TPM, 1W5-COV)
idx_control   <- 2:4
idx_condicion <- 5:7


### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original.
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

datos$Media_Control   <- apply(datos[, idx_control],   1, mean)
datos$DE_Control      <- apply(datos[, idx_control],   1, sd)
datos$Media_Condicion <- apply(datos[, idx_condicion], 1, mean)
datos$DE_Condicion    <- apply(datos[, idx_condicion], 1, sd)

head(datos)
dim(datos)   # 619 x 11


### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

# Las TPM tienen valores muy extremos (algunos genes > 10000), así que se representa
# log2(media + 1) para que el histograma sea legible. Ambos usan los mismos cortes y eje X.
log_media_control   <- log2(datos$Media_Control + 1)
log_media_condicion <- log2(datos$Media_Condicion + 1)
limite_x <- range(c(log_media_control, log_media_condicion))
cortes   <- seq(floor(limite_x[1]), ceiling(limite_x[2]), by = 1)

par(mfrow = c(1, 2))
hist(log_media_control, breaks = cortes, col = col_control,
     main = "Media del Control", xlab = "log2(media TPM + 1)", ylab = "Número de genes")
hist(log_media_condicion, breaks = cortes, col = col_condicion,
     main = "Media de la Condición", xlab = "log2(media TPM + 1)", ylab = "Número de genes")
par(mfrow = c(1, 1))


### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control.
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

boxplot(list(Control = log_media_control, Condicion = log_media_condicion),
        col = c(col_control, col_condicion),
        main = "Expresión media por gen: Control vs Condición",
        ylab = "log2(media TPM + 1)")
legend("topleft", legend = c("Control", "Condición"),
       fill = c(col_control, col_condicion), bty = "n")


### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control.
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

# Test de normalidad (Shapiro-Wilk) sobre la media de cada grupo (n = 619 < 5000, válido)
shapiro_control   <- shapiro.test(datos$Media_Control)
shapiro_condicion <- shapiro.test(datos$Media_Condicion)
shapiro_control
shapiro_condicion

# Comprobación visual
par(mfrow = c(1, 2))
qqnorm(datos$Media_Control, main = "QQ-plot Control", col = col_control)
qqline(datos$Media_Control)
qqnorm(datos$Media_Condicion, main = "QQ-plot Condición", col = col_condicion)
qqline(datos$Media_Condicion)
par(mfrow = c(1, 1))

# Si ambos grupos son normales -> t-test; si no -> test no paramétrico (Wilcoxon / Mann-Whitney)
es_normal <- shapiro_control$p.value > 0.05 && shapiro_condicion$p.value > 0.05
if (es_normal) {
  cat("Datos normales -> se usa t.test\n")
} else {
  cat("Datos NO normales (p < 0.05 en Shapiro-Wilk) -> se usa wilcox.test\n")
}

# Test por gen. suppressWarnings: con valores repetidos (muchos ceros) wilcox.test avisa de empates.
# tryCatch: si un gen tiene todos los valores idénticos devuelve NA en lugar de dar error.
calcular_pvalue <- function(fila) {
  control   <- as.numeric(fila[1:3])
  condicion <- as.numeric(fila[4:6])
  tryCatch({
    if (es_normal) t.test(condicion, control)$p.value
    else           wilcox.test(condicion, control)$p.value
  }, error = function(e) NA_real_)
}

datos$p_value <- suppressWarnings(
  apply(datos[, c(idx_control, idx_condicion)], 1, calcular_pvalue)
)

head(datos)
summary(datos$p_value)


### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados
### (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan,
### ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

# log2 Fold Change (Condición / Control), con +1 para evitar dividir entre 0
datos$log2FC <- log2((datos$Media_Condicion + 1) / (datos$Media_Control + 1))

degs <- datos[!is.na(datos$p_value) & datos$p_value < 0.05, ]
cat("Número de genes diferencialmente expresados (p < 0.05):", nrow(degs), "\n")

if (nrow(degs) > 0) {
  sobreexpresados <- degs[degs$Media_Condicion > degs$Media_Control, ]
  reprimidos      <- degs[degs$Media_Condicion < degs$Media_Control, ]
  cat("Sobreexpresados:", nrow(sobreexpresados), "\n")
  cat("Reprimidos:     ", nrow(reprimidos), "\n")
  print(degs[, c("Gene Name", "Media_Control", "Media_Condicion", "p_value", "log2FC")])
} else {
  cat("\nNingún gen tiene p-value < 0.05.\n",
      "JUSTIFICACIÓN: los datos no son normales, por lo que se usa wilcox.test, que se basa en rangos.\n",
      "Con solo 3 réplicas por grupo hay únicamente choose(6,3) = 20 ordenaciones posibles; incluso en\n",
      "el caso más extremo (las 3 réplicas de un grupo por encima de las 3 del otro) el p-value exacto\n",
      "bilateral es 2/20 = 0.1. Cuando hay empates (muchos ceros) R usa la aproximación normal y el\n",
      "mínimo p-value es aprox. 0.064. En ambos casos es > 0.05, por lo que es imposible alcanzar\n",
      "significación con este diseño y este test, por muy grande que sea la diferencia de expresión.\n")
}


### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

# Si hay genes diferenciales se usan esos. Si no (caso del apartado 6), se representan como
# alternativa los 50 genes con menor p-value (desempate: mayor |log2FC|) para ver el patrón.
if (nrow(degs) >= 2) {
  genes_hm <- degs
  titulo   <- "Heatmap de genes diferencialmente expresados"
} else {
  orden    <- order(datos$p_value, -abs(datos$log2FC))
  genes_hm <- datos[head(orden, 50), ]
  titulo   <- "Heatmap: top 50 genes por p-value (ningún gen con p < 0.05)"
}

matriz <- log2(as.matrix(genes_hm[, c(idx_control, idx_condicion)]) + 1)
rownames(matriz) <- genes_hm[["Gene Name"]]
matriz <- matriz[apply(matriz, 1, sd) > 0, , drop = FALSE]   # quita genes constantes

colores_columnas <- c(rep(col_control, 3), rep(col_condicion, 3))

heatmap(matriz,
        scale = "row",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        ColSideColors = colores_columnas,
        main = titulo,
        margins = c(8, 8),
        cexRow = 0.7)
legend("topright", legend = c("Control", "Condición"),
       fill = c(col_control, col_condicion), bty = "n", cex = 0.8)


### Guardar la tabla final con todas las columnas añadidas
write.table(datos, "resultados_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
