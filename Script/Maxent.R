#11/07/2024
#Modelamiento MaxENT para distribución espacial de Aedes aegypti en el Perú


# Cargar de paquetes -------------------------------------------------------
install.packages("rJava")
options(java.parameters = "-Xmx4g" )
library(sp)
library(raster)
library(dismo)
library(rJava) #https://www.java.com/en/download/manual.jsp
library(rgeos)
library(grDevices)
library(rgdal)
library(maptools)
library(sp)
library(readr)
library(corrplot)
library(xlsx)
library(MASS)
library(rasterVis)
library(RColorBrewer)
library(ggplot2)
library(foreign)
library(apaTables)
library(PerformanceAnalytics)
library(psych)
library(corrr)
library(pROC)


# Ruta de archivos --------------------------------------------------------
#Carapeta de trabajo
setwd('D:/BASE DE DATOS/BASE DE DATOS/VARIABLES/')
getwd()
#Datos de presencia
occs <- read.csv("Dengue1.csv")
#occs <- occs[,-1]
#write.csv(occs,"dengue_1_correr.csv")
#View(occs)

#Capas ambientales presente
layers.p <- stack(list.files('Wordlclim/',"*.asc$",full.names=T))
list(layers.p)
plot(layers.p)
names(layers.p) <- c('BIO1','BIO10','BIO11','BIO12','BIO13','BIO14','BIO15','BIO16','BIO17','BIO18','BIO19','BIO2','BIO3','BIO4',
                     'BIO5','BIO6','BIO7','BIO8','BIO9','elev')
names(layers.p)

#layers.f <- stack(list.files("Worlclim/capas ambientales futuras/","*.asc$",full.names=T))
#plot(layers.p)


# Limpieza y selección de variables ---------------------------------------
points(occs[,2:3], pch=7, cex=0.3)
datafra<- as.data.frame(occs, check.names = TRUE)
pres.covs<-terra::extract(layers.p, datafra[,2:3],cellnumbers=T,check.names = TRUE)
pres.covs<-na.omit(pres.covs)
pres.covs<-unique(pres.covs)
pres.covs<-pres.covs[,-1]
write.table(pres.covs, file = "ACTUAL/Dengue_process.csv",sep = ',', row.names = F) # guarda un archivo excel
View(pres.covs)

#Correlación de variables
datos.v <- read_csv('ACTUAL/Dengue_process.csv')
datos.v <- as.data.frame(datos.v)

#View(datos.v)
png("ACTUAL/Correlación.png", units = 'cm', width = 38, height = 30, res = 1200)
corr <- corrplot(cor(datos.v), type = "lower", method = "number")
dev.off()

#Exportar a word
apa.cor.table(datos.v, filename = 'ACTUAL/Tabla 1.doc', table.number = 1, show.conf.interval = FALSE, landscape = TRUE)

#regresar a caragr las variables linea 44
#Seleccion de variables a excluir (Agregar variables con correlación baja)
layers.p <- dropLayer(layers.p, c('BIO19','BIO18','BIO17', 'BIO16','BIO12','BIO9','BIO6')) #'BIO4','BIO6','BIO7','BIO12','BIO13','BIO14','BIO17'
names(layers.p)

#Ploteo de las variables y puntos de presencia
plot(layers.p[[1,]])
points(occs[,2:3], pch=7, cex=0.3)
pres.covs<-terra::extract(layers.p, datafra[,2:3],cellnumbers=T,check.names = TRUE)
pres.covs<-na.omit(pres.covs)
pres.covs<-unique(pres.covs)
pres.covs<-pres.covs[,-1]

write.table(pres.covs, file = "ACTUAL/dengue_final.csv",sep = ',', row.names = F)

#Correlación de variables
datos.v <- read_csv('ACTUAL/dengue_final.csv')
datos.v <- as.data.frame(datos.v)

#View(datos.v)
png("ACTUAL/Correlación_seleccionadas.png", units = 'cm', width = 38, height = 30, res = 1200)
corr <- corrplot(cor(datos.v), type = "lower", method = "number")
dev.off()

#Exportar a word
apa.cor.table(datos.v, filename = 'ACTUAL/Tabla 2.doc', table.number = 1, show.conf.interval = FALSE, landscape = TRUE)


#creamos la matriz con valores de pseudoausencias
bkg.covs<-sampleRandom(layers.p,10000,cells=T)
bkg.covs<-unique(bkg.covs)

#eliminamos la columna del numero de celdas
bkg.covs<-bkg.covs[,-1] 
#View(bkg.covs)

#Dividira los datos en datos de prueba y validacion mediante la tecnica de kfold partitioning
#Genera un indice aleatorio de los folds
fold <- kfold(pres.covs, k=8) 
occtest <- pres.covs[fold == 1, ]
occtrain <- pres.covs[fold != 1, ]

#View(occtrain)
y<-c(rep(1,nrow(occtrain)), rep(0,nrow(bkg.covs)))
env.values<-data.frame(rbind(occtrain, bkg.covs))

# Entrenar el modelo MaxEnt
me <- maxent(env.values, y, args=c("addsamplestobackground=true"), 
             path="D:/INFORME_PRO_FORESTAL/INVESTIGACIÓN/DENGUE/BASE DE DATOS/")

# Evaluar el modelo con los datos de Entrenamiento
J <- evaluate(me, p = data.frame(occtrain), a = data.frame(bkg.covs))
str(J)
auc_training <- round(J@auc, 2)

# Ploteo de la curva AUC - ROC para datos de Entrenamiento
png("AUC_training.png", units = 'cm', width = 38, height = 25, res = 1200)

# Ajustar los márgenes para centrar el gráfico
par(mar = c(5, 5, 4, 2) + 0.1)

## Graficar la curva ROC con grillas
plot((1 - J@TNR), J@TPR, type = "l",
     xlab = "False Positive Rate (1 - Specificity)", cex.lab = 1.5, cex.axis = 1.5,
     ylab = "True Positive Rate", cex.lab = 1.5, cex.axis = 1.5,
     xlim = c(0, 1), ylim = c(0, 1),
     col = "blue", xaxs = "i", yaxs = "i", lwd = 4)

grid(col = "gray", lty = "dotted")

# Calcular los valores de los ejes X e Y para la línea diagonal
x_values <- seq(0, 1, length.out = 100)
y_values <- x_values

# Trazar la línea diagonal
lines(x_values, y_values, col = "red", lty = 2, lwd = 4)

# Calcular la pendiente de la línea de regresión lineal
regression_slope <- round(lm(J@TPR ~ (1 - J@TNR))$coefficients[2], 2)

# Agregar leyenda
legend_text <- c(paste("ROC Curve (AUC =", auc_training, ")", sep = ""),
                 "Diagonal Line (Slope = 1)")

legend("bottomright", legend = legend_text, col = c("blue", "red"), lty = c(1, 2), lwd = c(4, 4), cex = 1.5)

# Agregar título
title("AUC - ROC Training", font.main = 2, cex.main = 2)

dev.off()

# Evaluar el modelo con los datos de Validación
e <- evaluate(me, p = data.frame(occtest), a = data.frame(bkg.covs))
str(e)
auc_testing <- round(e@auc, 2)

# Ploteo de la curva AUC - ROC para datos de Validación
png("AUC_testing.png", units = 'cm', width = 38, height = 25, res = 1200)

# Ajustar los márgenes para centrar el gráfico
par(mar = c(5, 5, 4, 2) + 0.1)

## Graficar la curva ROC con grillas
plot((1 - e@TNR), e@TPR, type = "l",
     xlab = "False Positive Rate (1 - Specificity)", cex.lab = 1.5, cex.axis = 1.5,
     ylab = "True Positive Rate", cex.lab = 1.5, cex.axis = 1.5,
     xlim = c(0, 1), ylim = c(0, 1),
     col = "blue", xaxs = "i", yaxs = "i", lwd = 4)

grid(col = "gray", lty = "dotted")

# Calcular los valores de los ejes X e Y para la línea diagonal
x_values <- seq(0, 1, length.out = 100)
y_values <- x_values

# Trazar la línea diagonal
lines(x_values, y_values, col = "red", lty = 2, lwd = 4)

# Calcular la pendiente de la línea de regresión lineal
regression_slope <- round(lm(e@TPR ~ (1 - e@TNR))$coefficients[2], 2)

# Agregar leyenda
legend_text <- c(paste("ROC Curve (AUC =", auc_testing, ")", sep = ""),
                 "Diagonal Line (Slope = 1)")

legend("bottomright", legend = legend_text, col = c("blue", "red"), lty = c(1, 2), lwd = c(4, 4), cex = 1.5)

# Agregar título
title("AUC - ROC Testing", font.main = 2, cex.main = 2)

dev.off()

# Ploteo de las curvas AUC - ROC para datos de Entrenamiento y Validación
png("AUC_combined.png", units = 'cm', width = 38, height = 25, res = 1200)

# Ajustar los márgenes para centrar el gráfico
par(mar = c(5, 5, 4, 2) + 0.1)

## Graficar la curva ROC con grillas para Entrenamiento
plot((1 - J@TNR), J@TPR, type = "l",
     xlab = "False Positive Rate (1 - Specificity)", cex.lab = 1.5, cex.axis = 1.5,
     ylab = "True Positive Rate", cex.lab = 1.5, cex.axis = 1.5,
     xlim = c(0, 1), ylim = c(0, 1),
     col = "blue", xaxs = "i", yaxs = "i", lwd = 4)

## Graficar la curva ROC para Validación
lines((1 - e@TNR), e@TPR, col = "green", lwd = 4)

grid(col = "gray", lty = "dotted")

# Calcular los valores de los ejes X e Y para la línea diagonal
x_values <- seq(0, 1, length.out = 100)
y_values <- x_values

# Trazar la línea diagonal
lines(x_values, y_values, col = "red", lty = 2, lwd = 4)

# Agregar leyenda
legend_text <- c(paste("Training ROC Curve (AUC =", auc_training, ")", sep = ""),
                 paste("Testing ROC Curve (AUC =", auc_testing, ")", sep = ""),
                 "Diagonal Line (Slope = 1)")

legend("bottomright", legend = legend_text, col = c("blue", "green", "red"), lty = c(1, 1, 2), lwd = c(4, 4, 4), cex = 1.5)

# Agregar título
title("AUC - ROC Training vs Testing", font.main = 2, cex.main = 2)

dev.off()

