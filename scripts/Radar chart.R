library(fmsb)

# Carica i dati
data <- read.delim("Radar chart.txt", header = TRUE, row.names = 1, check.names = FALSE)

max_vals <- apply(data, 2, max)
min_vals <- apply(data, 2, min)

data <- rbind(max_vals, min_vals, data)

# Seleziona solo i campioni 
data_selected <- data[rownames(data) %in% c("1","2","H-MA","C-MA"), ]

# Funzione personalizzata per il radar chart
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.9,
                                        caxislabels = NULL, title = NULL, ...) {
  radarchart(
    data, axistype = 1,
    # Personalizzazione del poligono
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Personalizzazione della griglia
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Personalizzazione dell'asse
    axislabcol = "grey", 
    # Etichette delle variabili
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

op <- par(mar = c(1, 2, 2, 2), xpd = TRUE)

# Crea il radar chart
create_beautiful_radarchart(
  data = data_selected, caxislabels = c(0, 5, 10, 15, 20),
  color = c("#00AFBB", "#E7B800")
)

# Aggiungi una legenda orizzontale
legend(
  x = "bottom", legend = rownames(data_selected[-c(1, 2), ]), horiz = TRUE,
  bty = "n", pch = 20, col = c("#00AFBB", "#E7B800"),
  text.col = "black", cex = 1, pt.cex = 1.5
)

# Ripristina i margini originali
par(op)

dev.copy(svg,'Radar_Chart_MANE.svg', width = 10, height = 7)
dev.off()


# Reduce plot margin using par()
op <- par(mar = c(1, 2, 2, 2))
# Create the radar charts
create_beautiful_radarchart(
  data = df, caxislabels = c(0, 5, 10, 15, 20),
  color = c("#00AFBB", "#E7B800", "#FC4E07")
)
# Add an horizontal legend
legend(
  x = "bottom", legend = rownames(df[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c("#00AFBB", "#E7B800", "#FC4E07"),
  text.col = "black", cex = 1, pt.cex = 1.5
)
par(op)

### VERSIONE NUOVA #####

data <- read.delim("Radar chart.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Seleziona i campioni da confrontare
samples_to_compare <- data[c("H-NGE", "C-NGE"), ]

# Imposta i limiti minimi e massimi per il radar chart
max_values <- apply(data, 2, max)
min_values <- apply(data, 2, min)

# Aggiungi le righe di limiti al dataframe
samples_to_compare <- rbind(max_values, min_values, samples_to_compare)

# Colori per i campioni
colors <- c("#E7B800","#00AFBB")

# Crea il radar chart
radarchart(samples_to_compare,
           axistype = 1,
           pcol = colors,
           pfcol = sapply(colors, function(c) adjustcolor(c, alpha.f = 0.5)),
           plwd = 2,
           cglcol = "grey", 
           cglty = 1, 
           axislabcol = "black", 
           caxislabels =  seq(min(min_values), max(max_values), by = 1))

# Aggiungi una legenda
legend("topright", legend = c("H-NGE", "C-NGE"), col = colors, pch = 15, pt.cex = 1.5)

dev.copy(svg,'Radar_Chart_NGE.svg', width = 10, height = 7)
dev.off()



#One column radard chart

# Caricamento del file con separatore tab
library(fmsb)
data_raw <- read.delim("Radar chart.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Seleziona solo MA-NS
sample<- data_raw["N-GE-WS-NS", , drop = FALSE]

# Calcola i valori minimi e massimi per ogni colonna
max_values <- apply(data_raw, 2, max)
min_values <- apply(data_raw, 2, min)

# Prepara il dataset per radarchart: max, min, e il campione
samples_to_plot <- rbind(max_values, min_values, sample)

# Colore per MA-NS
color <- "#00AFBB"

# Crea il radar chart
radarchart(samples_to_plot,
           axistype = 1,
           pcol = color,
           pfcol = adjustcolor(color, alpha.f = 0.5),
           plwd = 2,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "black",
           caxislabels = seq(min(min_values), max(max_values), by = 1),
           title = "Profilo sensoriale di N-GE-WS-NS")

# Legenda
legend("topright", legend = c("N-GE-WS-NS"), col = color, pch = 15, pt.cex = 1.5)

# Salva in SVG
svg("Radar_Chart_N-GE-WS-NS.svg", width = 10, height = 7)
radarchart(samples_to_plot,
           axistype = 1,
           pcol = color,
           pfcol = adjustcolor(color, alpha.f = 0.5),
           plwd = 2,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "black",
           caxislabels = seq(min(min_values), max(max_values), by = 1),
           title = "Profilo sensoriale di N-GE-WS-NS")
legend("topright", legend = c("N-GE-WS-NS"), col = color, pch = 15, pt.cex = 1.5)
dev.off()


