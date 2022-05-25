
#######################
##### ENTEROTYPES #####
#######################

## Proyecto Final, Genomica Funcional.
## 25 05 22

## Integrantes:
## Orozco, Paloma.
## Rendon, Fernanda.
## Villa, Adrian.

### BASE DE DATOS

## Dentro de phyloseq se tienen guardadas cuatro bases de datos,
## Una de ellas es "enterotype".

library (phyloseq) # Libreria principal
data (enterotype) # Base de datos

enterotype # Objeto tipo phyloseq.
# Contiene cuadro de OTU's en 553 taxa y 280 samples.
# Las muestras (180) por 9 variables.
# Y los taxa (553) agrupados en un solo taxon (enterobacterias).

## Librerias secundarias:

library (ggplot2)
library (igraph)

theme_set (theme_bw ())
# Fondo de las graficas blanco


## ENTEROTYPE CUCHAREADO ##

genus.sum <- tapply (taxa_sums (enterotype),
                     tax_table (enterotype) [ , "Genus"],
                     sum, na.rm = T)

# De los taxa en enterotype
# Se toma solamente a los generos
# La mayoria esta clasificada hasta generos

top10genus <- names (sort (genus.sum, T)) [1:10]
# Se eligen los primeros 10 generos mas abundantes


## Bar Plots
ent10 <- prune_species (top10genus, enterotype)
# Tomar solo los 10 generos mas abundantes de enterotype

pdf ("02_Graphs/BarPlot_Cuchareado_01.pdf")

p1 <- plot_bar (ent10, "Nationality", fill = "Enterotype",
                facet_grid = ~Genus)
p1 + geom_bar (aes (color = Enterotype, fill = Enterotype),
               stat = "identity", position = "stack")
# Grafica donde se indica la abundancia de los generos mas comunes
# De acuerdo a la nacionalidad
# Y coloreado por el tipo enterico
# Sin contorno en las barras.

dev.off ()


## Plot Microbiome Network

notype <- subset_samples (enterotype, !is.na (Enterotype))
# Remover todas aquellas muestras de enterotype
# Que tengan NA en Enterotype


eig <- make_network (notype, max.dist = 0.3)
# Objeto con la red de notype y de 0.3 de distancia
# A partir de este objeto se pueden hacer redes
# Pero la distancia nunca cambiaria

pdf ("02_Graphs/Network_Cuchareado_01.pdf")

plot_network (eig, enterotype,
              color = "Enterotype", shape = "Nationality")
# Red donde se tienen todas las muestras
# Las conexiones indican que tanto se parece el bacterioma entre ellas
# El color se basa en los enterotipos
# Y la forma en la nacionalidad
# La mayoria son NAs :(

dev.off ()


ig <- make_network (enterotype, "samples", max.dist = 0.2)
# Objeto con la red de muestras y 0.2 de distancia
# A partir de este objeto se pueden hacer redes
# Pero la distancia nunca cambiaria

pdf ("02_Graphs/Network_Cuchareado_02.pdf")

plot_network (ig, enterotype, color = "ClinicalStatus",
              shape = "Enterotype")

dev.off ()

# Red con los distintos enterotipos y condiciones medicas de los sujetos (sano, obeso o NA)
# La mayoria son NAs :(


pdf ("02_Graphs/Network_Cuchareado_03.pdf")

plot_network (ig, enterotype, color = "Nationality", shape = "Enterotype")

dev.off ()

# Red con los distintos enterotipos y nacionalidades de los suetos (Danes, Japones, Espanol, NA)
# La mayoria son NAs


pdf ("02_Graphs/Network_Cuchareado_04.pdf")

plot_network (ig, enterotype, color = "ClinicalStatus",
              shape = "Nationality")

dev.off ()

# Red con las nacionalidades de los sujetos y su estado clinico
# La mayoria son NAs


## Clusters
# Se indican tres grupos de los enterotipos
# Para ver si es cierto, agrupamos por diferentes metodos

cluster_edge_betweenness (eig)
# 26 gpos

cluster_fast_greedy (eig)
# 14 gpos

cluster_infomap (eig)
# 20 gpos

cluster_label_prop (eig)
# 18 gpos

cluster_leading_eigen (eig)
# 16 gpos

cluster_leiden (eig)
# 167 gpos

cluster_louvain (eig)
# 15 gpos

# En ninguno de estos metodos se forman tres grupos
# Ni cerca



## ENTEROTYPE ORIGINAL ##

sample_data (enterotype)
# Como se observa aqui, mas de la mitad de la base de datos
# Tiene NAs en las variables
# Son aquellas "no originales".

original <- subset_samples (enterotype, !is.na (ClinicalStatus))
# 39 muestras originales

genus.sum1 <- tapply (taxa_sums (original),
                      tax_table (original) [ , "Genus"],
                      sum, na.rm = T)
# De los taxa en enterotype
# Se toma solamente a los generos
# La mayoria esta clasificada hasta generos

top10genus1 <- names (sort (genus.sum1, T)) [1:10]
# Se eligen los primeros 10 generos mas abundantes


## Bar Plots
ori10 <- prune_species (top10genus1, original)
# Tomar solo los 10 generos mas abundantes de enterotype

pdf ("02_Graphs/BarPlot_Original_01.pdf")

p2 <- plot_bar (ori10, "Nationality", fill = "Enterotype",
                facet_grid = ~Genus)
p2 + geom_bar (aes (color = Enterotype, fill = Enterotype),
               stat = "identity", position = "stack")
# Grafica donde se indica la abundancia de los generos mas comunes
# De acuerdo a la nacionalidad
# Y coloreado por el tipo enterico
# Sin contorno en las barras.

dev.off ()


## Plot Microbiome Network
notype1 <- subset_samples (original, !is.na (Enterotype))
# Remover todas aquellas muestras de original
# Que tengan NA en Enterotype

oig <- make_network (notype1, max.dist = 0.3)
# Objeto con la red de notype y de 0.3 de distancia
# A partir de este objeto se pueden hacer redes
# Pero la distancia nunca cambiara

pdf ("02_Graphs/Network_Original_01.pdf")

plot_network (oig, original,
              color = "Enterotype", shape = "Nationality")
# Red de las muestras y su parecido
# Color segun el enterotipo
# Y forma de acuerdo a la nacionalidad

dev.off ()

og <- make_network (original, "samples")

pdf ("02_Graphs/Network_Original_02.pdf")

plot_network (og, original, color = "ClinicalStatus",
              shape = "Enterotype")

dev.off ()
# Red con el enterotipo y el estado clinico de los sujetos


pdf ("02_Graphs/Network_Original_03.pdf")

plot_network (og, original, color = "Nationality", shape = "Enterotype")

dev.off ()
# Red con la nacionalidad y el enterotipo


pdf ("02_Graphs/Network_Original_04.pdf")

plot_network (og, original, color = "ClinicalStatus",
              shape = "Nationality")

dev.off ()
# Red con la nacionalidad con el estado clinico y nacionalidad 


## Clusters
# Se indican tres grupos de los enterotipos
# Agrupamos por diferentes metodos otra vez

cluster_edge_betweenness (oig)
# 7 gpos

cluster_fast_greedy (oig)
# 6 gpos

cluster_infomap (oig)
# 7 gpos

cluster_label_prop (oig)
# 7 gpos

cluster_leading_eigen (oig)
# 6 gpos

cluster_leiden (oig)
# 20 gpos

cluster_louvain (oig)
# 7 gpos

# Super disminuidos los clusteres a diferencia de
# La base de datos cuchareada



