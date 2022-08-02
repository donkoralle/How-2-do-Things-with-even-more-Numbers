# Modell 3-1: ur_metatypus // Impfungen & FPÖ
# Eerweiterung: jetzt mit UR-Detailtypen zur Charakterisierung der Cluster

## ----DatenLaden, message=FALSE, warning=FALSE---------------
library(readxl)     # Excel-Dateien lesen
library(tidyverse)  # https://www.tidyverse.org/packages/

set.seed(1976) # for reproducibility with PAM

daten <- read_excel("data/corona_gem_cluster_v2.xlsx", sheet = "ex")

# qual. Variablen setzen
daten <- daten %>%
  mutate(ur_metatypus = factor(ur_metatypus, labels = c("Urbane Räume",
                                                "Regionale Zentren",
                                                "Ländliches Umland",
                                                "Ländlicher Raum")),
         ur_typus = factor(ur_typus, labels = c("Urbane Großzentren",
                                                "Urbane Mittelzentren",
                                                "Urbane Kleinzentren",
                                                "Regionale Zentren, zentral",
                                                "Regionale Zentren, intermediär",
                                                "Ländlicher Raum im Umland von Zentren, zentral",
                                                "Ländlicher Raum im Umland von Zentren, intermediär",
                                                "Ländlicher Raum im Umland von Zentren, peripher",
                                                "Ländlicher Raum, zentral",
                                                "Ländlicher Raum, intermediär",
                                                "Ländlicher Raum, peripher")),
         degurba_21 = factor(degurba_21))

# NA-Werte checken
colSums(is.na(daten)) %>%
  knitr::kable()

# Clustervariablen setzen
# myQuantVars = c("immun_100k", "anteil_60plus", "anteil_fpoe")
myQuantVars = c("immun_100k", "anteil_fpoe")
myQualVars = c("ur_metatypus")
myVars = c(myQuantVars, myQualVars)

# Deksription der Clustervariablen
daten %>%
  select(all_of(myVars)) %>%
  summary()

# Deskription graphisch
daten %>%
  select(all_of(myQuantVars)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Messwert") %>%
  ggplot(., aes(x = "", y = Messwert)) +
  geom_boxplot() +
  labs(x = "Variablen", y = "Messwerte\n") +
  facet_wrap(~ Variable, scales = "free_y")

daten %>%
  group_by(ur_metatypus) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(. , aes(x = "", y = n, fill = ur_metatypus)) +
  geom_bar(stat = "identity", position = 'fill') +
  labs(x = "ur_typus", y = "Anteile [%]")

# Korrelationen bei metrischen Variablen
library(GGally)
ggpairs(daten, myQuantVars)

# Daten transformieren
daten_trans <- daten %>%
  select(gem_id, gem_txt, all_of(myVars)) %>%
  mutate(across(all_of(myQuantVars), scale)) %>%
  mutate_if(is.matrix, as.numeric)  # scale liefert eine Matrix und keinen numeric > wieder zurückstellen

### Gower-Distanzen ermitteln

# Vorbereitung für besseren Output: Rownames setzen
daten_trans_fit0 <- daten_trans %>%
  mutate(gem_id_rows = gem_id) %>%
  tibble::column_to_rownames("gem_id_rows")

library(cluster)
gower_dist_fit0 <- daisy(daten_trans_fit0[myVars], metric = "gower")
summary(gower_dist_fit0)

# Single Linkage Clusterung
fit0 <- hclust(gower_dist_fit0, method = "single")
plot(fit0, main = "Single Linkage Clusterung", cex = 0.65)
# kann man nicht gut interpretieren: Einzelne Gemeinden ausmachen ->

# Agglomeration Schedule mittels Funktion (braucht Labels bereits bei hclust)
# Kudos: https://stat.ethz.ch/pipermail/r-help/2013-February/348148.html

f <- function(hc){
  data.frame(row.names=paste0("Cluster",seq_along(hc$height)),
             height=hc$height,
             components=ifelse(hc$merge<0, hc$labels[abs(hc$merge)],
                               paste0("Cluster",hc$merge)),
             stringsAsFactors=FALSE)
}

tail(f(fit0), 20)

# keine Filterung nötig

### V1 Clusterung: Hierarchisch
# vgl: https://stat.ethz.ch/education/semesters/ss2012/ams/slides/v9.pdf

# Daten für finale Clusterung filtern und mit Rownames versehen
daten_trans_fit1 <- daten_trans %>%
  mutate(gem_id_rows = gem_id) %>%
  tibble::column_to_rownames("gem_id_rows")

# Clusterung anstoßen
gower_dist_fit1 <- daisy(daten_trans_fit1[myVars], metric = "gower")
  fit1 <- hclust(gower_dist_fit1, method="ward.D2")
  plot(fit1, main = "Ward-Clusterung", cex = 0.5)

# Scree Plot erstellen
# Heterogenität der letzten 10 Schritte holen
height <- sort(fit1$height)
Schritt <- c(10:1)
height <- height[(length(height)-9):length(height)]
screeplot_data_1 <- data.frame(Schritt, height)
# plotten
ggplot(screeplot_data_1, aes(x=Schritt, y=height)) +
  geom_line(size=1) +
  scale_x_continuous(breaks=Schritt) +
  labs(title = "Elbow-Diagramm", x = "Anzhal Cluster") +
  geom_vline(xintercept=6, color = "red")

# Anzahl Cluster setzen
nCluster <- 6

# Dendorgramm erstellen
plot(fit1, cex = 0.75, main = "Ward Clusterung")
rect.hclust(fit1, k = nCluster, border="red")

# Bestaz der Cluster
table(cutree(fit1, k = nCluster))

###

### V2 Clusterung: PAM

summary(gower_dist_fit1)

## Check der Ähnlichkeiten
gower_mat <- as.matrix(gower_dist_fit1)

# 2 Ähnlichsten
daten_trans_fit1[
  which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]),
        arr.ind = TRUE)[1, ], ]
# 2 Unähnlichsten
daten_trans_fit1[
  which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]),
        arr.ind = TRUE)[1, ], ]


## Anzahl Clustr bestimmen: Silhuette Width Verfahren

sil_width <- c()

# Silhuette-Width für unterschiedliche Cluster-Lösungen (2 bis 10) ermitteln
for(i in 1:10){
  pam_fit <- pam(gower_dist_fit1,
                 diss = TRUE,
                 k = i)
  sil_width[i] <- pam_fit$silinfo$avg.width
}

ggplot(as_tibble(sil_width[-1]), aes(x=2:10, y=value)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept=4, color = "red") +
  scale_x_continuous(breaks=seq(2,10,1)) +
  labs(x = "Anzahl Cluster",
       y = "Avg. Silhouette Width",
       title = "Silhuette-Plot")

# Anzahl Cluster festlegen: 6
pam_fit1 <- pam(gower_dist_fit1, diss = TRUE, k = 6)
###

### Clusterzuordnung V1 abspeichern
# in transformierte Daten eintragen
daten_trans_fit1$fit1_cl6 <- as_factor(cutree(fit1, k = nCluster))
# auf die ungefilterten transformierten Records erweitern -> Ausreißer manuell als 6 Cluster Coden
  daten_trans <- daten_trans_fit1 %>%
    select(gem_id, fit1_cl6) %>%
    left_join(daten_trans, ., by = "gem_id")
  # Join auf die absoluten Daten
  daten <- daten_trans %>%
    select(gem_id, fit1_cl6) %>%
    left_join(daten, ., by = "gem_id")
# # check
# table(daten$fit1_cl6, useNA="always")


### Clusterzuordnung V2 abspeichern
# transformiert
daten_trans_fit1$pam_fit1_cl6 <- as_factor(pam_fit1$clustering)
# transformiert total (sollte bei Clusterung gefiltert worden sein)
daten_trans <- daten_trans_fit1 %>%
  select(gem_id, pam_fit1_cl6) %>%
  left_join(daten_trans, ., by = "gem_id")
# check: table(daten_trans$pam_fit1_cl6, useNA="always")
# join auf nicht-transformierte Daten
daten <- daten_trans %>%
  select(gem_id, pam_fit1_cl6) %>%
  left_join(daten, ., by = "gem_id")
# check: table(daten$pam_fit1_cl6, useNA="always")
###


### Beurteilung der Trennschärfe der gefundenen Cluster
clusplot(daten_trans_fit1[myVars], cutree(fit1, k = nCluster),
         labels=4, color=TRUE, shade=FALSE, lines = FALSE, col.p = cutree(fit1, k = nCluster),
         main = "Bivariater Clusterplot n = 6")
# mäßig

### v2 Visualisierung: t-SNE
library(Rtsne)
tsne_obj <- Rtsne(gower_dist_fit1, is_distance = TRUE)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = daten_trans_fit1$fit1_cl6,
         name = daten_trans_fit1$gem_id)

ggplot(tsne_data, aes(x = X, y = Y)) +
  geom_point(aes(color = cluster))
###


### V1 Charakterisierung der Cluster - numerisch

daten %>%
  group_by(fit1_cl6) %>%
  summarise(across(all_of(myQuantVars), median), .groups = "keep")

### Charakterisierung graphisch: Balkendiagramme transformiert
daten_trans %>%
  group_by(fit1_cl6) %>%
  summarise(across(all_of(myQuantVars), mean), .groups="keep") %>%
  pivot_longer(all_of(myQuantVars), names_to = "variable", values_to ="wert") %>%
  ggplot(., aes(x=variable, y=wert, fill=variable)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ fit1_cl6)
# Boxplots
daten_trans %>%
  select(fit1_cl6, all_of(myQuantVars)) %>%
  pivot_longer(all_of(myQuantVars), names_to = "variable", values_to ="wert") %>%
  ggplot(., aes(x=variable, y=wert, fill=variable)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ fit1_cl6)

## charakterisierung UR-Typologie

# Grobe Sichtweise
# numerisch
table(daten$fit1_cl6, daten$ur_metatypus)
# graphisch
ggplot(daten, aes(x = fit1_cl6, fill = ur_metatypus)) +
  geom_bar(position = "fill")

# Feienere Sichtweise
table(daten$ur_typus, daten$fit1_cl6)
ggplot(daten, aes(x = fit1_cl6, fill = ur_typus)) +
  geom_bar(position = "fill")

## Charakterisierung variablenweise

ggplot(daten, aes(x = fit1_cl6, y = immun_100k, fill = fit1_cl6)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "immun_100k", x = "Cluster") +
  geom_jitter(width=0.2,alpha=0.25, color="black") +
  theme(legend.position = "none")

ggplot(daten, aes(x = fit1_cl6, y = anteil_fpoe, fill = fit1_cl6)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "anteil_fpoe", x = "Cluster") +
  geom_jitter(width=0.2,alpha=0.25, color="black") +
  theme(legend.position = "none")
###

### V2 Charakterisierung
daten_trans %>%
  group_by(pam_fit1_cl6) %>%
  summarise(across(all_of(myQuantVars), mean), .groups="keep") %>%
  pivot_longer(all_of(myQuantVars), names_to = "variable", values_to ="wert") %>%
  ggplot(., aes(x=variable, y=wert, fill=variable)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ pam_fit1_cl6)
# Boxplots
daten_trans %>%
  select(pam_fit1_cl6, all_of(myQuantVars)) %>%
  pivot_longer(all_of(myQuantVars), names_to = "variable", values_to ="wert") %>%
  ggplot(., aes(x=variable, y=wert, fill=variable)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ pam_fit1_cl6)

## charakterisierung UR-Typologie

# Grobe Sichtweise
# numerisch
table(daten$pam_fit1_cl6, daten$ur_metatypus)
# graphisch
ggplot(daten, aes(x = pam_fit1_cl6, fill = ur_metatypus)) +
  geom_bar(position = "fill")

# Feienere Sichtweise
table(daten$ur_typus, daten$pam_fit1_cl6)
ggplot(daten, aes(x = pam_fit1_cl6, fill = ur_typus)) +
  geom_bar(position = "fill")

## Charakterisierung variablenweise

ggplot(daten, aes(x = pam_fit1_cl6, y = immun_100k, fill = pam_fit1_cl6)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "immun_100k", x = "Cluster") +
  geom_jitter(width=0.2,alpha=0.25, color="black") +
  theme(legend.position = "none")

ggplot(daten, aes(x = fit1_cl6, y = anteil_fpoe, fill = fit1_cl6)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "anteil_fpoe", x = "Cluster") +
  geom_jitter(width=0.2,alpha=0.25, color="black") +
  theme(legend.position = "none")

###


### thematische Karte vorbereiten
library(sf)
library(tmap)

# Geometriedaten laden
gem <- read_sf("data/gem/STATISTIK_AUSTRIA_GEM_20210101.shp")
gem$id <- as.integer(gem$id)

# checken wie Wien abgebildet wird
# tail(gem, 30)

# Wien muss aus Union gebildet werden:
# gem_id richtig setzen
gem %>%
  filter(id > 90000)

# tail(gem_clean,30)

gem_clean <- gem %>%
  mutate(id = ifelse(id > 90000, 90001, id))

tmap::qtm(gem_clean)

# # jetzt nach id gruppieren und union laufen lassen

gem_clean <- gem_clean %>%
  group_by(id) %>%
  summarise(geometry = st_union(geometry))

# tmap::qtm(gem_clean)

# Clusterlösungen joinen: V1  V2 über daten
gem_clean_join <- left_join(gem_clean, daten,
                            by = c("id" = "gem_id"))

### Output erzeugen für V1 & V2
# V1 Karte: hierarchische Lösung
map_V31_Cl6 <- tm_shape(gem_clean_join) +
  tm_polygons(col = "fit1_cl6",
              title = "Hie 6-er Cluster-\nLösung",
              palette = "Set3",
              # palette = wes_palette("Zissou1", 5),
              legend.hist = TRUE) +
  # tm_text("id", size = 0.45, alpha = 0.5, remove.overlap = TRUE) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_legend(outside = TRUE,
            legend.outside.size = 0.15,
            hist.width = 1,
            outer.margins = 0)
map_V31_Cl6

# Karte speichern
tmap_save(map_V31_Cl6, filename = "output/hillbilly_these_V31_Cl6.png",
          units = "px", dpi = 300,
          width = 2000)

# V2 Karte: PAM-Lösung
map_V31_pamCl6 <- tm_shape(gem_clean_join) +
  tm_polygons(col = "pam_fit1_cl6",
              title = "PAM 6-er Cluster-\nLösung",
              palette = "Set3",
              # palette = wes_palette("Zissou1", 5),
              legend.hist = TRUE) +
  # tm_text("id", size = 0.45, alpha = 0.5, remove.overlap = TRUE) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_legend(outside = TRUE,
            legend.outside.size = 0.15,
            hist.width = 1,
            outer.margins = 0)
map_V31_pamCl6

# Karte speichern
tmap_save(map_V31_pamCl6, filename = "output/hillbilly_these_V31_pamCl6.png",
          units = "px", dpi = 300,
          width = 2000)
