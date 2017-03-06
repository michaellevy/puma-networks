library(tidyverse)
library(stringr)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)

layers = 
  str_split(list.files("data/HomeRanges/"), "\\.") %>%
  map_chr(`[`, 1) %>% unique %>% sort
exclude = c(grep("M62", layers), length(layers))

proj = 
  lapply(layers[-exclude], function(x) {
    pdata = readOGR("data/HomeRanges/", layer = x)
    pdata@data$id = rownames(pdata@data)
    df = fortify(pdata, region = "id")
    df = merge(df, pdata@data)
    spl = toupper(str_split(x, "_")[[1]])
    df$cat = spl[1]
    df$year = if (is.null(spl[2])) NA else spl[2]
    df$sex = str_sub(df$cat, 1, 1)
    df = df[, setdiff(names(df), "Area")]  # Area appears in only one 
    
    # Calcuate centroids
    # Just on 95% isopleth. Not sure what it does without the filtering
    cent = data.frame(gCentroid(pdata[pdata$ISOPLETH == .95, ]), cat = spl[1])
    
    list(polygon = df, center = cent)
  })

# Just use most older year for two-year cats
keep = setdiff(1:15, c(2, 6, 15))

oneeach = do.call(rbind, transpose(proj)[[1]][keep]) %>%
  mutate(cat = str_replace(cat, "(M|F)0", "\\1"))

centers = 
  do.call(rbind, transpose(proj)[[2]][keep])  %>%
  mutate(cat = str_replace(cat, "(M|F)0", "\\1"))


# Get the network and coerce ggraph layout to use centroid locations for nodes
library(ggraph)
library(igraph)
uncloc = 1:0  # quantiles of lat/long to place uncollared F's. e.g. 0, 1 puts at min/max
net = readRDS("data/derived/network.RDS")
inet = intergraph::asIgraph(net)

lay = createLayout(inet, lay = "igraph", alg = "kk")
ats = attributes(lay)
lay = full_join(select(lay, -x, -y), centers, by = c("vertex.names" = "cat"))
lay = lay[, c((ncol(lay) - 1):ncol(lay), 1:(ncol(lay) - 2))]
lay$x[is.na(lay$x)] = quantile(lay$x, uncloc, na.rm = T)
lay$y[is.na(lay$y)] = quantile(lay$y, uncloc, na.rm = T)
attributes(lay) = ats

# Manually jitter some locations to clarify the visual:
move = function(cat, x, y = x, scl = diff(range(lay$x))) {
  theCat = lay$vertex.names == cat
  lay$x[theCat] = lay$x[theCat] + x * scl
  lay$y[theCat] = lay$y[theCat] + y * scl
  lay
}

# lay = move("M62", .1, .15)
lay = move("F47", .2, -.3)
lay = move("F57", -.1, 0)
lay = move("F61", -.1, 0)
lay = move("F51", 0, .09)
lay = move("M21", -.2, .4)
lay = move("M85", .04, .1)
lay = move("F49", 0, -.1)
lay = move("UncF", .3, -.3)
lay = move("UncF2", .1, .95)


baseplot = 
  filter(oneeach, ISOPLETH == .95, cat %in% c("M62", "M68", "M29", "M85")) %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = interaction(group, cat), 
                   fill = cat), alpha = .7) +
  scale_fill_brewer(palette = "Set1", name = "Territorial\nMales' Ranges") +
  coord_equal() +
  ggforce::theme_no_axes(base.theme = theme_bw(base_size = 18)) +
  geom_edge_fan(data = gEdges()(lay), aes(color = ..index..), 
                spread = 1.8, n = 300
                , edge_width = .85
                ) +
  scale_edge_color_gradient(breaks = c(.15, .85), labels = c("Sharer", "Receiver"), 
                            name = "Prey Shares", low = "black", high = "white") +
  guides(edge_colour = guide_edge_colorbar(barwidth = 0.2))  # Doesn't work!?
baseplot + 
  scale_shape_manual(values = c("M" = 24, "F" = 25), guide = "none") +
  geom_point(data = lay, mapping = aes(x = x, y = y, size = weight_kg)) +
  scale_size_continuous(range = c(4, 7), name = "Weight", breaks = c(40, 70), labels = c("40 kg", "70 kg")) + 
  geom_text(data = lay, mapping = aes(x = x, y = y, label = vertex.names)
            , size = 4, fontface = "bold"
            , color = "white" ) +
  theme(panel.background = element_rect(fill = "darkgray"))

ggsave("results/mapWithNetwork-Style1.png", height = 8, width = 10)

# OR 
baseplot + 
  geom_point(data = lay, mapping = aes(x = x, y = y, shape = sex),
             size = 2, color = "darkgray", fill = "gray") +
  scale_shape_manual(values = c("M" = 24, "F" = 25), guide = "none") +
  geom_text(data = lay, mapping = aes(x = x, y = y, label = vertex.names)
            # , vjust = 0, nudge_y = -2e3
            , vjust = "outward", hjust = "outward", 
            , size = 5, fontface = "bold")
ggsave("results/mapWithNetwork-Style2.pdf", height = 8, width = 10)


# Retired
#########

# Missing F108 ("polygon2"), M21 (M.E. omitted), M62 (the isolate) -- have them now.

ranges = 
  filter(oneeach, !cat %in% c("M62", "M68", "M85")) %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = interaction(group, cat), fill = cat), alpha = .5) +
  scale_fill_brewer(palette = "Set1", name = "") +
  # scale_fill_manual(name = "", values = cols) +
  geom_point(data = centers, mapping = aes(x = x, y = y)) +
  geom_text(data = centers, mapping = aes(x = x, y = y, label = cat)
            , vjust = "outward", hjust = "outward"
  ) +
  coord_equal() +
  ggthemes::theme_map() +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1))
ggsave("results/rangeMap.png", ranges, width = 7, height = 7)

multi = do.call(rbind, transpose(proj)[[1]])
multi$catyear = paste0(multi$cat, multi$year)
filter(multi, str_detect(cat, "(F047)|(F061)")) %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = interaction(group, catyear), fill = cat), alpha = .5) +
  facet_wrap(~catyear) +
  ggthemes::theme_map() +
  theme(legend.position = "right")
