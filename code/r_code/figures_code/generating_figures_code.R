#################################################
# Code to replicate Fig 3b, Fig 4a, Fig 4b
# Extended Data Fig 1, Extended Fig 6b and Fig 6c
#################################################

# Required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, magrittr, data.table, ggplot2, ggtern, viridis, rgdal, broom, scales, RcolorBrewer)

# Directory (all subsequent paths will be relative to the one defined below)
directory <- "/home/ltmb/Dropbox/Sebastien-Leo-Sol/time writing/code/r_code/figures_code/"



################
# Fig 3b
################


integral.kernels <- fread(paste0(directory, "KAT_share_CI.csv"))
integral.kernels <- integral.kernels[, 1:4]
setnames(integral.kernels,  c("id", "Past", "Present", "Future"))

# Create a ratio variable (future/past)
integral.kernels[, ratio := Future/Past]

# Look at the distribution 
integral.kernels %>% 
  ggplot( aes(x=ratio)) + 
  geom_histogram(bins=20, fill='skyblue', color='#69b3a2') + theme_classic()

# Keep only values that fall between 0 and 1
integral.kernels <- integral.kernels[data.table::between(Past, 0, 1) & 
                                       data.table::between(Present, 0, 1) & 
                                       data.table::between(Future, 0, 1), ]

# Plot the Ternary diagram showing the density
ggtern(integral.kernels[id == "POOL", ], aes(Past, Present, Future)) + 
  stat_density_tern(data = integral.kernels[id != "POOL", ], geom = 'polygon',
                    n         = 200,
                    aes(fill  = ..level..,
                        alpha = ..level..), bdl = 0.01) +
  geom_point() +
  scale_T_continuous(breaks= seq(0.1, 1, 0.05),labels=seq(0.1, 1, 0.05)) +
  scale_L_continuous(breaks=seq(0.1, 1, 0.05),labels=seq(0.1, 1, 0.05)) +
  scale_R_continuous(breaks=seq(0.1, 1, 0.05),labels=seq(0.1, 1, 0.05)) +
  scale_fill_viridis(name = "Frequency") +
  theme_bw() +
  theme_showsecondary() +
  theme_showarrows() +
  theme_clockwise() + 
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_text(size=15, face="bold"))

################
# Fig 4a
################


## Read the world shapefile from Natural Earth
world.shp  <- readOGR(
  dsn = paste0(directory, "ne_10m_admin_0_countries/"),
  layer = "ne_10m_admin_0_countries")

# Recode "-99" iso 2 codes 
world.shp@data[world.shp$NAME == "Norway", ]$ISO_A2 <- "NO"
world.shp@data[world.shp$NAME == "France", ]$ISO_A2 <- "FR"

# Remove Antarctica
world.shp <- world.shp[world.shp$NAME != "Antarctica", ]


# Convert world shapefile to dataframe
world.df <- broom::tidy(world.shp, region = "ISO_A2")


# Join the map dataframe with the integral kernels
world.df <- left_join(world.df, integral.kernels, by = "id")


## First, we focus on the map diplaying the ratio of future to past awareness

# find the extremes
(min.ratio <- min(world.df$ratio, na.rm = T))
(max.ratio <- max(world.df$ratio, na.rm = T))

# Get cutoffs
(breaks <- c(min.ratio, seq(0.5, 1.75, 0.25), max.ratio) %>% round(2))

# Palette of colors
red.tones <- colorRampPalette(brewer.pal(3, "Reds"))(6)
blue.tones <- colorRampPalette(brewer.pal(3, "Blues"))(8)

# Define a new variable on the data set just as above
world.df$breaks_ratio <- cut(world.df$ratio, 
                             breaks = breaks, 
                             include.lowest = TRUE)

# Add a missing observation category
levels(world.df$breaks_ratio) <- c(levels(world.df$breaks_ratio), "Not Available")
world.df$breaks_ratio[is.na(world.df$breaks_ratio)] <- "Not Available"

ratio.breaks.scale <- levels(world.df$breaks_ratio)

# Plot the map
ratio.map.plot <- ggplot() +
  geom_polygon(data = world.df, aes(fill = breaks_ratio, x = long, y = lat, group = group)) +
  geom_path(data = world.df, aes(x = long, y = lat, group = group), size = 0.1) +
  theme_void() +
  coord_equal() +
  scale_fill_manual(values = c("grey20", rev(blue.tones)[1:4], red.tones[3:5]), 
                    breaks = rev(ratio.breaks.scale), name = "Ratio of future to past attention",
                    drop = FALSE, labels = rev(ratio.breaks.scale), 
                    guide = guide_legend(direction = "horizontal", keyheight = unit(3, units = "mm"), keywidth = unit(10, units = "mm"), 
                                         title.position = 'top',
                                         nrow = 1, byrow = T, reverse = T, label.position = "bottom")) + 
  theme(legend.position = "bottom", legend.title.align = 0.5, legend.text = element_text(size = 8),
        legend.key.size = unit(80, "mm"))

# Show the ratio of attention (future/past) map
ratio.map.plot

################
# Fig 4b
################


## Now, we focus on the map displaying awareness to the present

# find the extremes
(min.present <- min(world.df$Present, na.rm = T))
(max.present <- max(world.df$Present, na.rm = T))

# Get cutoffs
(breaks <- c(min.present, seq(0.05, 0.45, 0.05), max.present) %>% round(2))

# Palette of colors
library(RColorBrewer)
ylgn.tones <- colorRampPalette(brewer.pal(9, "YlGn"))(10)


# Define a new variable on the data set just as above
world.df$breaks_present <- cut(world.df$Present, 
                               breaks = breaks, 
                               include.lowest = TRUE)

# Add a missing observation category
levels(world.df$breaks_present) <- c(levels(world.df$breaks_present), "Not Available")
world.df$breaks_present[is.na(world.df$breaks_present)] <- "Not Available"

present.breaks.scale <- levels(world.df$breaks_present)


present.map.plot <- ggplot() +
  geom_polygon(data = world.df, aes(fill = breaks_present, x = long, y = lat, group = group)) +
  geom_path(data = world.df, aes(x = long, y = lat, group = group), size = 0.1) +
  theme_void() +
  coord_equal() +
  # labs(caption = paste(uniqueN(world.df$ratio[!is.na(world.df$Present)]), "countries represented")) +
  scale_fill_manual(values = c("grey20", rev(ylgn.tones)), breaks = rev(present.breaks.scale), 
                    name = "Present attention",
                    drop = FALSE, labels = rev(present.breaks.scale), 
                    guide = guide_legend(direction = "horizontal", keyheight = unit(3, units = "mm"), 
                                         keywidth = unit(10, units = "mm"), 
                                         title.position = 'top',
                                         nrow = 1, byrow = T, reverse = T, label.position = "bottom")) + 
  theme(legend.position = "bottom", legend.title.align = 0.5, legend.text = element_text(size = 8),
        legend.key.size = unit(80, "mm"))

# Show the attention to the present map
present.map.plot

################
# ED Fig 1
################

# Read and format the data
placebo.dt <- fread(paste0(directory, "placebo_comparison.txt"))
placebo.dt[, date :=  substr(date, 1, 10) %>% as.Date()]
placebo.dt[hits == "<1", hits := "0"]
placebo.dt[, hits := as.numeric(hits)]



# Plot the time series
ggplot(data = placebo.dt[(keyword %in% c("sunscreen", "donald trump", "cat")), ], aes(x = date, y = hits)) + 
  geom_line(aes(color = keyword), alpha = 1, size = 0.6, linetype = 1)  +
  xlab("Date") +
  scale_x_date(breaks = pretty_breaks(10)) + 
  ylab("Search Volume") +
  facet_wrap(~keyword, ncol = 1) + 
  theme_classic() + 
  theme(legend.position = 'none')

  
################
# ED Fig 6b
################

# Define special points for the United States, France and Peru
special.points <- integral.kernels[id %in% c("US", "FR", "PE"), ]

ggtern(data= special.points, aes(Past, Present, Future)) + 
  geom_mask() + 
  geom_point(data = integral.kernels, color = "grey75", size = 3, shape = 16, alpha = 1) +
  geom_point(size = 3.1, color = c("dodgerblue", "firebrick", "darkolivegreen" ), shape = 1) +
  scale_T_continuous(breaks= seq(0.1, 1, 0.1),labels=seq(0.1, 1, 0.1)) +
  scale_L_continuous(breaks=seq(0.1, 1, 0.1),labels=seq(0.1, 1, 0.1)) +
  scale_R_continuous(breaks=seq(0.1, 1, 0.1),labels=seq(0.1, 1, 0.1)) +
  theme_bw() +
  theme_showarrows() +
  theme_clockwise() + 
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_text(size=15, face="bold"))

################
# ED Fig 6c
################

# Define the region each country belongs to
region.dt <- world.shp@data[, c("ISO_A2", "REGION_UN")] %>% as.data.table()
setnames(region.dt, "ISO_A2", "id", skip_absent = T)

# Merge regions with integrals
integral.kernels <- left_join(integral.kernels, region.dt, by = "id") %>% as.data.table()
integral.kernels <- integral.kernels[!is.na(REGION_UN), ]

# Attribute one color for each region
color.blind.palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

# Plot by region
for(a in 1:5){
  assign(paste0("ternary_", a), 
         ggtern(data= special.points, aes(Past, Present, Future)) + 
           geom_mask() + 
           geom_point(data = integral.kernels, color = "grey75", size = 3, shape = 16, alpha = 0.8) +
           geom_point(data = integral.kernels[REGION_UN == unique(integral.kernels$REGION_UN)[a], ], 
                      color = color.blind.palette[a], size = 3, shape = 16, alpha = 1) +
           scale_T_continuous(breaks= seq(0.1, 1, 0.1),labels=seq(0.1, 1, 0.1)) +
           scale_L_continuous(breaks=seq(0.1, 1, 0.1),labels=seq(0.1, 1, 0.1)) +
           scale_R_continuous(breaks=seq(0.1, 1, 0.1),labels=seq(0.1, 1, 0.1)) +
           labs(subtitle = unique(integral.kernels$REGION_UN)[a]) +
           theme_bw() +
           theme_showarrows() +
           theme_clockwise() + 
           theme(legend.position = "bottom",
                 axis.text=element_text(size=10),
                 plot.subtitle=element_text(size=18, face = "bold"),
                 axis.title=element_text(size=15, face="bold"))
  )
}

# Print all 5 region plots
ternary_1
ternary_2
ternary_3
ternary_4
ternary_5

