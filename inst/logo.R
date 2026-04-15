library(ggplot2)
library(hexSticker)
library(sysfonts)
library(showtext)
library(ggtext)

font_add_google("Sora", "sora")
showtext_auto()
showtext_opts(dpi = 300)


n_points <- 200
r_min <- 1.5
r_max <- 5.5
n_turns <- 12  # many turns = looks like circles

radii <- seq(r_min, r_max, length.out = n_points)

# Uniform angular spacing (constant angle increment)
angles <- seq(0, n_turns * 2 * pi, length.out = n_points)

df <- data.frame(
  x = radii * cos(angles),
  y = radii * sin(angles),
  radius = radii,
  size = rep(1 * c(0.5, 1, 2, 1.5, 2.5, 1, 3, 2), length.out = n_points)
)
df$col_val <- df$radius

p <- ggplot(df, aes(x = x, y = y, size = size, colour = col_val)) +
  geom_point(alpha = 0.65, stroke = 0) +
  coord_fixed() +
  scale_colour_gradient(low = "#F59E0B", high = "#60A5FA")+
  scale_size_identity() +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA)
  ) +
  guides(size = "none", colour = "none")

p

# --- Sticker ---
library(ggtext)

s <- sticker(
  p,
  package = "",
  s_x = 1,
  s_y = 1,
  s_width = 1.8,
  s_height = 1.8,
  h_fill = "#020617",
  h_color = "#3B82F6",
  h_size = 1.4,
  filename = "resemble_hex.png",
  dpi = 300, 
  url = "https://github.com/l-ramirez-lopez/resemble", 
  u_color = "white", 
  u_size = 1.1
)

s <- s +
  geom_richtext(
    aes(x = 1, y = 1),
    label = "<span style='color:#F8FAFC;'>rese</span><span style='color:#F59E0B;'>mbl</span><span style='color:#F8FAFC;'>e</span>",
    family = "outfit",
    fontface = "bold",
    size = 9,
    fill = NA,
    label.color = NA
  )

ggsave("man/figures/logo.png", s, width = 43.94, height = 50.8, units = "mm", dpi = 300, bg = "transparent")
