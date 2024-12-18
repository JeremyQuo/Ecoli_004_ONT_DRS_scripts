library(waffle)
library(ggplot2)

data <- c(IN = 47, OUT = 28, other = 6)

# Create the waffle plot
temp <- waffle(data, rows = 9)

# Save the plot to a PDF file
ggsave("G.pdf", temp, width = 4, height = 4)

