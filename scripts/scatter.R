# Scatter plot
plot(x, y,
     pch = 19,
     col = factor(group))

# Legend
legend("topleft",
       legend = levels(factor(group)),
       pch = 19,
       col = factor(levels(factor(group))))

# Color selection
colors <- c("#FDAE61", # Orange
            "#D9EF8B", # Light green
            "#66BD63") # Darker green

# Scatter plot
plot(x, y,
     pch = 19,
     col = colors[factor(group)])

# Legend
legend("topleft",
       legend = c("Group 1", "Group 2", "Group 3"),
       pch = 19,
       col = colors)