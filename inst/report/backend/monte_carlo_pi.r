# Monte carlo simulation of pi

# Number of points (~5 sec on my system)
n <- 100000000

# Generate random points
x <- runif(n, 0, 1)
y <- runif(n, 0, 1)

# Check if the points are inside the circle
inside <- x^2 + y^2 < 1

# Calculate pi
pi_estimate <- 4 * sum(inside) / n


save( pi_estimate, file = here::here("inst", "report","backend", "pi.Rdata"))
