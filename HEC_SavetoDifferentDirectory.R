a = runif(n = 10, min = 0, max = 100)
b = runif(n = 10, min = 1000, max = 1000000)
c = matrix(data = c(a,b), nrow = 10, ncol = 2)
write.csv(c, file = "test.csv")