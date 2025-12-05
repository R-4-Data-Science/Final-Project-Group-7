set.seed(2)
n <- 200; p <- 6
Xb <- matrix(rnorm(n*p), n, p)
linpred <- 1.2*Xb[,1] - 1*Xb[,2] + 0.8*Xb[,5]
prob <- 1 / (1 + exp(-linpred))
ybin <- rbinom(n, 1, prob)
colnames(Xb) <- paste0("x", 1:p)
dfb <- as.data.frame(cbind(y = ybin, Xb))

# Compare intercept-only vs. adding one variable by AIC:
fit0 <- glm(y ~ 1, family = binomial(), data = dfb)



# predict(fit0, )
# we will end up with a subset of models that we will use in this function to compare to the


# Assuming 'my_dataframe' is the dataframe you want to save
write.csv(wdbc, "/Users/srascon/Downloads/6210/git/Final Project/finalproject6210/breastcancer.csv", row.names = FALSE)




