LLreg_value = function(beta, X, Y){
	sigma2 = beta[length(beta)]
	beta = beta[-length(beta)]
	
	sum(log(dnorm(Y - X %*% beta, 0, sqrt(sigma2))))
	}

LLreg = function(predictors, dependant, D){
	X = as.matrix(D[,predictors])
	Y = as.vector(D[,dependant])
	
	beta_hat = solve(t(X) %*% X) %*% t(X) %*% Y
	res = Y - X %*% beta_hat
	# maxLL for variance is "/ N" not "/ N - k"
	sigma2_hat = sum(res*res) / length(Y)
	s = sqrt(sigma2_hat)
	sd = c( sqrt(diag(s**2 * solve(t(X) %*% X))), sqrt(2 * s**4 / length(Y)))
	betah = c(beta_hat, sigma2_hat)
	results = data.frame(names = c(predictors, "sigma2"), estimate = betah, sd = sd,
			t = betah / sd, p = 2*pnorm(-abs(betah / sd)))
	list(results = results, maxLL = LLreg_value(betah,X,Y))
	}
	
# example
if(FALSE){
	n = 1000
	npar = 5
	beta = rep(1,npar)
	s = 1
	X = matrix(rexp(npar * n), n, npar)
	Y = X %*% as.matrix(beta) + rnorm(n, 0, s)

	D = data.frame(X,Y)

	predictors = c("X1", "X2", "X3", "X4", "X5")
	dependant = "Y"
	LLreg(predictors, dependant, D)
	}
