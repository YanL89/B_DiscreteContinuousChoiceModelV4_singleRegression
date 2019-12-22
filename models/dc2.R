source("models/oprobit.R")

dc2 = list(
	LLVec = function(b,args){
		params = dc2$getParams(b,args)
		betaOp = params$betaOp
		betaReg = params$betaReg
		gamma = params$gamma
		sigma2 = params$sigma2
		rho = params$rho

		Yhat = args$Xreg %*% betaReg
		Zhat = args$Xop %*% betaOp
		#print(t(args$Yreg - Yhat))	
		p_reg = dnorm(args$Yreg, Yhat, sqrt(sigma2))

		# conditional mean and var of the residual of the oprobit	
		m = 0 + rho * (args$Yreg - Yhat) * 1 / sqrt(sigma2)
		s = 1 * (1 - rho**2)

		#for(i in 1:length(args$Yop))
		#	cat("Y = ", args$Yop[i], "  ", gamma[args$Yop[i]+2], "  ",
		#		gamma[args$Yop[i] + 1], "\n")
		p_pr = (pnorm(gamma[args$Yop + 2] - Zhat,m,s) 
				- pnorm(gamma[args$Yop + 1] - Zhat,m,s))

		#print(cbind(p_reg, p_pr))
		p_reg * p_pr
		},

	computeArgs = function(spec,D){
		Xop = as.matrix(D[,spec$op])
		Xreg = as.matrix(D[,spec$reg])
		Yop = as.vector(D[,spec$Yop])
		Yop = Yop - min(Yop)
		Yreg = as.vector(D[,spec$Yreg])
		list(Xop = Xop, Xreg = Xreg, Yop = Yop, Yreg = Yreg)
		},

	computeStart = function(spec,D){	
		startOP = rep(0,length(spec$op))
		ngamma = length(unique(D[,spec$Yop])) - 2
		startAlpha = 1:ngamma
		f = as.formula(paste(spec$Yreg,paste(spec$reg, collapse=" + ")
				, sep = " ~ 0 +"))
		R = lm(f,D)
		startReg = R$coefficients
		startS2 = mean(R$residuals**2)
		startRho = 0.5	
		c(startOP,startAlpha,startReg,startS2, startRho)
		},

	computeOther = function(spec, D){
		# cx names
		ngamma = length(unique(D[,spec$Yop])) - 2
		alphaNames = paste("alpha",1:ngamma,sep="_")
		paramNames = c(spec$op,alphaNames,spec$reg,"sigma2","rho")

		# bounds
		big_num = 999999
		lb=c(rep(-big_num,length(spec$op)) 
			, rep(0, ngamma)
			, rep(-big_num,length(spec$reg))
			, 0
			, -1)
		ub = c(rep(big_num,length(spec$op))
			, rep(big_num, ngamma)
			, rep(big_num, ngamma)
			, big_num
			, 1)
		list(names = paramNames, lb = lb, ub = ub)		
		},

	getParams = function(b,args){	
		start = c(1)
		start = c(start,start[1] + ncol(args$Xop))
		start = c(start,start[2] + max(args$Yop) - 1)
		start = c(start,start[3] + ncol(args$Xreg))
		start = c(start,start[4] + 1)

		betaOp = b[start[1]:(start[2] -1)]
		gamma = b[start[2]:(start[3] -1)]
		betaReg = b[start[3]:(start[4] - 1)]
		sigma2 = b[start[4]]
		rho = b[start[5]]
	
		# order gamma
		if(length(gamma) > 1)
			for(i in 2:length(gamma))
				gamma[i] = gamma[i] + gamma[i-1]
		
		gamma = c(-100,0,gamma,100)

		list(betaOp = betaOp, gamma = gamma, betaReg = betaReg,
				sigma2 = sigma2, rho = rho)
		},

	apply = function(b,spec,D){
		args = dc2$computeArgs(spec,D)
		params = dc2$getParams(b,args)
		betaOp = params$betaOp
      betaReg = params$betaReg
      gamma = params$gamma
      sigma2 = params$sigma2
      rho = params$rho
		S = matrix(c(1,rho*sqrt(sigma2),rho*sqrt(sigma2),sigma2),2,2)
		L = t(chol(S))

		n = nrow(D)
		
		err = L %*% matrix(rnorm(2*n),2,n)
		Yreg = args$Xreg %*% betaReg + err[2,]
		Z = args$Xop %*% betaOp + err[1,]

		Yop = rep(0,n)
		for(i in 2:(length(gamma) - 1))
			Yop[Z > gamma[i] ] = i -1
		
		list(Yreg = Yreg, Ydisc = Yop)
		},

	computeFitted = function(b, spec, D){
		args = dc2$computeArgs(spec,D)
		params = dc2$getParams(b,args)
		betaOp = params$betaOp
      betaReg = params$betaReg
      gamma = params$gamma
      sigma2 = params$sigma2
      rho = params$rho

		n = nrow(D)
		nalt = length(unique(args$Yop))
		
		# compute op alone
		predictors = spec$op
		dependant = spec$Yop
		OPalone = op(dependant, predictors, D)
		OPLL = OPalone$loglikelihood

		# compute regression alone
		f = paste(spec$Yreg, paste(spec$reg, collapse = " + "), sep = " ~ 0 +")
		REGalone = lm(f,data = D)
		s2 = var(REGalone$residuals)
		REGLL = sum(log(dnorm(REGalone$residuals,sd = sqrt(s2))))
	
		# predicted and fitted values
		Yhat = args$Xreg %*% betaReg
		res = D[,spec$Yop] - Yhat

		# (conditional) probabilities of each alternative and of chosen alternative
		probs = matrix(0,nrow = n, ncol = nalt)
		Z = args$Xop %*% betaOp	
		m = 0 + rho * (args$Yreg - Yhat) * 1 / sqrt(sigma2)
		s = 1 * (1 - rho**2)
		for(i in 0:(nalt-1))
			probs[,(i + 1)] =  (pnorm(gamma[i+2] - Z,m,s)
				- pnorm(gamma[i+1] - Z,m,s))
		probs = round(1000 * probs) / 1000
		chosenProb = rep(0,n)
		for(i in 1:n)
			chosenProb[i] = probs[i,args$Yop[i] + 1]

		# add the chosen alternative to the table to make it easier
		# to interpret
		probs = data.frame(probs)
		probs$chosen = args$Yop + 1

		list(fitted = Yhat, res = res, probTabe = probs, chosenProb = chosenProb,
				OPaloneLL = OPLL, REGaloneLL = REGLL)
		}
	)     
