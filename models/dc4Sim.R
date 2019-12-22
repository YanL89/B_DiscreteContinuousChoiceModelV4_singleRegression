
dc4Sim = list(
	LLVec = function(b, args){
		attach(args)
		attach(dc4Params(b, nparPr = ncol(XPr), nparReg = ncol(XReg)))
		S = L %*% t(L)
		
		# LL for the regression
		epsilon = YReg - (XReg %*% as.matrix(betaReg))
		sigma2 = S[nrow(S), ncol(S)]
		LLreg = sum(log(dnorm(epsilon, 0, sigma2)))
		
		# utilities
		U = matrix(XPr %*% betaPr, nrow = n, ncol = nalt, byrow = TRUE)
		
		# we simulate for the LL of the probit, 
		# we need a loop because at every line we need to reparametrize
		LLPr = 0
		for(i in 1:ncol(args$YPr)){
			SigmaDiffCond = condCov(S, args$nalt)
			muDiffCond = condMean(S, mu = rep(0, args$nalt, posObs = nalt
					, obs = epsilon[i]))
			
			if(args$method = "simulation"){
				LLPr = LLPr + log(simProbitL(U[i,], args$YPr[i], SigmaDiffCond
						, muDiffConf, args$Z[[i]], args$nalt))
				} else {
				# call pmvnorm
				}
				
			}
		LLReg + LLPr 
		},
		

		
	computeStart = function(spec, D){
		dc4$computeStart(spec, D)
		},
		
	computeOther = function(spec, D){
		dc4$computeOther(spec,D)
		}
		
	)	
