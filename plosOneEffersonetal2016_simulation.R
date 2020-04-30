simulateFourDim <- function(numModels,d,gamma,phi,alphaMin,alphaMax,num_gamma,num_alpha,num_q1,num_phi,signalMean,sigma,maxGeneration,critOne,critTwo)  {
	
	# numModels is the number of models sampled for social learning from the previous generation
	# d is the fitness advantage associated with choosing the local optimum
	# gamma is the probability that the env changes states from one generation to the next
	# phi is the probability that the signal correctly indicates when models and learners have the same optimum and also the probability that the signal correctly indicates when models and learners have different optima
	# alphaMin is the minimum value of alpha considered when creating the genotype space
	# alphaMax is analogous
	# num_gamma is the number of gamma values considered => creating a lattice for genotype space
	# num_alpha is the number of alpha values considered => creating a lattice for genotype space
	# num_q1 is the number of q1 values considered => creating a lattice for genotype space
	# signalMean is the mean for the signal distribution when the state is 1, -signalMean holds when state is 0
	# maxGeneration is the maximum possible generation, which only comes into play if a stopping criteria is not met sooner
	# critOne is a threshold frequency such that, if a single genotype has a frequency >= critOne, the simulation is stopped
	# critTwo is a threshold number of generations such that, if a single model genotype has persisted as the single modal genotype for >= critTwo generations, the simulation is stopped

	# Create genotype space
	alpha_temp <- seq(from = alphaMin, to = alphaMax, length.out = num_alpha)
	gamma_temp <- seq(from = 0.0001, to = 0.9999, length.out = num_gamma)
	q1_temp <- seq(from = 0.0001, to = 0.9999, length.out = num_q1)
	phi_temp <- seq(from = 0.0001, to = 0.9999, length.out = num_phi)

	temp <- expand.grid(alpha_temp,gamma_temp,q1_temp,phi_temp)

	geno <- data.frame(temp[ ,1],temp[ ,2],temp[ ,3],temp[ ,4])
	names(geno) <- c('alpha_g','gamma_g','q1_g','phi_g')

	rm(alpha_temp,gamma_temp,q1_temp,phi_temp,temp)
	
	# Construct matrices for relevant matrix operations.
	critMatrix_aZero <- rep(0,nrow(geno) * (numModels + 1)) # The matrix for the case when a = 0 => signal indicates models and learners do not have the same optimum
	dim(critMatrix_aZero) <- c(nrow(geno),(numModels + 1))
	
	critMatrix_aOne <- rep(0,nrow(geno) * (numModels + 1)) # The matrix for the case when a = 1 => signal indicates models and learners have the same optimum
	dim(critMatrix_aOne) <- c(nrow(geno),(numModels + 1))
	
	alpha_g <- geno$alpha_g; dim(alpha_g) <- c(nrow(geno),1)
	gamma_g <- geno$gamma_g; dim(gamma_g) <- c(nrow(geno),1)
	q1_g <- geno$q1_g; dim(q1_g) <- c(nrow(geno),1)
	phi_g <- geno$phi_g; dim(phi_g) <- c(nrow(geno),1)
	probChooseOne <- rep(-99,nrow(geno)); dim(probChooseOne) <- c(nrow(geno),1)
	fitness <- rep(-99,nrow(geno)); dim(fitness) <- c(nrow(geno),1)
	frequency <- rep(-99,nrow(geno))
	frequency[1:(nrow(geno) - 1)] <- 1 / nrow(geno)
	frequency[nrow(geno)] <- 1 - sum(frequency[1:(nrow(geno) - 1)])
	dim(frequency) <- c(nrow(geno),1)

	rm(geno)

	# Population-level quantities
	gen <- 1:maxGeneration # Generation
	t <- 1
	finalGeneration <- 0
	newState <- 0 # A dummy indicating if the state has changed w.r.t. the last generation

	criterionOne <- rep(0,length(gen)) # Criteria for stopping the simulation
	criterionTwo <- rep(0,length(gen))

	population <- data.frame(gen,criterionOne,criterionTwo)

	rm(gen,criterionOne,criterionTwo)

	# Initialize state and other population-level variables
	population$state <- c(round(runif(1)),rep(-99,nrow(population) - 1))
	population$meanChoice <- rep(-99,nrow(population))
	population$meanFitness <- rep(-99,nrow(population))
	population$mean_alpha_g <- rep(-99,nrow(population))
	population$mean_gamma_g <- rep(-99,nrow(population))
	population$mean_q1_g <- rep(-99,nrow(population))
	population$mean_phi_g <- rep(-99,nrow(population))

	idxModalGenotype <- -99 # Will only take a positive index value if a single genotype has the highest frequency, => no ties admitted
	numGenModalGenotype <- 0 # Tracks the number of generations with the same modal genotype

	while (t <= maxGeneration)
	{
		if (population$criterionOne[t] == 0 & population$criterionTwo[t] == 0)
		{
			# Probabilities of choosing 1
			if (t == 1) # Only individual learning
			{
				probChooseOne[ ,1] <- population$state[t] * (1 - pnorm(0,mean = signalMean,sd = sigma)) + (1 - population$state[t]) * (1 - pnorm(0,mean = -signalMean,sd = sigma))
			}
			else
			{
				for (i in 0:numModels)
				{
					# The matrix for when the signal indicates models and learners have a different state (A = 0)
					critMatrix_aZero[ ,i+1] <- 1 - pnorm(alpha_g * (log(q1_g^(numModels - i) * (1 - q1_g)^i * (1 - gamma_g) * (1 - phi_g) + q1_g^i * (1 - q1_g)^(numModels - i) * gamma_g * phi_g) - log(q1_g^(numModels - i) * (1 - q1_g)^i * gamma_g * phi_g + q1_g^i * (1 - q1_g)^(numModels - i) * (1 - gamma_g) * (1 - phi_g))), mean = (1 - population$state[t]) * (-1) * signalMean + population$state[t] * signalMean, sd = sigma)
					
					# The matrix for when the signal indicates models and learners have the same state (A = 1)
					critMatrix_aOne[ ,i+1] <- 1 - pnorm(alpha_g * (log(q1_g^(numModels - i) * (1 - q1_g)^i * (1 - gamma_g) * phi_g + q1_g^i * (1 - q1_g)^(numModels - i) * gamma_g * (1 - phi_g)) - log(q1_g^(numModels - i) * (1 - q1_g)^i * gamma_g * (1 - phi_g) + q1_g^i * (1 - q1_g)^(numModels - i) * (1 - gamma_g) * phi_g)), mean = (1 - population$state[t]) * (-1) * signalMean + population$state[t] * signalMean, sd = sigma)
				}
				if (newState == 0) # Models and learners have the same state
				{
					probChooseOne <- phi * (critMatrix_aOne %*% dbinom((0:numModels),numModels,population$meanChoice[t-1])) + (1 - phi) * (critMatrix_aZero %*% dbinom((0:numModels),numModels,population$meanChoice[t-1]))		
				}	
				else # Models and learners have a different state
				{
					probChooseOne <- (1 - phi) * (critMatrix_aOne %*% dbinom((0:numModels),numModels,population$meanChoice[t-1])) + phi * (critMatrix_aZero %*% dbinom((0:numModels),numModels,population$meanChoice[t-1]))		
				}	
			}
			
			fitness <- 1 + d * (population$state[t] * probChooseOne + (1 - population$state[t]) * (1 - probChooseOne))
			
			population$meanChoice[t] <- t(frequency) %*% probChooseOne
			population$meanFitness[t] <- t(frequency) %*% fitness
			population$mean_alpha_g[t] <- t(frequency) %*% alpha_g
			population$mean_gamma_g[t] <- t(frequency) %*% gamma_g
			population$mean_q1_g[t] <- t(frequency) %*% q1_g
			population$mean_phi_g[t] <- t(frequency) %*% phi_g
			
	
			# Update genotype frequencies
			frequency <- frequency * fitness / (population$meanFitness[t]) 
			
			# Check for stopping criteria
			if (t == 1)
			{
				if (sum(frequency == max(frequency)) == 1)
				{
					idxModalGenotype <- which.max(frequency)
					numGenModalGenotype <- 1
				}
			}
			else
			{
				if (which.max(frequency) == idxModalGenotype & sum(frequency == max(frequency)) == 1)
				{
					idxModalGenotype <- which.max(frequency)
					numGenModalGenotype <- numGenModalGenotype + 1
				}
				else if (which.max(frequency) != idxModalGenotype & sum(frequency == max(frequency)) == 1)
				{
					idxModalGenotype <- which.max(frequency)
					numGenModalGenotype <- 1
				}
				else
				{
					idxModalGenotype <- -99
					numGenModalGenotype <- 0
				}
			}
			
			if (max(frequency) >= critOne)
			{
				population$criterionOne[t] <- 1
			}
			if (numGenModalGenotype >= critTwo)
			{
				population$criterionTwo[t] <- 1
			}
			if (population$criterionOne[t] == 1 | population$criterionTwo[t] == 1) # stop simulation and record final state in geno data frame
			{
				finalGeneration <- t
				t <- maxGeneration + 1
			}
	
			# State changes
			if (runif(1) <= gamma & t < maxGeneration)
			{
				population$state[t+1] <- abs(population$state[t] - 1)
				newState <- 1
			} 
			else
			{
				if (t < maxGeneration)
				{
					population$state[t+1] <- population$state[t]
					newState <- 0
				}
			}
			
			if (t <= maxGeneration)
			{
				t <- t + 1
				print(t)
			}
		}
		else
		{
			t <- maxGeneration + 1
		}
	}
	
	if (finalGeneration == 0)
	{
		finalGeneration <- maxGeneration
	}
	
	geno <- data.frame(alpha_g,gamma_g,q1_g,phi_g,probChooseOne,fitness,frequency)
	geno$finalGen <- rep(finalGeneration,nrow(geno))

	population <- population[population$state > -99, ]
	
	out <- list(geno,population)
	names(out) <- c('geno','population')

	out

}