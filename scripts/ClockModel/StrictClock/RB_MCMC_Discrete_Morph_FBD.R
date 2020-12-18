## ----global_options, eval = TRUE, include=TRUE---------------------------
name = "mk_gamma_FBD"

moves = VectorMoves()
monitors = VectorMonitors()

## ---- include=TRUE, eval = TRUE------------------------------------------
    morpho <- readDiscreteCharacterData("data/Cinctans_for_RevBayes.nex")


## ---- include=TRUE, eval = TRUE------------------------------------------
    taxa <- morpho.names()
    num_taxa <- morpho.size()
    num_branches <- 2 * num_taxa - 2

    source("scripts/Basic_FBD_model.R")


## ---- include=TRUE, eval = TRUE------------------------------------------
    alpha_morpho ~ dnUniform( 0, 1E6 )
    rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )
    #Moves on the parameters to the Gamma distribution.
    moves.append(mvScale(alpha_morpho, lambda=1, weight=2.0))
    clock_morpho ~ dnExponential(1.0)

    moves.append( mvScale(clock_morpho, lambda=0.01, weight=4.0) )
    moves.append( mvScale(clock_morpho, lambda=0.1,  weight=4.0) )
    moves.append( mvScale(clock_morpho, lambda=1,    weight=4.0) )



## ---- include=TRUE, eval = TRUE------------------------------------------
n_max_states <- 4
idx = 1
morpho_bystate <- morpho.setNumStatesVector()
for (i in 2:n_max_states) {
    # make local tmp copy of data
    # only keep character blocks with state space equal to size i
	# get number of characters per character size wth i-sized states
    nc = morpho_bystate[i].nchar()
    # for non-empty character blocks
    if (nc > 0) {
        # make i-by-i rate matrix
        q[idx] <- fnJC(i)
# create model of evolution for the character block
        m_morph[idx] ~ dnPhyloCTMC( tree=tau,
                                    Q=q[idx],
                                    nSites=nc,
                                    siteRates=rates_morpho,
                                    branchRates=clock_morpho,
                                    type="Standard")

        # attach the data
	    m_morph[idx].clamp(morpho_bystate[i])

        # increment counter
        idx = idx + 1
idx
}
}

## ---- include=TRUE, eval = TRUE------------------------------------------
    mymodel = model(tau)


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors.append(mnModel(filename="output/" +name + "log", printgen=10))


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors.append(mnFile(filename="output/" + name + "trees", printgen=10, tau))


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors.append(mnScreen(printgen=100))


## ---- include=TRUE, eval = TRUE------------------------------------------
#    mymcmc = mcmc(mymodel, monitors, moves, nruns=2, combine="mixed")

ss_analysis = powerPosterior(mymodel, monitors, moves, "output/mk_gamma_FBD/" + name + "/ss", cats=20, alpha=0.3)
ss_analysis.burnin(generations=1000,tuningInterval=100)
ss_analysis.run(generations=50000)

ss = steppingStoneSampler("output/" + name + "/mk_gamma_FBD", "power", "likelihood", TAB)
ss.marginal()
### ---- include=TRUE, eval = TRUE------------------------------------------
    mymcmc.run(generations=100000, tuningInterval=200)


## ---- include=TRUE, eval = TRUE------------------------------------------
    q()


## ------------------------------------------------------------------------
