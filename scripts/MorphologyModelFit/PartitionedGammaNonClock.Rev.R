## ----global_options, eval = TRUE, include=TRUE---------------------------
name = "f_nf"
getwd()

## ---- include=TRUE, eval = TRUE------------------------------------------
    morpho_f <- readDiscreteCharacterData("data/feedingCharacters.nex")
    morpho_nf <- readDiscreteCharacterData("data/nonFeedingCharacters.nex")


## ---- include=TRUE, eval = TRUE------------------------------------------
    taxa <- morpho_f.names()
    num_taxa <- morpho_f.size()
    num_branches <- 2 * num_taxa - 2


## ---- include=TRUE, eval = TRUE------------------------------------------
    mvi = 1
    mni = 1


## ---- include=TRUE, eval = TRUE------------------------------------------
    br_len_lambda ~ dnExp(0.2)
    moves[mvi++] = mvScale(br_len_lambda, weight=2)


## ---- include=TRUE, eval = TRUE------------------------------------------

    phylogeny ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExponential(br_len_lambda))
    moves[mvi++] = mvNNI(phylogeny, weight=num_branches/2.0)
    moves[mvi++] = mvSPR(phylogeny, weight=num_branches/10.0)
    moves[mvi++] = mvBranchLengthScale(phylogeny, weight=num_branches)
    tree_length := phylogeny.treeLength()


## ---- include=TRUE, eval = TRUE------------------------------------------
    alpha_morpho ~ dnUniform( 0, 1E6 )
    alpha_morpho2 ~ dnUniform( 0, 1E6 )

    rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )
    rates_morpho2 := fnDiscretizeGamma( alpha_morpho2, alpha_morpho2, 4 )

    #Moves on the parameters to the Gamma distribution.
    moves[mvi++] = mvScale(alpha_morpho, lambda=1, weight=2.0)
    moves[mvi++] = mvScale(alpha_morpho2, lambda=1, weight=2.0)


## ---- include=TRUE, eval = TRUE------------------------------------------
n_max_states <- 2
idx = 1
morpho_f_bystate <- morpho_f.setNumStatesVector()
for (i in 1:n_max_states) {
    nc = morpho_f_bystate[i].nchar()
    # for non-empty character blocks
    if (nc > 0) {
        # make i-by-i rate matrix
        q[idx] <- fnJC(i)
# create model of evolution for the character block
        m_morph[idx] ~ dnPhyloCTMC( tree=phylogeny,
                                    Q=q[idx],
                                    nSites=nc,
                                    siteRates=rates_morpho,
                                    type="Standard")

        # attach the data
	    m_morph[idx].clamp(morpho_f_bystate[i])

        # increment counter
        idx = idx + 1

}
}

n_max_states <- 4
morpho_nf_bystate <- morpho_nf.setNumStatesVector()
for (i in 1:n_max_states) {
    nc = morpho_nf_bystate[i].nchar()
    # for non-empty character blocks
    if (nc > 0) {
        # make i-by-i rate matrix
        q[idx] <- fnJC(i)
# create model of evolution for the character block
        m_morph[idx] ~ dnPhyloCTMC( tree=phylogeny,
                                    Q=q[idx],
                                    nSites=nc,
                                    siteRates=rates_morpho2,
                                    type="Standard")

        # attach the data
	    m_morph[idx].clamp(morpho_nf_bystate[i])
        # increment counter
        idx = idx + 1

}
}




## ---- include=TRUE, eval = TRUE------------------------------------------
    mymodel = model(phylogeny)


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors[mni++] = mnModel(filename="output/" + name + ".log", printgen=10)


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors[mni++] = mnFile(filename="output/" + name + ".trees", printgen=10, phylogeny)


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors[mni++] = mnScreen(printgen=100)


## ---- include=TRUE, eval = TRUE------------------------------------------
   mymcmc = mcmc(mymodel, monitors, moves, nruns=2, combine="mixed")
  mymcmc.run(generations=100000, tuningInterval=200)

#ss_analysis = powerPosterior(mymodel, monitors, moves, "output/" + name + "/ss", cats=20, alpha=0.3)
#ss_analysis.burnin(generations=25000,tuningInterval=100)
#ss_analysis.run(generations=25000)

#ss = steppingStoneSampler("output/" + name + "/ss", "power", "likelihood", TAB)
#ss.marginal()
### ---- include=TRUE, eval = TRUE------------------------------------------


## ---- include=TRUE, eval = TRUE------------------------------------------
#    q()


## ------------------------------------------------------------------------
