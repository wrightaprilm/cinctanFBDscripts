
taxa <- readTaxonData("data/cincta_fas_fuzzy.tsv")
taxa
morpho_f <- readDiscreteCharacterData("data/feedingCharacters.nex")
morpho_nf <- readDiscreteCharacterData("data/nonfeedingCharacters.nex")
moves = VectorMoves()
monitors = VectorMonitors()

n_taxa <- taxa.size()

num_branches <- 2 * n_taxa - 2


fbd_indicator ~ dnCategorical(simplex(1, 1, 1))
moves.append( mvRandomGeometricWalk(fbd_indicator, weight=10.0, tune=FALSE) )
timeline1 <- v(10.8);
timeline2 <- v(.5, 4.5, 9)
timeline3 <- v(0.9, 2.3, 3.9, 6.3, 8.4, 10.8);
fbd_vec := v(timeline1, timeline2, timeline3)
fbd_indicator.setValue(2)
timeline := fbd_vec[fbd_indicator]
print(timeline)

for(i in 1:(timeline.size()+1)){

    speciation_rate[i] ~ dnExponential(1.471);
    moves.append(mvScale(speciation_rate[i], lambda=0.01, weight=5));
    moves.append(mvScale(speciation_rate[i], lambda=0.10, weight=3));
    moves.append(mvScale(speciation_rate[i], lambda=1.00, weight=1));
    speciation_rate[i]
    turnover[i] ~ dnLnorm(ln(0.945), 0.6564);                   # dnUnif(0.9, 1.05);
    moves.append(mvSlide(turnover[i], delta=0.01, weight=5));
    moves.append(mvSlide(turnover[i], delta=0.10, weight=3));
    moves.append(mvSlide(turnover[i], delta=1.00, weight=1));

    extinction_rate[i] := turnover[i]*speciation_rate[i]
    diversification[i] := speciation_rate[i] - extinction_rate[i]

    psi[i] ~ dnExponential(3.892);
    moves.append( mvScale(psi[i], lambda = 0.01) )
    moves.append( mvScale(psi[i], lambda = 0.1) )
    moves.append( mvScale(psi[i], lambda = 1) )
}

      # Proportional Taxon Sampling of Youngest Time Slice
      rho <- 0.506; # 'extant' sampling.

     # Establish Basal Divergence Time
     origin_time ~ dnUnif(7.3, 12.11);
     moves.append(mvSlide(origin_time, delta=0.01, weight=5));
     moves.append(mvSlide(origin_time, delta=0.10, weight=3));
     moves.append(mvSlide(origin_time, delta=1.00, weight=1));

     fbd_dist = dnFBDP(originAge=origin_time, lambda=speciation_rate, mu=extinction_rate, psi=psi, rho=rho, timeline=timeline, taxa=taxa, condition="sampling")


     outgroup = clade("Ctenocystis");
     ingroup = clade("Gyrocystis_platessa","Gyrocystis_testudiformis","Gyrocystis_cruzae","Gyrocystis_badulesiensis","Gyrocystis_erecta","Progyrocystis_disjuncta","Protocinctus_mansillaensis","Elliptocinctus_barrandei","Elliptocinctus_vizcainoi","Sucocystis_theronensis","Sucocystis_bretoni","Lignanicystis_barriosensis","Undatacinctus_undata","Sucocystis_acrofera","Undatacinctus_quadricornuta","Undatacinctus_melendezi","Asturicystis_jaekeli","Sotocinctus_ubaghsi","Trochocystites_bohemicus","Trochocystoides_parvus","Ludwigicinctus_truncatus","Graciacystis_ambigua","Asturicystis_havliceki","Nelegerocystis_ivantzovi","Rozanovicystis_triangularis","Davidocinctus_pembrokensis");

    constraints = v(ingroup)


     fbd_tree ~ dnConstrainedTopology(fbd_dist, constraints=constraints)


     moves.append(mvFNPR(fbd_tree, weight=15.0))

     moves.append(mvCollapseExpandFossilBranch(fbd_tree, origin_time, weight=6.0))

     moves.append(mvNodeTimeSlideUniform(fbd_tree, weight=40.0))

     moves.append(mvRootTimeSlideUniform(fbd_tree, origin_time, weight=5.0))


 # Setup the fossil tip sampling #

 # Use a for loop to create a uniform distribution on the occurence time for each fossil #

 # The boundaries of the uniform distribution are specified in the tsv file #

 fossils = fbd_tree.getFossils()
 for(i in 1:fossils.size())
 {
     t[i] := tmrca(fbd_tree, clade(fossils[i]))

     a_i = fossils[i].getMinAge()
     b_i = fossils[i].getMaxAge()

     F[i] ~ dnUniform(t[i] - b_i, t[i] - a_i)
     F[i].clamp( 0 )
 }

 # Add a move to sample the fossil times #
 moves.append( mvFossilTimeSlideUniform(fbd_tree, origin_time, weight=5.0) )


     num_samp_anc := fbd_tree.numSampledAncestors()

     pruned_fbd_tree := fnPruneTree(fbd_tree, prune=v("Asturicystis_havliceki","Nelegerocystis_ivantzovi","Rozanovicystis_triangularis","Davidocinctus_pembrokensis"))

clock_indicator ~dnCategorical(simplex(1,1))
moves.append( mvRandomGeometricWalk(clock_indicator, weight=10.0, tune=FALSE) )

    source("rjMCMC/model_ACLN.Rev")
    source("rjMCMC/model_UCLN.Rev")

    clock_vec := v(acln_rates, branch_rates)
    clock_rates := clock_vec[clock_indicator]
    clock_indicator.setValue(2)



###########################################
###########################################
#SHDM Model
#source("rjMCMC/SHDM_8.Rev")
source("rjMCMC/SHDM_7.Rev")
source("rjMCMC/SHDM_6.Rev")
source("rjMCMC/SHDM_5.Rev")
source("rjMCMC/SHDM_4.Rev")
source("rjMCMC/SHDM_3.Rev")
source("rjMCMC/SHDM_2.Rev")

###########################################

#################################################
# Define the model of among-site rate variation #
#################################################

alpha ~ dnGamma(1E8, 0.5)
moves.append(mvScale(alpha, weight=10.0))
moves.append(mvScale(alpha, weight=10.0))


site_rates := fnDiscretizeGamma(alpha, alpha, 4)
###########################################

model_indicator_f ~ dnCategorical(simplex(1,1,1,1,1,1))
moves.append( mvRandomGeometricWalk(model_indicator_f, weight=10.0, tune=FALSE) )
model_indicator_f.setValue(1)

model_indicator_nf ~ dnCategorical(simplex(1,1,1,1,1,1))
moves.append( mvRandomGeometricWalk(model_indicator_nf, weight=10.0, tune=FALSE) )
model_indicator_nf.setValue(1)

morpho_f_bystate[1] <- morpho_f
morpho_nf_bystate[1] <- morpho_nf

n_max_states <- 4
i<- 1
for (i in 1:n_max_states) {
    morpho_f_bystate[i] <- morpho_f
    morpho_nf_bystate[i] <- morpho_nf

    morpho_f_bystate[i].setNumStatesPartition(i)
    morpho_nf_bystate[i].setNumStatesPartition(i)

    nc = morpho_f_bystate[i].nchar()
    if (nc > 0) {

    Q_vec_f := v(Q_SHDM_f7, Q_SHDM_f6, Q_SHDM_f5, Q_SHDM_f4, Q_SHDM_f3, Q_SHDM_f2)
    Q_vec_nf := v(Q_SHDM_nf7, Q_SHDM_nf6, Q_SHDM_nf5, Q_SHDM_nf4, Q_SHDM_nf3, Q_SHDM_nf2)

    Q_f := Q_vec_f[model_indicator_f]
    Q_nf := Q_vec_nf[model_indicator_nf]

    matrix_f_probs_vec := v(matrix_probs_f7, matrix_probs_f6, matrix_probs_f5, matrix_probs_f4, matrix_probs_f3, matrix_probs_f2)
    matrix_probs_vec_nf := v(matrix_probs_nf7, matrix_probs_nf6, matrix_probs_nf5, matrix_probs_nf4, matrix_probs_nf3, matrix_probs_nf2)

    matrix_probs_f := matrix_f_probs_vec[model_indicator_f]
    matrix_probs_nf := matrix_probs_vec_nf[model_indicator_nf]

    seq[i] ~ dnPhyloCTMC(tree=fbd_tree, Q = Q_f, type="Standard", siteRates=site_rates, branchRates = clock_rates, siteMatrices = matrix_probs_f)
    nf_seq[i] ~ dnPhyloCTMC(tree=fbd_tree, Q = Q_nf, type="Standard", siteRates=site_rates, branchRates = clock_rates, siteMatrices = matrix_probs_nf)

    seq[i].clamp(morpho_f_bystate[i])
    nf_seq[i].clamp(morpho_nf_bystate[i])

    print(i)
    i = i+1

    }
}


     mymodel = model(fbd_tree)



     monitors.append(mnModel(filename="output/cinc6_dated.log", printgen=10))



     monitors.append(mnFile(filename="output/cinc6_dated.trees", printgen=10, pruned_fbd_tree))



     monitors.append(mnScreen(printgen=10, num_samp_anc,fbd_indicator, clock_indicator, model_indicator_f, model_indicator_nf))



     mymcmc = mcmc(mymodel, monitors, moves)


     mymcmc.run(generations=1000000)



     q()
