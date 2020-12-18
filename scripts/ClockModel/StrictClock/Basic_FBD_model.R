########################################################################
# Set up appropriate parameters for speciation, extinction & sampling  #
#   "Seed" numbers based on analyses of Paleobiology Database data.    #
########################################################################
# clock

mean_ra <- 7.0
stdv_ra <- 0.25
mu_ra <- ln(mean_ra) - ((stdv_ra*stdv_ra) * 0.5)

root_time ~ dnLognormal(mu_ra, stdv_ra, offset=7.3)

# Diversification Rates based on Echinodermata
speciation_rate ~ dnExponential(1.471);
# NOTE: If it gets stuck in this script, then set origination & extinction to 1.0
moves.append(mvScale(speciation_rate, lambda=0.01, weight=5));
moves.append(mvScale(speciation_rate, lambda=0.10, weight=3));
moves.append(mvScale(speciation_rate, lambda=1.00, weight=1));

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# NOTE: FBD scripts often allow extinction to vary independently of speciation;     #
# However, empirical studies show that these two rates usually are close to equal   #
#               and they definitely are not independent.                            #
# So, here we'll make turnover (ext/orig) an independent variable and use it        #
#               to scale extinction relative to origination                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
turnover ~ dnUnif(0.9, 1.05);
moves.append(mvSlide(turnover, delta=0.01, weight=5));
moves.append(mvSlide(turnover, delta=0.10, weight=3));
moves.append(mvSlide(turnover, delta=1.00, weight=1));
extinction_rate := turnover*speciation_rate;
diversification := speciation_rate - extinction_rate;

# old extinction stuff. We should not use this, as extinction should not be independent of origination!
#extinction_rate ~ dnExponential(1.471);
#moves.append(mvScale(extinction_rate, lambda=0.01, weight=5));
#moves.append(mvScale(extinction_rate, lambda=0.10, weight=3));
#moves.append(mvScale(extinction_rate, lambda=1.00, weight=1));
#turnover := extinction_rate/speciation_rate;

# Fossil Sampling Rates based on collection occupied by Echinodermata
psi ~ dnExponential(3.892);
completeness := psi/(extinction_rate+psi);
moves.append(mvScale(psi, lambda=0.01, weight=5));
moves.append(mvScale(psi, lambda=0.10, weight=3));
moves.append(mvScale(psi, lambda=1.00, weight=1));

# Proportional Taxon Sampling of Youngest Time Slice
rho <- 0.506;	# 'extant' sampling.

# Establish Basal Divergence Time
origin_time ~ dnUnif(7.3, 12.11);
moves.append(mvSlide(origin_time, delta=0.01, weight=5));
moves.append(mvSlide(origin_time, delta=0.10, weight=3));
moves.append(mvSlide(origin_time, delta=1.00, weight=1));

fbd_dist = dnFBDRP(origin=origin_time, lambda=speciation_rate, mu=extinction_rate, psi=psi, rho=rho, taxa=taxa);

############################################################################
#                               Set up tree                                #
############################################################################
# create the vector of clade constraints
outgroup = clade("Ctenocystis");
ingroup = clade("Gyrocystis_platessa","Gyrocystis_testudiformis","Gyrocystis_cruzae","Gyrocystis_badulesiensis","Gyrocystis_erecta","Progyrocystis_disjuncta","Protocinctus_mansillaensis","Elliptocinctus_barrandei","Elliptocinctus_vizcainoi","Sucocystis_theronensis","Sucocystis_bretoni","Lignanicystis_barriosensis","Undatacinctus_undata","Sucocystis_acrofera","Undatacinctus_quadricornuta","Undatacinctus_melendezi","Asturicystis_jaekeli","Sotocinctus_ubaghsi","Trochocystites_bohemicus","Trochocystoides_parvus","Ludwigicinctus_truncatus","Graciacystis_ambigua","Asturicystis_havliceki","Nelegerocystis_ivantzovi","Rozanovicystis_triangularis","Davidocinctus_pembrokensis");
unscored_taxa <- v(24,25,26,27);
constraints = v(ingroup);
fbd_tree ~ dnConstrainedTopology(fbd_dist,constraints=constraints);

moves.append(mvFNPR(fbd_tree , weight=num_branches/2));                              # time-tree pruning & grafting
moves.append(mvNNI(fbd_tree , weight=num_branches/2));                               # nearest-neighbor interchanges
moves.append(mvCollapseExpandFossilBranch(fbd_tree ,origin_time,weight=num_taxa/4)); # consider ancestor-descendant rather than sister species
moves.append(mvNodeTimeSlideUniform(fbd_tree , weight=num_branches/2));              # adjust divergence times
moves.append(mvRootTimeSlideUniform(fbd_tree , origin_time, weight=5));            # adjust basal divergence time.

num_samp_anc := fbd_tree.numSampledAncestors();
for (bn in 1:num_branches)	{
	divergence_dates[bn]:=fbd_tree.nodeAge(bn)                   # this is when a hypothesized ancestor diverges or an OTU is first seen;
	branch_lengths[bn]:=fbd_tree.branchLength(bn);               # this is branch *duration* not expected change!
	origin_dates[bn]:=fbd_tree.branchLength(bn)+fbd_tree.nodeAge(bn); # this is when a lineage diverged from its ancestor
	}
summed_gaps := sum(branch_lengths);


intervals = readDataDelimitedFile(file="data/cincta_fossil_intervals_FA.tsv", header=true)
# Setup the fossil tip sampling #
# Use a for loop to create a uniform distribution on the occurence time for each fossil #
# The boundaries of the uniform distribution are specified in the tsv file #
for(i in 1:intervals.size())
{
    taxon  = intervals[i][1]
    a_i = intervals[i][2]
    b_i = intervals[i][3]
    t[i] := tmrca(fbd_tree, clade(taxon))

    fossil[i] <- a_i
#    fossil[i] ~ dnSoftBoundUniformNormal(t[i] - b_i, t[i] - a_i, sd = 2, p = 0.025)
#    fossil[i].clamp(0)
}

clade_extant = clade("Sucocystis_acrofera");
#age_extant := tmrca(fbd_tree, clade_extant);	# There is no particularly good reason to keep this!

pruned_fbd_tree := fnPruneTree(fbd_tree, prune=v("Asturicystis_havliceki","Nelegerocystis_ivantzovi","Rozanovicystis_triangularis","Davidocinctus_pembrokensis"))
#### for easy printing on screen ####
fbd_p := speciation_rate;
fbd_q := extinction_rate;
fbd_r := psi;
