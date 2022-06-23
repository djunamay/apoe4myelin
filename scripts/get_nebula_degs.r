library(nebula)


sce = readRDS('../data/single_cell_data/single_cell_experiment_object.rds')

# for all
oli = sce[,sce$cell.type=='Oli']
ace$fraction_mito = ace$MitoCounts/ace$TotalCounts
ace$APOE4 = ifelse(ace$apoe_genotype%in%c('34', '44', '24'), 1, 0)

model = ~ amyloid + nft + msex + age_death + pmi + fraction_mito + APOE4
counts = counts(ace)
meta = colData(ace)
mod = model.matrix(model, data=meta)

# run nebula
re = nebula(counts,meta$projid, pred=mod, offset=meta$total_counts, cpc=0.1, model = "NBGMM", kappa = 800, cutoff_cell = 20) # use default cpc?

# remove genes that dont pass the QC

# compare nebula to the pseudo bulk results and compare nebula to the wilcox results
