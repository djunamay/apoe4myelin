print('making metadata file')
all_metadata = read.csv('../data/single_cell_data/qc_counts_data/qc_column_metadata.csv')
meta = all_metadata[!duplicated(all_metadata$projid),]
meta = meta[,c('projid', 'batch', 'amyloid', 'nft', 'msex', 'age_death','pmi', 'AD','apoe_genotype' ,'APOE4')]
write.csv(meta, '../data/single_cell_data/metadata_by_individual.csv')
print('done.')
