# DEPENDENCIES ####
library(data.table)
library(ggplot2)
library(cowplot)

# DATA LOAD ####

experiments <- fread('/home/nlim/MDE/GemmaPaper/Output/EE_Export.TSV')
tags <- fread('/home/nlim/MDE/GemmaPaper/Output/EETag_Export.TSV')
ontologies <- fread('/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED.TSV')
ontologies.defs <- fread('/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED_DEF.TSV')

experiments <- experiments[!(ad.IsTroubled) & ee.IsPublic, ]
tags <- tags[ee.ID %in% experiments$ee.ID, ]

# TISSUES ####

# Generate the map...
mMap.tissues_cells <- getMap(experiments, tags, ontologies, ontologies.defs, c('UBERON', 'CL'))

# Count across every region...
tissues_cells.all <- mMap.tissues_cells %>% getCounts()
tissues_cells.all.aggregate <- tissues_cells.all %>% getAggregateCounts()
tissues_cells.all.top <- c('brain', 'hematopoietic cell', 'leukocyte', 'blood', 'neural cell', 'liver', 'embryo', 'lymphocyte', 'skeletal system', 'epithelial cell', 'lung', 'macrophage', 'glial cell', 'spleen', 'fibroblast', 'heart', 'muscle structure', 'skin of body', 'T cell', 'retina', 'intestine', 'breast', 'kidney')

# Limit the scope to only terms reachable from 'brain'.
mGraph <- getGraph(ontologies, c('UBERON', 'CL'))
tissues_cells.brain <- mMap.tissues_cells %>% getCounts(ontoScope = c('http://purl.obolibrary.org/obo/UBERON_0000955', # brain
                                                                      'http://purl.obolibrary.org/obo/CL_0002319'), # neural cell
                                                        mGraph = mGraph)
tissues_cells.brain.aggregate <- tissues_cells.brain %>% getAggregateCounts()
tissues_cells.brain.top <- c('medulla oblongata', 'pons', 'cerebellum', 'thalamic complex', 'hypothalamus', 'pituitary gland', "Ammon's horn", 'amygdala', 'basal ganglion', 'olfactory lobe', 'frontal lobe', 'parietal lobe', 'occipital lobe', 'temporal lobe', 'cingulate cortex', 'midbrain', 'substantia nigra', 'motor neuron', 'sensory neuron', 'interneuron', 'epithalamus', 'astrocyte', 'Schwann cell', 'microglial cell')
rm(mGraph)

write.csv(tissues_cells.all, 'output/counts_tissues.GLOBAL.csv', row.names = F)
write.csv(tissues_cells.brain, 'output/counts_tissues.BRAIN.csv', row.names = F)
write.table(tissues_cells.all.top, 'output/top_tissues.GLOBAL.csv', sep = ',', quote = F, col.names = F, row.names = F)
write.table(tissues_cells.brain.top, 'output/top_tissues.BRAIN.csv', sep = ',', quote = F, col.names = F, row.names = F)

# DISEASES ####

# Generate the map...
mMap.diseases <- getMap(experiments, tags, ontologies, ontologies.defs, 'DO')

# Count!
diseases.all <- mMap.diseases %>% getCounts()
diseases.all.aggregate <- diseases.all %>% getAggregateCounts()
diseases.all.top <- c('neurodegenerative disease', 'disease of mental health', 'leukemia', 'breast cancer', 'lung cancer', 'bone cancer', 'intestinal cancer', 'autoimmune hypersensitivity disease', "Alzheimer's disease", 'non-Hodgkin lymphoma', 'brain cancer', 'liver cancer', "Parkinson's disease", 'prostate cancer', 'diabetes mellitus', 'mood disorder', 'schizophrenia', 'arthritis', "Huntington's disease", 'multiple sclerosis', 'autism spectrum disorder', 'epilepsy', 'bipolar disorder', 'asthma', 'tuberculosis', 'influenza', 'pneumonia', 'Down syndrome')

write.csv(diseases.all, 'output/counts_diseases.GLOBAL.csv', row.names = F)
write.table(diseases.all.top, 'output/top_diseases.GLOBAL.csv', sep = ',', quote = F, col.names = F, row.names = F)

# BRAIN TISSUE x DISEASE ####

# Extract only rows that relate to a brain disease...
terms.brain <- names(subcomponent(graph_from_data_frame(ontologies[OntologyScope == 'DO', .(ChildNode_Long, ParentNode_Long)]),
                                  'http://purl.obolibrary.org/obo/DOID_863', # nervous system disease
                                  'in'))
mMap.brain_diseases <- mMap.diseases[parents %in% terms.brain]

# Extract only rows that relate to a brain tissue...
terms.brain <- unlist(lapply(c('http://purl.obolibrary.org/obo/UBERON_0000955', # brain
                               'http://purl.obolibrary.org/obo/CL_0002319'), # neural cell
                             function(term)
                               names(subcomponent(graph_from_data_frame(ontologies[OntologyScope %in% c('UBERON', 'CL'), .(ChildNode_Long, ParentNode_Long)]),
                                                  term, 'in'))))
mMap.brain_tissues <- mMap.tissues_cells[parents %in% terms.brain]
rm(terms.brain)

# Link brain tissues with brain diseases by experiment ID and count.
diseases.x.tissues <- merge(mMap.brain_diseases, mMap.brain_tissues, by = c('ee.ID', 'ee.Taxon')) %>% getCounts()
diseases.x.tissues.aggregate <- diseases.x.tissues %>% getAggregateCounts()
diseases.x.tissues.top <- c("Alzheimer's disease:hippocampal formation", "Alzheimer's disease:limbic lobe", "Huntington's disease:basal ganglion", 'epilepsy:hippocampal formation', "Alzheimer's disease:Ammon's horn", "Parkinson's disease:substantia nigra", "Alzheimer's disease:glial cell", 'epilepsy:limbic lobe', 'brain cancer:cerebellum', "epilepsy:Ammon's horn", "Parkinson's disease:dopaminergic neuron", "Alzheimer's disease:temporal lobe", 'amyotrophic lateral sclerosis:motor neuron')

write.csv(diseases.x.tissues, 'output/counts_intersect_diseases_tissues.BRAIN.csv', row.names = F)
write.table(diseases.x.tissues.top, 'output/top_intersect_diseases_tissues.BRAIN.csv', sep = ',', quote = F, col.names = F, row.names = F)

# DISEASE x DRUGS ####

# Generate the map...
mMap.drugs <- getMap(experiments, tags, ontologies, ontologies.defs, 'CHEBI')

# Link drugs with diseases by experiment ID and count.
diseases.x.drugs <- merge(mMap.drugs, mMap.diseases, by = c('ee.ID', 'ee.Taxon'), allow.cartesian = T) %>% getCounts()
diseases.x.drugs.aggregate <- diseases.x.drugs %>% getAggregateCounts()

write.csv(diseases.x.drugs, 'output/counts_intersect_diseases_drugs.csv', row.names = F)

# GENES ####

genes.all <- setorder(merge(tags[et.Cat == 'genotype', ee.ID, et.Val], experiments[, .(ee.ID, ee.Taxon)], by = 'ee.ID')[, .N, .(et.Val, ee.Taxon)], -N)
write.csv(genes.all, 'output/counts_genes.GLOBAL.csv', row.names = F)

# getOrderedDefinitions ####
# Yields the Definitions inside of a @arg{countMap} (@see{getCounts}) as a factor with levels sorted by
# total (aggregate, @see{getAggregateCounts}) counts.
getOrderedDefinitions <- function(countMap) {
  factor(countMap[, Definition], levels = setorder(getAggregateCounts(countMap)[, .(Definition, V1)], V1)$Definition)
}

# PLOT ####

theme_set(theme_cowplot())

tissues_cells.all[Definition %in% tissues_cells.all.top] %>%
  ggplot(aes(getOrderedDefinitions(.[]), N, fill = ee.Taxon)) +
  geom_bar(stat = 'identity') + coord_flip() + scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2500, 250)) + theme(legend.position = c(0.8, 0.2)) + ylab('# of Datasets') +
  ggtitle('Top Tissues/Cell Types') + xlab('Tissue or Cell Type') + labs(fill = 'Taxon'); ggsave2('output/figures/plot_tissues.GLOBAL.pdf', width = 11, height = 8.5)

diseases.all[Definition %in% diseases.all.top] %>%
  ggplot(aes(getOrderedDefinitions(.[]), N, fill = ee.Taxon)) +
  geom_bar(stat = 'identity') + coord_flip() + scale_y_continuous(expand = c(0, 0), breaks = seq(0, 400, 50)) + theme(legend.position = c(0.8, 0.2)) + ylab('# of Datasets') +
  ggtitle('Top Diseases') + xlab('Disease') + labs(fill = 'Taxon'); ggsave2('output/figures/plot_diseases.GLOBAL.pdf', width = 11, height = 8.5)

tissues_cells.brain[Definition %in% tissues_cells.brain.top] %>%
  ggplot(aes(getOrderedDefinitions(.[]), N, fill = ee.Taxon)) +
  geom_bar(stat = 'identity') + coord_flip() + scale_y_continuous(expand = c(0, 0), breaks = seq(0, 450, 50)) + theme(legend.position = c(0.8, 0.2)) + ylab('# of Datasets') +
  ggtitle('Top Brain Tissues/Cell Types') + xlab('Tissue or Cell Type') + labs(fill = 'Taxon'); ggsave2('output/figures/plot_tissues.BRAIN.pdf', width = 11, height = 8.5)

diseases.x.tissues[, Definition := paste(Definition.x, Definition.y, sep = ':')][Definition %in% diseases.x.tissues.top] %>%
  ggplot(aes(getOrderedDefinitions(.[]), N, fill = ee.Taxon)) +
  geom_bar(stat = 'identity') + coord_flip() + scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30, 5)) + theme(legend.position = c(0.8, 0.2)) + ylab('# of Datasets') +
  ggtitle('Top Disease/Tissue Overlaps') + xlab('Disease:Tissue') + labs(fill = 'Taxon'); ggsave2('output/figures/plot_intersect_diseases_tissues.BRAIN.pdf', width = 11, height = 8.5)
