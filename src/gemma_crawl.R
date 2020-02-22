# DEPENDENCIES ####
library(data.table)
library(igraph)
library(dplyr)

# getGraph ####
# Generates an igraph object from the given @arg{ontologies}, including only those in @arg{ontoScope}
getGraph <- function(ontologies, ontoScope = ontologies[, unique(OntologyScope)]) {
  suppressWarnings(graph_from_data_frame(ontologies[OntologyScope %in% ontoScope, .(ChildNode_Long, ParentNode_Long)]))
}

# getMap ####
# Generates a data table with a row for every parent term that can be inferred from the explicit terms in the experiment @arg{tags}.
# Associates taxon information from the paired @arg{experiments} and uses @arg{ontologies} (structure) and @arg{ontologies.defs} (definitions)
# to make inferrences, within the scope of @arg{ontoScope}.
getMap <- function(experiments, tags, ontologies, ontologies.defs, ontoScope = ontologies.defs[, unique(OntologyScope)]) {
  # Restrict focus to only terms in the desired ontological scope...
  mTags <- tags[!is.na(et.ValUri)] %>% .[et.ValUri %in% ontologies.defs[Node_Long %in% et.ValUri & OntologyScope %in% ontoScope, Node_Long]]
  
  # Generate the graph...
  mGraph <- getGraph(ontologies, ontoScope)
  
  # Generate the mapping table...
  do.call(rbind, lapply(1:nrow(mTags), function(indx) {
    mParents <- subcomponent(mGraph, mTags[indx, et.ValUri], 'out')
    data.table(ee.ID = mTags[indx, ee.ID], et.ID = mTags[indx, et.ID], parents = names(mParents))
  }))[parents != 'NA'] -> mMap
  
  # Merge in taxon information based on experiment IDs...
  mMap <- merge(mMap, experiments[, .(ee.ID, ee.Taxon)], by = 'ee.ID', all.x = T)
  
  # Merge in definitions...
  mMap <- merge(mMap, ontologies.defs[OntologyScope %in% ontoScope & Node_Long %in% mMap[, unique(parents)], unique(Definition), Node_Long],
                by.x = 'parents', by.y = 'Node_Long', all.x = T)[, c('Definition', 'V1') := list(V1, NULL)]
}

# getCounts ####
# Perform a taxon-level count of unique terms from a @arg{countMap} (@see{getMap}). Optionally restrict focus to
# the given @arg{experimentScope} (by ee.ID), @arg{taxaScope} or @arg{ontoScope}. If specifying an ontoScope, an
# igraph must be provided in @arg{mGraph} (@see{getGraph}).
getCounts <- function(countMap, experimentScope = countMap[, ee.ID], taxaScope = c('human', 'mouse', 'rat'), ontoScope = NULL, mGraph = NULL) {
  ## If we're doing interactions, we count a little differently.
  getLinkedCounts <- function(countMap, experimentScope, taxaScope) {
    if(!is.null(ontoScope)) warning('Counts were linked; ignoring ontoScope!')
    countMap[ee.ID %in% experimentScope & ee.Taxon %in% taxaScope, unique(data.table(Definition.x, Definition.y)), .(ee.Taxon, ee.ID)
             ][, .N, .(ee.Taxon, Definition.x, Definition.y) # Counts...
               ] %>% setorder(-N) # Ordering.
  }
  
  # Do linked counts if our count map looks like it's linked.
  if(ncol(countMap) == 8) return(getLinkedCounts(countMap, experimentScope, taxaScope))
  
  # If we're limiting our ontology scope, we only want to count tags that are reachable from the given vertices.
  if(!is.null(ontoScope)) {
    if(is.null(mGraph)) warning('mGraph was not supplied; ignoring ontoScope!')
    else countMap <- countMap[parents %in% unique(unlist(lapply(ontoScope, function(term) names(subcomponent(mGraph, term, 'in')))))]
  }
  
  # Finally, count!
  countMap[ee.ID %in% experimentScope & ee.Taxon %in% taxaScope, unique(Definition), .(ee.Taxon, ee.ID) # Taxon/experiment-level uniqueness of values...
           ][, .N, .(ee.Taxon, V1) # Counts...
             ][, c('Definition', 'V1') := list(V1, NULL)] %>% setorder(-N) # Sane names and ordering.
}

# getAggregateCounts ####
# Get taxon-independent counts from a previously computed @arg{countMap} (@see{getCounts}).
getAggregateCounts <- function(countMap) {
  if(ncol(countMap) == 4) countMap[, sum(N), .(Definition.x, Definition.y)] %>% setorder(-V1)
  else countMap[, sum(N), Definition] %>% setorder(-V1)
}

