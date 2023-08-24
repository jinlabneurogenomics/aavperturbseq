#' Reduce redundancy of human GO terms
#'
#' @description This function will reduce GO redundancy first by creating a
#'  semantic similarity matrix (using
#'  \code{GOSemSim::\link[GOSemSim:mgoSim]{mgoSim}}), which is then passed
#'  through \code{rrvgo::\link[rrvgo:reduceSimMatrix]{reduceSimMatrix()}},
#'  which will reduce a set of GO terms based on their semantic similarity and
#'  scores (in this case, a default score based on set size is assigned.)
#'
#' @details Semantic similarity is calculated using the "Wang" method, a
#'   graph-based strategy to compute semantic similarity using the topology of
#'   the GO graph structure. \code{GOSemSim::\link[GOSemSim:mgoSim]{mgoSim}}
#'   does permit use of other measures (primarily information-content measures),
#'   but "Wang" is used as the default in GOSemSim (and was, thus, used as the
#'   default here). If you wish to use a different measure, please refer to the
#'   [GOSemSim
#'   documentation](https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html).
#'
#'  \code{rrvgo::\link[rrvgo:reduceSimMatrix]{reduceSimMatrix()}} creates a
#'  distance matrix, defined as (1-simMatrix). The terms are then hierarchically
#'  clustered using complete linkage (an agglomerative, or "bottom-up"
#'  clustering approach), and the tree is cut at the desired threshold. The term
#'  with the highest "score" is used to represent each group.
#'
#' @param pathway_df a `data.frame` or [tibble][tibble::tbl_df-class] object,
#'  with the following columns:
#'  \itemize{
#'  \item `go_type`: the sub-ontology the
#'  GO term relates to. Should be one of `c("BP", "CC", "MF")`.
#'  \item `go_id`: the gene ontology identifier (e.g. GO:0016209)
#'  }
#' @param orgdb `character()` vector, indicating name of the org.* Bioconductor
#'   package to be used
#' @param threshold `numeric()` vector. Similarity threshold (0-1) for
#'  \code{rrvgo::\link[rrvgo:reduceSimMatrix]{reduceSimMatrix()}}. Default
#'  option is 0.7. Some guidance:
#'  \itemize{
#'  \item For large term groupings, use `threshold = 0.9`
#'  \item For medium term groupings, use `threshold = 0.7`
#'  \item For small term groupings, use `threshold = 0.5`
#'  \item For tiny term groupings, use `threshold = 0.4`
#'  }
#' @param scores \strong{named} vector, with scores (weights) assigned to each
#'  term. Higher is better. Can be NULL (default, means no scores. In this case,
#'  a default score based on set size is assigned, thus favoring larger sets).
#'  Note: if you have p-values as scores, consider log-transforming them
#'  (`-log10(p)`).
#' @param measure `character()` vector, indicating method to be used to
#'   calculate semantic similarity measure. Must be one of the methods supported
#'   by GOSemSim: c("Resnik", "Lin", "Rel", "Jiang", "Wang"). Default is "Wang".
#'
#' @return a [tibble][tibble::tbl_df-class] object of pathway results, a
#'  "reduced" parent term to which pathways have been assigned. New columns:
#'  \itemize{
#'  \item `parent_id`: the GO ID of the parent term
#'  \item `parent_term`: a description of the GO ID
#'  \item `parent_sim_score`: the similarity score between the child GO term and
#'  its parent term
#'  }
#' @export
#'
#' @importFrom stats setNames
#' @importFrom tidyselect contains
#'
#' @family GO-related functions
#' @seealso \code{\link{go_plot}} for plotting the output of `go_reduce`,
#'   \code{GOSemSim::\link[GOSemSim:mgoSim]{mgoSim}} for calculation of semantic
#'   similarity and
#'   \code{rrvgo::\link[rrvgo:reduceSimMatrix]{reduceSimMatrix()}} for reduction
#'   of similarity matrix
#'
#' @references \itemize{
#'  \item Yu et al. (2010) GOSemSim: an R package for measuring semantic
#'  similarity among GO terms and gene products \emph{Bioinformatics} (Oxford,
#'  England), 26:7 976--978, April 2010.
#'  \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
#'  PMID: 20179076
#'  \item Yu (2021) Biomedical Knowledge Mining using GOSemSim and
#'  clusterProfiler
#'  \url{https://yulab-smu.top/biomedical-knowledge-mining-book/index.html}
#'  \item Sayols S (2020). rrvgo: a Bioconductor package to
#'  reduce and visualize Gene Ontology terms. \url{https://ssayols.github.io/rrvgo}
#'  }
#'
#'
#' @examples
#' file_path <-
#'     system.file(
#'         "testdata",
#'         "go_test_data.txt",
#'         package = "rutils",
#'         mustWork = TRUE
#'     )
#'
#' pathway_df <-
#'     readr::read_delim(file_path,
#'         delim = "\t"
#'     )
#'
#' go_reduce(
#'     pathway_df = pathway_df,
#'     orgdb = "org.Hs.eg.db",
#'     threshold = 0.9,
#'     scores = NULL,
#'     measure = "Wang"
#' )
go_reduce <- 
  function(
    pathway_df,
    orgdb = "org.Hs.eg.db",
    threshold = 0.7,
    scores = NULL,
    measure = "Wang"
    ) {
    if (!measure %in% c("Resnik", "Lin", "Rel", "Jiang", "Wang")) {
        stop('Chosen measure is not one of the recognised measures, c("Resnik", "Lin", "Rel", "Jiang", "Wang").')
    }

    # Based on chosen measure, set whether or not to compute information content
    if (measure == "Wang") {
        computeIC <- FALSE
    } else {
        computeIC <- TRUE
    }

    # Get ontology
    ont <-
        pathway_df %>%
        .[["go_type"]] %>%
        unique()

    terms.not.reduc = c()

    if (any(!ont %in% c("BP", "CC", "MF"))) {
        stop('Column go_type does not contain the recognised sub-ontologies, c("BP", "CC", "MF")')
    }

    go_similarity <-
        setNames(
            object =
                vector(
                    mode = "list",
                    length = length(ont)
                ),
            nm = ont
        )

    for (i in 1:length(ont)) {
        print(stringr::str_c("Reducing sub-ontology: ", ont[i]))

        hsGO <-
            GOSemSim::godata(
                OrgDb = orgdb,
                ont = ont[i],
                computeIC = computeIC
            )

        # Get unique ont terms from pathway_df
        terms <-
            pathway_df %>%
            dplyr::filter(.data$go_type == ont[i]) %>%
            .[["go_id"]] %>%
            unique()

	if (length(terms)<2) {
		terms.not.reduc = c(terms.not.reduc, terms)
		next
		#cur.row = which(pathway_df$go_id %in% terms)
		#go_similarity[[i]] <- df[cur.row,] %>%
		#tibble::as_tibble() %>%
		#dplyr::mutate(
                #	parent_id = .data$go_id,
	        #        parent_term = .data$Description,
        	#        parent_sim_score = NA
		#)
		#next
	}

        # Calculate semantic similarity
        sim <-
            GOSemSim::mgoSim(
                GO1 = terms,
                GO2 = terms,
                semData = hsGO,
                measure = measure,
                combine = NULL
            )

        # Reduce terms as based on scores or set size assigned
        go_similarity[[i]] <-
            rrvgo::reduceSimMatrix(
                simMatrix = sim,
                threshold = threshold,
                orgdb = orgdb,
                scores = scores
            ) %>%
            tibble::as_tibble() %>%
            dplyr::rename(
                parent_id = .data$parent,
                parent_term = .data$parentTerm,
                parent_sim_score = .data$parentSimScore
            )
    }

    go_similarity <- Filter(Negate(is.null), go_similarity)

    go_sim_df <-
        go_similarity %>%
        qdapTools::list_df2df(col1 = "go_type")

    pathway_go_sim_df <-
        pathway_df %>%
        dplyr::inner_join(
            go_sim_df %>%
                dplyr::select(.data$go_type,
                    go_id = .data$go,
                    contains("parent")
                ),
            by = c("go_type", "go_id")
        ) %>%
        dplyr::arrange(.data$go_type, .data$parent_id, -.data$parent_sim_score)

    if (length(terms.not.reduc)>0) {
	# add back terms not clustered
	df.not.reduc = pathway_df[which(pathway_df$go_id %in% terms.not.reduc),] %>%
			dplyr::mutate(
				parent_id = .data$go_id,
				parent_sim_score = NA,
				parent_term = .data$Description
			)
    
	pathway_go_sim_df = rbind(pathway_go_sim_df, df.not.reduc)
    }

    return(pathway_go_sim_df)
}
