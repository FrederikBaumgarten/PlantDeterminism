===============================================================================
Classification of tree species for their degree of determinacy
===============================================================================

DOI of this archive:   10.5281/zenodo.18637216
License:               CC0-1.0 (public domain dedication)

-------------------------------------------------------------------------------
1. AUTHOR / CONTACT
-------------------------------------------------------------------------------
Frederik Martin Baumgarten
ORCID: 0000-0002-8284-8384
University of Basel, Department of Environmental Sciences
Contact: frederik.baumgarten@unibas.ch

Funding: Swiss National Science Foundation, postdoc mobility P500PB_210943.

-------------------------------------------------------------------------------
2. ASSOCIATED PUBLICATION
-------------------------------------------------------------------------------
This archive contains the data and code underlying the determinacy phylogeny
figure in:

Baumgarten et al. (2026). (In)determinacy in woody plants: limits and opportunities for timing growth in a changing climate. Ecology Letters. [DOI]


-------------------------------------------------------------------------------
3. PROJECT OVERVIEW
-------------------------------------------------------------------------------
The data record a qualitative classification of tree species by the degree of
determinacy of their apical shoot growth (determinate / intermediate /
indeterminate), compiled from the literature. The accompanying R script maps
this trait onto a pruned seed-plant phylogeny to examine whether indeterminate
growth is concentrated in earlier-diverging lineages.

Note: the compiled data are deliberately not exhaustive; they were assembled
to support a perspective paper and reflect the qualitative reporting of apical
shoot growth in the source literature.

-------------------------------------------------------------------------------
4. FILES IN THIS ARCHIVE
-------------------------------------------------------------------------------
  det_phylogeny_fig.R        R script: reads the data, encodes the determinacy
                             trait, prunes the megatree, and produces the
                             figure. Runs from the folder it sits in.

  Determinism_data.csv       The data table (UTF-8 CSV; see column dictionary
                             below and Determinism_metadata.csv). 105 rows.

  Determinism_metadata.csv   Column-by-column data dictionary for
                             Determinism_data.csv.

  phyloIntColor.pdf          The output figure (provided for reference; the
                             script regenerates it).

  read_me.txt                This file.

  Determinism.xlsx           (Optional) original spreadsheet, retained for
                             convenience. The CSV files are the archival
                             versions of record.

NOT INCLUDED — must be downloaded separately (see section 6):
  ALLMB.tre                  Smith & Brown (2018) seed-plant megatree.

-------------------------------------------------------------------------------
5. DATA DICTIONARY (Determinism_data.csv)
-------------------------------------------------------------------------------
Each row is one species record from one source study. A species may appear in
more than one row; the script averages the determinacy trait per species.
Blank cells denote "not applicable" or "not reported".

  paperID    Source study: first-author surname + year (e.g. "Halle1978").
  source     DOI or other identifier of the source study.
  genus      Genus name.
  species    Species epithet. "spp."/"ssp." marks unidentified or sub-specific
             entries; these rows are dropped from the analysis.
  var        Variety, if applicable.
  continent  Continent(s) of distribution; comma-separated if more than one.
  deciduousness  Leaf habit: D = deciduous, E = evergreen; "D, E" if both occur.
  distribution   Climatic region(s): e.g. temperate, subtropical, tropical,
             boreal, mediterranean; comma/slash-separated if more than one.
             A trailing "?" marks an uncertain assignment.
  Determinate_all_preformed        "x" = all leaves preformed in the previous
             year, no further foliage in the current season. -> trait value 3.
  Intermediate_preformed_neogrown  "x" = part of the canopy formed in the
             previous AND part in the current year.           -> trait value 2.
  Indeterminate_all_neogrown       "x" = (almost) no preformation; all leaves
             formed in the current season.                    -> trait value 1.
  sec_flush_y_no    "y" = a second mid-season flush occurs ("lammas growth" /
             "Johannistrieb"); blank = not reported.
  multiple_flushes  "y" = several flushes occur in the current season
             ("polycyclic flushing"); blank = not reported.
  remarks    Free-text remarks.

Determinacy trait encoding used in the analysis:
  1 = indeterminate, 2 = intermediate (mixed), 3 = determinate.

-------------------------------------------------------------------------------
6. EXTERNAL DATA REQUIRED: ALLMB.tre
-------------------------------------------------------------------------------
The script prunes the seed-plant megatree ALLMB.tre. This file is large and is
distributed under its own terms, so it is NOT bundled here. Download
"ALLMB.tre v1.0" and place it in the SAME folder as the script:

  https://doi.org/10.6084/m9.figshare.9747638

Cite:
  Smith, S.A. & Brown, J.W. (2018) Constructing a broadly inclusive seed plant
  phylogeny. American Journal of Botany 105, 302-314.
  https://doi.org/10.1002/ajb2.1019

-------------------------------------------------------------------------------
7. HOW TO REPRODUCE THE FIGURE
-------------------------------------------------------------------------------
  1. Put all archive files and ALLMB.tre in a single folder.
  2. Open a CLEAN R session (no pre-loaded objects).
  3. Source the script from that folder, e.g.:
        source("det_phylogeny_fig.R")
     or run it with Rscript from within the folder:
        Rscript det_phylogeny_fig.R
  4. Outputs are written to the same folder:
        determinacyPhylogeny.tre
        phyloIntColor.pdf

The script uses paths relative to its own location; no paths or working
directories need to be edited.

-------------------------------------------------------------------------------
8. SOFTWARE ENVIRONMENT (original analysis)
-------------------------------------------------------------------------------
  R version:    R version 4.3.1 (2023-06-16)
  Platform:     aarch64-apple-darwin20

  R packages used (record the versions actually used):
    stringr    1.5.1
    ape        5.8 
    phytools   2.3-0
    geiger     2.0.11
    pez        1.2-4
    caper      1.0.3
    phangorn   2.11.1


