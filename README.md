# Coester-2026

**Hypothalamic projection architecture links genetic risk for common obesity to Lmx1a/b neurons**

Bernd Coester, Charlotte Høy Kruse, Oliver Christensen, Mathias Mathiasen, Dylan Belmont-Rausch, Tune H. Pers

Novo Nordisk Foundation Center for Basic Metabolic Research, University of Copenhagen, Copenhagen, Denmark

Corresponding author:
Tune H. Pers
Novo Nordisk Foundation Center for Basic Metabolic Research, Faculty of Health and Medical Sciences, University of Copenhagen
Maersk Tower, Building: 07-8-78
Blegdamsvej 3B, 2200 Copenhagen, Denmark
Email: tune.pers@sund.ku.dk

**Abstract**
Genome-wide association studies implicate the brain in common obesity, yet the neuronal populations mediating this risk remain functionally unvalidated. We present a projection-resolved atlas of the mouse mediobasal hypothalamus, combining retrograde tracing with single-nucleus RNA sequencing of 193,989 nuclei and spatial transcriptomics of 289,888 cells to resolve 43 neuronal and 13 glial cell types, each anatomically localised and assigned predicted projections to eight brain regions. Integrating this atlas with obesity GWAS identified three enriched populations spanning the ventromedial, dorsomedial and posterior hypothalamus. The most enriched, Lmx1a/b-expressing neurons of the premammillary and posterior hypothalamus, projected to the amygdala and lateral hypothalamus confirmed by whole-brain light-sheet imaging. Silencing these neurons caused progressive weight gain, exacerbated by an energy-dense diet. Corresponding LMX1A-expressing human neurons independently carried BMI enrichment. This browsable atlas links hypothalamic transcriptional identity to long-range output and nominates Lmx1a/b neurons as a conserved regulator of body weight.

HTML files are plotly-rendered spatial coordinates of 24 brain sections, where a jitter is added to blend the individual planes. The rendered cells are centroids and downsampled by cell type.
- FULL contains all transferred single-cell celltypes with the same 0.4 score threshold used throughout the paper
- REGIONS contains the exact manual annotation of Paxinos regions, left-right merged
- PROJECTIONS contains the projection label for each cell based on spatial statistics within each brain section (those are calls for individual cells, not the cell-type based consensus described in the paper)

CODE directory is populated with all the processing code. This is NOT curated, streamlined, or making sense in any way. I recommend feeding this to an LLM if you want to get any information out of it.

Conda_envs has all the details on the environments used.

The 8-plex IPYNBs contain the code used to create the figures in the paper. Any deviation from these plots was manually done in Inkscape.
