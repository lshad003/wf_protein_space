# PROJECT_LOG.md - one line per chat, append-only
2026-08-21 | catalog build | rebuilt as one three-cohort catalog, 89 metagenomes, three tiers
2026-08-22 | annotation | eggNOG complete, 48.27% annotated, dark fraction 55.35% after filtering
2026-08-23 | taxonomy and Pfam | annotated fraction 92.6% bacterial; Pfam switched hmmscan to hmmsearch
2026-08-23 | [ANALYSIS] abundance groundwork + union staging | supported set, WF22 CDS, 1.75Gbp reference+index, reads_89 manifest join, WF24 code fix, pilot+union jobs queued
2026-08-23 | [ANALYSIS] union + rep classes + GU/EU screen | supported union-dark 51.42% (from 55.35%), rep_classes verified, sprot mislabel fixed, nr_cluster search launched
2026-08-24 | [ANALYSIS] abundance + dark fraction | 89/89 mapped, matrices built, U = 28.33% of genes but 19.28% of copies; U genes short (402bp) and 66.9% complete, so small proteins not fragments
2026-08-24 | [ANALYSIS] dark fraction characterized | U = 28.33% of genes, 19.28% of copies; 402bp and 66.9% complete so small proteins not fragments; AntiFam 0.25% spurious; 497,044 U genes in all 3 years carrying 91.18% of U reads; WF22 found to be 15 animals not 44, v1 p-values inflated
2026-08-24 | [ANALYSIS] four-way split and treatment null | EU 430,831 genes (20.82%) absent from NCBI nr at E<=1e-10 qcov/scov>=50; unknowns do not form families; EU contigs 91.5% unclassified; no Basidiobolus effect in 9 contrasts across 3 cohorts and 7 strains; LinDA replaced by limma-voom
