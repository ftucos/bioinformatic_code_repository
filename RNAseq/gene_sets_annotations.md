# Gene Sets Annotations

Useful annotated genesets for GSEA/ORA

### MSigDB

Main resource used for every annotation

https://www.gsea-msigdb.org/gsea/msigdb/

To retreive the annotation set within R you can use:

```R
library(msigdbr)
signature = msigdbr(species = "human")
```

or 

```R
library(msigdf)
msigdf::msigdf.human
# or 
msigdf::msigdf.mouse
```

there is some discrepancy within the two libraries, this is probabile due to the presence in the msigdf retreived annotation sets of the hortologs of mouse annotation sets

### Mitocarta for metabolic annotation

https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways

### YAP-TAZ signature

Hippo/YAP-TAZ pathway is difficoult to retreive from common GSEA (even absent in most annotation sets) because it's regulation is mainly post-translational. 

A custom signature of specific downstream targets has been defined:

> Wang, Yumeng et al. “Comprehensive Molecular Characterization of the Hippo Signaling Pathway in Cancer.” Cell reports vol. 25,5 (2018): 1304-1317.e5. doi:10.1016/j.celrep.2018.10.001

```R
YAPTAZ_targets_signature <- c("MYOF", "AMOTL2", "LATS2", "CTGF", "CYR61", "ANKRD1", "ASAP1", "AXL", "F3", "IGFBP3", "CRIM1", "FJX1", "FOXF2", "GADD45A", "CCDC80", "NT5E", "DOCK5", "PTPN14", "ARHGEF17", "NUAK2", "TGFB2", "RBMS3")
```

### Collection of pathway Gene Sets

http://baderlab.org/GeneSets

They keep updating genesets periodically from different resources, offering also an "Human_allpathways*.gmt" with all the resources grouped together. This resource is particularry convenient for KEGG annotation. Some years ago KEGG changed the redistribution policy for it's gene sets and now manual download from theri website si required. Barderlab does that (don't think it's legit).