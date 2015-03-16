# diffReps #
| **Column Name** | **Explanation** |
|:----------------|:----------------|
| **Chrom** |Chromosome name|
| **Start** |1-based start coordinate|
| **End** |1-based end coordinate (inclusive)|
| **Length** |Length of the differential site in bps|
| **Treatment.cnt** |Normalized read counts of the treatment group (separated by semi-colon)|
| **Control.cnt** |Normalized read counts of the control group (separated by semi-colon)|
| **Treatment.avg** |Avg. count (normalized) of the treatment group|
| **Control.avg** |Avg. count (normalized) of the control group|
| **Treatment.enr** |Fold enrichment vs. input based on avg. count (normalized) for treatment group|
| **Control.enr** |Fold enrichment vs. input based on avg. count (normalized) for control group|
| **Event** |Direction of enrichment change using the control group as reference|
| **log2FC** |Log2 fold change|
| **pval** |P-value|
| **padj** |BH-adjusted p-value (FDR)|
| **winSta** |Start coordinate of the core window which is the window with most significant p-value among the site|
| **winEnd** |End coordinate of the core window|
| **winFC** |Log2 fold change of the core window|
| **winP** |P-value of the core window|
| **winQ** |Adjusted p-value of the core window|
| **GName** |Associated gene name|
| **TName** |Associated transcript name|
| **Strand** |Associated transcription strand|
| **TSS** |Associated transcript start site|
| **TES** |Associated transcript end site|
| **Feature** |Classification of the site's location|
| **D2TSS** |Distance from the site center to the associated transcript's TSS|

# findHotspots #
| **Column Name** | **Explanation** |
|:----------------|:----------------|
| **Chrom** |Chromosome name|
| **Start** |1-based start coordinate|
| **End** |1-based end coordinate (inclusive)|
| **Length** |Length of the hotspot in bps|
| **enrich** |Enrichment ratio (compared with background)|
| **pval** |P-value of enrichment|
| **padj** |Adjusted p-value|
| **nsite** |Number of differential sites in the hotspot|
| **Sites** |One-string description of the hotspots. The differential sites are separated by semi-colon. Each differential site is described by filename:(info columns separated by comma):line\_number.|
| **ntype** |Number of different types of chromatin marks in the hotspot|
| **MarkType** |Unique chromatin marks in the hotspot. In order for findHotspots to recognize mark names, you must make the diffReps output file names contain "standard" chromatin mark names. Such as, diffout.h3k4me3.txt, diffout.h3k9ac.txt, diffout.polII.txt. If it cannot recognize the marks, it will put the whole file names here.|