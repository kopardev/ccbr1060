## Scripts used and documentation for ccbr1060 data analysis.


**Data type:** Total RNASeq

**Details:**

*Cells:* 596-7 clone from HT-29 cellline ( Human colon adenicarcinoma cellline) 

*Groups:* 

Control: test1-test4 cells with shRNA construct but not activated with Dox

EFNB2 silenced: test5-test8 cells with shRNA construct activated with Dox


**TODOs:**
* Detect and identify location of shRNA integration sites
* Quantify reads representing shRNA integration
* Detect and quantify opposite strand transcription (OST = detectable expression of RNA belonging the strand opposite of annotated gene at the given genomic location) at some shRNA integration sites
* Why are you seeing EFNB2 transcripts even though no protein product is detected (lab-confirmed)? Exon-level analysis to confirm EFNB2 silencing
* Write generic scripts to detect OST in publicly available total RNAseq datasets