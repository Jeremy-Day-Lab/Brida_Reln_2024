# Differentially Expressed Gene (DEG) Analyses

(LRT = likelihood-ratio test)

### DEGs differentiating among cell types

**DEG_cellType.r** = DEGs for each cell type

### DEGs differentiating Reln knockdown from lacZ control

**DEG-Reln_ctrl-sex.r** = LRT controlling for sex

Others (haven't followed up on)
- DEG-Reln_LRT.r = LRT, reducing by intercept
- DEG-Reln_Wald.r = Wald test, no reducing factor
- DEG-Reln_bySex-notusing.r = LRT comparing males vs. females
- DEG-Reln_ctrl-group.r = LRT controlling for FACS group
- DEG-Reln_ctrl-groupAndSex.r = LRT controlling for FACS group and sex
