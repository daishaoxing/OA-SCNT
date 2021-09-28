bcftools merge -O b -o merged_DEL_new.bcf  \
                08431GFP_DELnew.bcf \
                08431neg_DELnew.bcf \
                08431pos_DELnew.bcf \
                08431WT_DELnew.bcf \
                090202body_DELnew.bcf \
                NC1_DELnew.bcf \
                NC2_DELnew.bcf \
                NC3_DELnew.bcf \
                NC4_DELnew.bcf \
                NC5_DELnew.bcf \
                NC6_DELnew.bcf \
                NC7_DELnew.bcf \
                NC8_DELnew.bcf \
                NC9_DELnew.bcf \
                PC1_DELnew.bcf \
                PC2_DELnew.bcf \
                PC3_DELnew.bcf \
                PC4_DELnew.bcf \
                PC5_DELnew.bcf \
                PC6_DELnew.bcf \
                PC7_DELnew.bcf \
                PC8_DELnew.bcf \
                090202head_DELnew.bcf \
                cloning-JR_DELnew.bcf \
                cloning-PF_DELnew.bcf \
                cloning-W_DELnew.bcf
bcftools view merged_DEL_new.bcf >merged_DEL_new.vcf

bcftools merge -O b -o merged_BND_new.bcf  \
                08431GFP_BNDnew.bcf \
                08431neg_BNDnew.bcf \
                08431pos_BNDnew.bcf \
                08431WT_BNDnew.bcf \
                090202body_BNDnew.bcf \
                NC1_BNDnew.bcf \
                NC2_BNDnew.bcf \
                NC3_BNDnew.bcf \
                NC4_BNDnew.bcf \
                NC5_BNDnew.bcf \
                NC6_BNDnew.bcf \
                NC7_BNDnew.bcf \
                NC8_BNDnew.bcf \
                NC9_BNDnew.bcf \
                PC1_BNDnew.bcf \
                PC2_BNDnew.bcf \
                PC3_BNDnew.bcf \
                PC4_BNDnew.bcf \
                PC5_BNDnew.bcf \
                PC6_BNDnew.bcf \
                PC7_BNDnew.bcf \
                PC8_BNDnew.bcf \
                090202head_BNDnew.bcf \
                cloning-JR_BNDnew.bcf \
                cloning-PF_BNDnew.bcf \
                cloning-W_BNDnew.bcf
bcftools view merged_BND_new.bcf >merged_BND_new.vcf

bcftools merge -O b -o merged_DUP_new.bcf  \
                08431GFP_DUPnew.bcf \
                08431neg_DUPnew.bcf \
                08431pos_DUPnew.bcf \
                08431WT_DUPnew.bcf \
                090202body_DUPnew.bcf \
                NC1_DUPnew.bcf \
                NC2_DUPnew.bcf \
                NC3_DUPnew.bcf \
                NC4_DUPnew.bcf \
                NC5_DUPnew.bcf \
                NC6_DUPnew.bcf \
                NC7_DUPnew.bcf \
                NC8_DUPnew.bcf \
                NC9_DUPnew.bcf \
                PC1_DUPnew.bcf \
                PC2_DUPnew.bcf \
                PC3_DUPnew.bcf \
                PC4_DUPnew.bcf \
                PC5_DUPnew.bcf \
                PC6_DUPnew.bcf \
                PC7_DUPnew.bcf \
                PC8_DUPnew.bcf \
                090202head_DUPnew.bcf \
                cloning-JR_DUPnew.bcf \
                cloning-PF_DUPnew.bcf \
                cloning-W_DUPnew.bcf
bcftools view merged_DUP_new.bcf >merged_DUP_new.vcf

bcftools merge -O b -o merged_INS_new.bcf  \
                08431GFP_INSnew.bcf \
                08431neg_INSnew.bcf \
                08431pos_INSnew.bcf \
                08431WT_INSnew.bcf \
                090202body_INSnew.bcf \
                NC1_INSnew.bcf \
                NC2_INSnew.bcf \
                NC3_INSnew.bcf \
                NC4_INSnew.bcf \
                NC5_INSnew.bcf \
                NC6_INSnew.bcf \
                NC7_INSnew.bcf \
                NC8_INSnew.bcf \
                NC9_INSnew.bcf \
                PC1_INSnew.bcf \
                PC2_INSnew.bcf \
                PC3_INSnew.bcf \
                PC4_INSnew.bcf \
                PC5_INSnew.bcf \
                PC6_INSnew.bcf \
                PC7_INSnew.bcf \
                PC8_INSnew.bcf \
                090202head_INSnew.bcf \
                cloning-JR_INSnew.bcf \
                cloning-PF_INSnew.bcf \
                cloning-W_INSnew.bcf
bcftools view merged_INS_new.bcf >merged_INS_new.vcf

bcftools merge -O b -o merged_INV_new.bcf  \
                08431GFP_INVnew.bcf \
                08431neg_INVnew.bcf \
                08431pos_INVnew.bcf \
                08431WT_INVnew.bcf \
                090202body_INVnew.bcf \
                NC1_INVnew.bcf \
                NC2_INVnew.bcf \
                NC3_INVnew.bcf \
                NC4_INVnew.bcf \
                NC5_INVnew.bcf \
                NC6_INVnew.bcf \
                NC7_INVnew.bcf \
                NC8_INVnew.bcf \
                NC9_INVnew.bcf \
                PC1_INVnew.bcf \
                PC2_INVnew.bcf \
                PC3_INVnew.bcf \
                PC4_INVnew.bcf \
                PC5_INVnew.bcf \
                PC6_INVnew.bcf \
                PC7_INVnew.bcf \
                PC8_INVnew.bcf \
                090202head_INVnew.bcf \
                cloning-JR_INVnew.bcf \
                cloning-PF_INVnew.bcf \
                cloning-W_INVnew.bcf
bcftools view merged_INV_new.bcf >merged_INV_new.vcf
