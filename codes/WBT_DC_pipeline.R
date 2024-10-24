pacman::p_load(tidyverse, tidymodels, doParallel,GSVA,ggplot2,ROSE,pROC,vip)


# CD --------------------------------------------------------
cd_1016_deg <- readRDS("data/cd_1016_deg.rds")
cd_expr_matrix <- readRDS("data/cd_expr_list.rds")
cd_annot <- readRDS("data/cd_annot.rds")
cd_tpm_list <- readRDS("data/cd_tpm_list.rds")

#  WBT-DC pipeline
cd_gene_list <- list()
set.seed(1234)
for (i in 1:100)
  cd_gene_list[[paste0("up_feature_",i)]] <- cd_1016_deg$DEG_res$sig_up$gene_id[1:200][sample(140,sample(61:100,1))]
for (i in 1:100)
  cd_gene_list[[paste0("down_feature_",i)]] <- cd_1016_deg$DEG_res$sig_down$gene_id[1:200][sample(140,sample(61:100,1))]


## oversample

resample_cd_data <- list()
resample_cd_data$data_1016 <- get_resample_data(gene_expression_file = cd_expr_matrix$data_1016,
                                                annot_file = cd_annot$data_1016,seed = 5)

resample_cd_data$data_87 <- get_resample_data(gene_expression_file = cd_expr_matrix$data_87,
                                              annot_file = cd_annot$data_87,seed = 5)

resample_cd_data$data_235 <- get_resample_data(gene_expression_file = cd_expr_matrix$data_235,
                                               annot_file = cd_annot$data_array_235,seed = 5)

cd_gsva_score <- list()

cd_gsva_score$data_1016_over <- do_gsva_array(resample_cd_data$data_1016$expr_over,
                                              cd_gene_list,
                                              resample_cd_data$data_1016$annot_over,"CD","Control")

cd_gsva_score$data_87_over <- do_gsva_array(resample_cd_data$data_87$expr_over,
                                            cd_gene_list,
                                            resample_cd_data$data_87$annot_over,"CD","Control")

cd_gsva_score$data_235_over <- do_gsva_array(resample_cd_data$data_235$expr_over,
                                             cd_gene_list,
                                             resample_cd_data$data_235$annot_over,"CD","Control")


cd_ml_res <- list()

cd_ml_res$rf_over <- do_ml_rf(training_data = cd_gsva_score$data_1016_over,
                              testing_data_1 = cd_gsva_score$data_87_over,
                              testing_data_2 = cd_gsva_score$data_235_over,
                              preprocess_step = "NON",
                              disease_name = "CD",seed = 13)


# conventional methods

cd_deg_top100 <- c(cd_1016_deg$DEG_res$sig_up$gene_id[1:100],cd_1016_deg$DEG_res$sig_down$gene_id[1:100])

cd_1016_expr_top100 <- resample_cd_data$data_1016$expr_over[rownames(resample_cd_data$data_1016$expr_over) %in% cd_deg_top100,] %>% t %>% 
  as.data.frame %>% rownames_to_column("ID") %>% left_join(resample_cd_data$data_1016$annot_over) %>% 
  relocate(class,.after = ID) %>% mutate(class = factor(class))

cd_87_expr_top100 <- resample_cd_data$data_87$expr_over[rownames(resample_cd_data$data_87$expr_over) %in% cd_deg_top100,] %>% t %>% 
  as.data.frame %>% rownames_to_column("ID") %>% left_join(resample_cd_data$data_87$annot_over) %>% 
  relocate(class,.after = ID) %>% mutate(class = factor(class))

cd_235_expr_top100 <- resample_cd_data$data_235$expr_over[rownames(resample_cd_data$data_235$expr_over) %in% cd_deg_top100,] %>% t %>% 
  as.data.frame %>% rownames_to_column("ID") %>% left_join(resample_cd_data$data_235$annot_over) %>% 
  relocate(class,.after = ID) %>% mutate(class = factor(class))


cd_ml_res$rf_conventional <- do_ml_rf(training_data = cd_1016_expr_top100,
                                      testing_data_1 = cd_87_expr_top100,
                                      testing_data_2 = cd_235_expr_top100,
                                      preprocess_step = "NON",
                                      disease_name = "CD",seed = 123)



# ROC versus converntional
cd_auc <- list()
cd_auc$actural_87 <- ifelse(cd_ml_res$rf_over$rf_pred_1$class == "CD", 1, 0)
cd_auc$actural_235 <- ifelse(cd_ml_res$rf_over$rf_pred_2$class == "CD", 1, 0)
cd_auc$pre_new_87 <- cd_ml_res$rf_over$rf_pred_1$.pred_CD
cd_auc$pre_old_87 <- cd_ml_res$rf_conventional$rf_pred_1$.pred_CD
cd_auc$pre_new_235 <- cd_ml_res$rf_over$rf_pred_2$.pred_CD
cd_auc$pre_old_235 <- cd_ml_res$rf_conventional$rf_pred_2$.pred_CD

## plot CD cohort 87
pdf("results/roc_comparison_cd_85.pdf", width = 6, height = 6)
rocobj1 <- plot.roc(cd_auc$actural_87, cd_auc$pre_new_87,
                    main="Crohn's Disease", percent=TRUE, col="#DB9D85",identity = FALSE)
rocobj2 <- lines.roc(cd_auc$actural_87, cd_auc$pre_old_87 , percent=TRUE, col="black")
abline(a = 100, b = -1, col = "grey", lty = 2)
testobj <- roc.test(rocobj1, rocobj2,method = "delong")
text(50, 50, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))
text(30, 25, labels = paste("AUC =", round(auc(rocobj1), 2)), col = "#DB9D85",cex = 0.8)
text(30, 20, labels = paste("AUC =", round(auc(rocobj2), 2)), col = "black",cex = 0.8)
legend("bottomright", legend=c("our method", "conventional method"), col=c("#DB9D85", "black"), lwd=2)
dev.off()


### plot CD cohort 235
pdf("results/roc_comparison_cd_235.pdf", width = 6, height = 6)
rocobj1 <- plot.roc(cd_auc$actural_235, cd_auc$pre_new_235,
                    main="Crohn's Disease", percent=TRUE, col="#DB9D85")
rocobj2 <- lines.roc(cd_auc$actural_235, cd_auc$pre_old_235 , percent=TRUE, col="black")
abline(a = 100, b = -1, col = "grey", lty = 2)
testobj <- roc.test(rocobj1, rocobj2,method = "delong")
text(50, 50, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))
text(30, 25, labels = paste("AUC =", round(auc(rocobj1), 2)), col = "#DB9D85",cex = 0.8)
text(30, 20, labels = paste("AUC =", round(auc(rocobj2), 2)), col = "black",cex = 0.8)
legend("bottomright", legend=c("our method", "conventional method"), col=c("#DB9D85", "black"), lwd=2)
dev.off()


# UC --------------------------------------------------------
uc_1016_deg <- readRDS("data/uc_1016_deg.rds")
uc_expr_matrix <- readRDS("data/uc_expr_list.rds")
uc_annot <- readRDS("data/uc_annot.rds")
uc_tpm_list <- readRDS("data/uc_tpm_list.rds")

#  WBT-DC pipeline
uc_gene_list <- list()
set.seed(12345)
for (i in 1:100)
  uc_gene_list[[paste0("up_feature_",i)]] <- uc_1016_deg$DEG_res$sig_up$gene_id[1:200][sample(200,sample(71:100,1))]
for (i in 1:200)
  uc_gene_list[[paste0("down_feature_",i)]] <- uc_1016_deg$DEG_res$sig_down$gene_id[1:200][sample(200,sample(71:100,1))]


## oversample

resample_uc_data <- list()
resample_uc_data$data_1016 <- get_resample_data(gene_expression_file = uc_expr_matrix$data_1016,
                                                annot_file = uc_annot$data_1016,seed = 5)

resample_uc_data$data_87 <- get_resample_data(gene_expression_file = uc_expr_matrix$data_87,
                                              annot_file = uc_annot$data_87,seed = 5)

resample_uc_data$data_235 <- get_resample_data(gene_expression_file = uc_expr_matrix$data_235,
                                               annot_file = uc_annot$data_array_235,seed = 5)


uc_gsva_score <- list()

uc_gsva_score$data_1016_over <- do_gsva_array(resample_uc_data$data_1016$expr_over,
                                              uc_gene_list,
                                              resample_uc_data$data_1016$annot_over,"UC","Control")

uc_gsva_score$data_87_over <- do_gsva_array(resample_uc_data$data_87$expr_over,
                                            uc_gene_list,
                                            resample_uc_data$data_87$annot_over,"UC","Control")

uc_gsva_score$data_235_over <- do_gsva_array(resample_uc_data$data_235$expr_over,
                                             uc_gene_list,
                                             resample_uc_data$data_235$annot_over,"UC","Control")


uc_ml_res <- list()

uc_ml_res$rf_over <- do_ml_rf(training_data = uc_gsva_score$data_1016_over,
                              testing_data_1 = uc_gsva_score$data_87_over,
                              testing_data_2 = uc_gsva_score$data_235_over,
                              preprocess_step = "NON",
                              disease_name = "UC",seed = 123)


## conventional methods
uc_deg_top100 <- c(uc_1016_deg$DEG_res$sig_up$gene_id[1:100],uc_1016_deg$DEG_res$sig_down$gene_id[1:100])

uc_1016_expr_top100 <- resample_uc_data$data_1016$expr_over[rownames(resample_uc_data$data_1016$expr_over) %in% uc_deg_top100,] %>% t %>% 
  as.data.frame %>% rownames_to_column("ID") %>% left_join(resample_uc_data$data_1016$annot_over) %>% 
  relocate(class,.after = ID) %>% mutate(class = factor(class))

uc_87_expr_top100 <- resample_uc_data$data_87$expr_over[rownames(resample_uc_data$data_87$expr_over) %in% uc_deg_top100,] %>% t %>% 
  as.data.frame %>% rownames_to_column("ID") %>% left_join(resample_uc_data$data_87$annot_over) %>% 
  relocate(class,.after = ID) %>% mutate(class = factor(class))

uc_235_expr_top100 <- resample_uc_data$data_235$expr_over[rownames(resample_uc_data$data_235$expr_over) %in% uc_deg_top100,] %>% t %>% 
  as.data.frame %>% rownames_to_column("ID") %>% left_join(resample_uc_data$data_235$annot_over) %>% 
  relocate(class,.after = ID) %>% mutate(class = factor(class))


uc_ml_res$rf_conventional <- do_ml_rf(training_data = uc_1016_expr_top100 ,
                                      testing_data_1 = uc_87_expr_top100,
                                      testing_data_2 = uc_235_expr_top100,
                                      preprocess_step = "NON",
                                      disease_name = "UC",seed = 12)


## ROC vs conventional
uc_auc <- list()
uc_auc$actural_87 <- ifelse(uc_ml_res$rf_over$rf_pred_1$class == "UC", 1, 0)
uc_auc$actural_235 <- ifelse(uc_ml_res$rf_over$rf_pred_2$class == "UC", 1, 0)
uc_auc$pre_new_87 <- uc_ml_res$rf_over$rf_pred_1$.pred_UC
uc_auc$pre_old_87 <- uc_ml_res$rf_conventional$rf_pred_1$.pred_UC
uc_auc$pre_new_235 <- uc_ml_res$rf_over$rf_pred_2$.pred_UC
uc_auc$pre_old_235 <- uc_ml_res$rf_conventional$rf_pred_2$.pred_UC

## plot UC cohort 87
pdf("results/roc_comparison_uc_85.pdf", width = 6, height = 6)
rocobj1 <- plot.roc(uc_auc$actural_87, uc_auc$pre_new_87,
                    main="Ulcerative Colitis ", percent=TRUE, col="#9DB469",identity = FALSE)
rocobj2 <- lines.roc(uc_auc$actural_87, uc_auc$pre_old_87 , percent=TRUE, col="black")
abline(a = 100, b = -1, col = "grey", lty = 2)
testobj <- roc.test(rocobj1, rocobj2,method = "delong")

text(50, 50, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))
text(30, 25, labels = paste("AUC =", round(auc(rocobj1), 2)), col = "#9DB469",cex = 0.8)
text(30, 20, labels = paste("AUC =", round(auc(rocobj2), 2)), col = "black",cex = 0.8)
legend("bottomright", legend=c("our method", "conventional method"), col=c("#9DB469", "black"), lwd=2)
dev.off()


## plot UC cohort 235
pdf("results/roc_comparison_uc_235.pdf", width = 6, height = 6)
rocobj1 <- plot.roc(uc_auc$actural_235, uc_auc$pre_new_235,
                    main="Ulcerative Colitis", percent=TRUE, col="#9DB469",identity = FALSE)
rocobj2 <- lines.roc(uc_auc$actural_235, uc_auc$pre_old_235 , percent=TRUE, col="black")
abline(a = 100, b = -1, col = "grey", lty = 2)
testobj <- roc.test(rocobj1, rocobj2,method = "delong")
text(50, 50, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))
text(30, 25, labels = paste("AUC =", round(auc(rocobj1), 2)), col = "#9DB469",cex = 0.8)
text(30, 20, labels = paste("AUC =", round(auc(rocobj2), 2)), col = "black",cex = 0.8)
legend("bottomright", legend=c("our method", "conventional method"), col=c("#9DB469", "black"), lwd=2)
dev.off()


# ALS --------------------------------------------------------

als_741_deg <- readRDS("data/als_741_DEGs.rds")
als_expr_matrix <- readRDS("data/als_expr_list.rds")
als_annot <- readRDS("data/als_annot.rds")

#  WBT-DC pipeline
resample_als_data <- list()
resample_als_data$data_741 <- get_resample_data(gene_expression_file = als_expr_matrix$data_741,
                                                annot_file = als_annot$data_741,seed = 5)

resample_als_data$data_85 <- get_resample_data(gene_expression_file = als_expr_matrix$data_85,
                                               annot_file = als_annot$data_85,seed = 5)

als_gene_list <- list()
set.seed(1235)
for (i in 1:100)
  als_gene_list[[paste0("up_feature_",i)]] <- als_741_deg$up_gene[1:160][sample(160,sample(51:100,1))]
for (i in 1:100)
  als_gene_list[[paste0("down_feature_",i)]] <- als_741_deg$down_gene[1:160][sample(160,sample(51:100,1))]

als_gsva_score <- list()


als_gsva_score$data_741_over <- do_gsva_array(resample_als_data$data_741$expr_over,
                                              als_gene_list,
                                              resample_als_data$data_741$annot_over,"ALS","Control")

als_gsva_score$data_85_over <- do_gsva_array(resample_als_data$data_85$expr_over,
                                             als_gene_list,
                                             resample_als_data$data_85$annot_over,"ALS","Control")


als_ml_res <- list()

als_ml_res$rf_over <- do_ml_rf(training_data = als_gsva_score$data_741_over,
                               testing_data_1 = als_gsva_score$data_85_over,
                               testing_data_2 = NULL,
                               preprocess_step = "NON",
                               disease_name = "ALS",seed = 123)


## conventional methods
als_deg_top100 <- c(als_741_deg$up_gene[1:100],als_741_deg$down_gene[1:100])

als_741_expr_top100 <- als_expr_matrix$data_741[rownames(als_expr_matrix$data_741) %in% als_deg_top100,] %>% t %>%
  as.data.frame %>% rownames_to_column("ID") %>% left_join(als_annot$data_741 %>% dplyr::select(ID,class)) %>%
  relocate(class,.after = ID) %>% mutate(class = factor(class))

als_85_expr_top100 <- als_expr_matrix$data_85[rownames(als_expr_matrix$data_85) %in% als_deg_top100,] %>% t %>%
  as.data.frame %>% rownames_to_column("ID") %>% left_join(als_annot$data_85 %>% dplyr::select(ID,class)) %>%
  relocate(class,.after = ID) %>% mutate(class = factor(class))

als_ml_res$rf_convention <- do_ml_rf(training_data = als_741_expr_top100,
                                     testing_data_1 = als_85_expr_top100,
                                     testing_data_2 = NULL,
                                     preprocess_step = "NON",
                                     disease_name = "ALS",seed = 123)

## auc vs conventional methods

als_auc <- list()
als_auc$actural_85_1 <- ifelse(als_ml_res$rf_over$rf_pred_1$class == "ALS", 1, 0)
als_auc$pre_new_85 <- als_ml_res$rf_over$rf_pred_1$.pred_ALS
als_auc$actural_85_2 <- ifelse(als_ml_res$rf_convention$rf_pred_1$class == "ALS", 1, 0)
als_auc$pre_old_85 <- als_ml_res$rf_convention$rf_pred_1$.pred_ALS

### plot
pdf("results/roc_comparison_als_85.pdf", width = 6, height = 6)
rocobj1 <- plot.roc(als_auc$actural_85_1, als_auc$pre_new_85,
                    main="Amyotrophic Lateral Sclerosis ", percent=TRUE, col="#BB9FE0")
rocobj2 <- lines.roc(als_auc$actural_85_2, als_auc$pre_old_85 , percent=TRUE, col="black")
abline(a = 100, b = -1, col = "grey", lty = 2)
testobj <- roc.test(rocobj1, rocobj2,method = "delong")
text(50, 50, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))
text(30, 25, labels = paste("AUC (Our Method) =", round(auc(rocobj1), 2)), col = "#BB9FE0",cex = 0.8)
text(30, 20, labels = paste("AUC (Conventional Method) =", round(auc(rocobj2), 2)), col = "black",cex = 0.8)  
legend("bottomright", legend=c("our method", "conventional method"), col=c("#BB9FE0", "black"), lwd=2)
dev.off()



# RA --------------------------------------------------------
ra_sza_deg <- readRDS("data/ra_sza_deg.rds")
ra_expr_matrix <- readRDS("data/ra_expr_list.rds")
ra_annot <- readRDS("data/ra_annot.rds")

#  WBT-DC pipeline
ra_gene_list <- list()
set.seed(1234)
for (i in 1:100)
  ra_gene_list[[paste0("up_feature_",i)]] <- ra_sza_deg$DEG_res$sig_up$gene_id[1:200][sample(200,sample(61:100,1))]
for (i in 1:100)
  ra_gene_list[[paste0("down_feature_",i)]] <- ra_sza_deg$DEG_res$sig_down$gene_id[1:200][sample(200,sample(61:100,1))]


## oversample
resample_ra_data <- list()
resample_ra_data$sza <- get_resample_data(gene_expression_file = ra_expr_matrix$sza,
                                          annot_file = ra_annot$sza,seed = 5)

resample_ra_data$data_275 <- get_resample_data(gene_expression_file = ra_expr_matrix$data_275,
                                               annot_file = ra_annot$data_275,seed = 5)


ra_gsva_score <- list()
ra_gsva_score$data_sza_over <- do_gsva_array(resample_ra_data$sza$expr_over,
                                             ra_gene_list,
                                             resample_ra_data$sza$annot_over,"RA","Control")

ra_gsva_score$data_275_over <- do_gsva_array(resample_ra_data$data_275$expr_over,
                                             ra_gene_list,
                                             resample_ra_data$data_275$annot_over,"RA","Control")

ra_ml_res <- list()
ra_ml_res$rf_over <- do_ml_rf(training_data = ra_gsva_score$data_sza_over,
                              testing_data_1 = ra_gsva_score$data_275_over,
                              testing_data_2 = NULL,
                              preprocess_step = "NON",
                              disease_name = "RA",seed = 13)

## conventional methods
ra_deg_top200 <- c(ra_sza_deg$DEG_res$sig_up$gene_id[1:200],ra_sza_deg$DEG_res$sig_down$gene_id[1:200])

ra_sza_expr_top200 <- ra_expr_matrix$sza[rownames(ra_expr_matrix$sza) %in% ra_deg_top200,] %>% t %>%
  as.data.frame %>% rownames_to_column("ID") %>% left_join(ra_annot$sza %>% dplyr::select(ID,class)) %>%
  relocate(class,.after = ID) %>% mutate(class = factor(class))

ra_275_expr_top200 <- ra_expr_matrix$data_275[rownames(ra_expr_matrix$data_275) %in% ra_deg_top200,] %>% t %>%
  as.data.frame %>% rownames_to_column("ID") %>% left_join(ra_annot$data_275 %>% dplyr::select(ID,class)) %>%
  relocate(class,.after = ID) %>% mutate(class = factor(class))


ra_ml_res$rf_convention <- do_ml_rf(training_data = ra_sza_expr_top200,
                                    testing_data_1 = ra_275_expr_top200,
                                    testing_data_2 = NULL,
                                    preprocess_step = "NON",
                                    disease_name = "RA",seed = 123)

## plot

ra_auc <- list()
ra_auc$actural_275 <- ifelse(ra_ml_res$rf_over$rf_pred_1$class == "RA", 1, 0)
ra_auc$pre_275 <- ra_ml_res$rf_over$rf_pred_1$.pred_RA

pdf("results/result/roc_ra_275.pdf", width = 6, height = 6)
rocobj <- plot.roc(ra_auc$actural_275 , ra_auc$pre_275,
                   main = "Confidence intervals", 
                   percent=TRUE,
                   ci = TRUE,                  
                   print.auc = TRUE)           
ciobj <- ci.se(rocobj,                         
               specificities = seq(0, 100, 5)) 
plot(ciobj, type = "shape", col = "#1c61b6AA")    
plot(ci(rocobj, of = "thresholds", thresholds = "best")) 
dev.off()


