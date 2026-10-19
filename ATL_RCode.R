################################################################################
##  ATLANTA CHAPTER — FINAL MODELS AND FIGURES
##  
##
##  OUTPUT  (each figure is written as .pdf for submission and .tiff for Word)
##    Figure2-2_TreeSize_Community.tiff
##    Figure2-3_TreeSize_Species.tiff
##    Figure2-4_PestMetrics_Community.tiff
##    Figure2-5_SiteInteractions.tiff
##    Figure2-6_EcosystemServices.tiff
##    (Figure 2-1 is the study map and is not produced here.)
##
##  PACKAGES TO CITE: glmmTMB (Brooks et al. 2017), ordinal::clmm
##    (Christensen 2023), car::Anova (Fox & Weisberg 2019), emmeans (Lenth 2024),
##    DHARMa (Hartig 2024), performance (Lüdecke 2021), ggplot2 (Wickham 2016).
################################################################################


#### PACKAGES ####
library(tidyverse)
library(glmmTMB)       # Gaussian / negative-binomial / beta GLMMs
library(ordinal)       # clmm() for ordinal responses (Condition, HerbRate)
library(car)           # Anova() Type II omnibus tests
library(emmeans)       # marginal means and slopes
library(DHARMa)        # simulated residual diagnostics (glmmTMB only)
library(performance)   # r2(), check_collinearity()
library(patchwork)     # figure assembly
library(scales)        # axis formatting


#### DATA IMPORT ####
d2024 <- read_csv("Atlanta2024.csv", show_col_types = FALSE)
d2025 <- read_csv("Atlanta2025.csv", show_col_types = FALSE)

rd <- bind_rows(
  d2024 %>% mutate(Year = 2024),
  d2025 %>% mutate(Year = 2025)
)


#### VARIABLE TYPES ####
rd$TreeID  <- factor(rd$TreeID)
rd$Year    <- factor(rd$Year, levels = c("2024", "2025"))
rd$DCI     <- factor(rd$DCI,  levels = c("Zone_One", "Zone_Five"))  # ref = Zone 1
rd$Species <- factor(rd$Species)                                    # ref = Acer

rd$Impervious25   <- as.numeric(rd$Impervious25)
rd$SoilMoisture   <- as.numeric(rd$SoilMoisture)
rd$SoilCompaction <- as.numeric(rd$SoilCompaction)

rd$Height       <- as.numeric(rd$Height)
rd$DBH          <- as.numeric(rd$DBH)
rd$CanopyVolume <- as.numeric(rd$CanopyVolume)
rd$Condition    <- factor(rd$Condition, levels = 1:5, ordered = TRUE)
rd$HerbRate     <- factor(rd$HerbRate,  levels = 0:5, ordered = TRUE)
rd$AphidCount   <- as.integer(rd$AphidCount)
rd$ScaleCount   <- as.integer(rd$ScaleCount)
rd$SpeciesRichness <- as.integer(rd$SpeciesRichness)   # arthropod FAMILY richness
rd$Simpson_D    <- as.numeric(rd$Simpson_D)


#### PREDICTOR STANDARDIZATION ####

std2 <- function(x) as.numeric((x - mean(x, na.rm = TRUE)) / (2 * sd(x, na.rm = TRUE)))

rd$Imperv25_s       <- std2(rd$Impervious25)
rd$SoilMoisture_s   <- std2(rd$SoilMoisture)
rd$SoilCompaction_s <- std2(rd$SoilCompaction)


#### 2025-ONLY SUBSET ####
# DBH, canopy volume, arthropod richness, Simpson diversity and the ecosystem
# services were measured in 2025 only, so those models carry no Year term and no
# random intercept.
rd_2025 <- rd %>% filter(Year == "2025")


#### PEST SUBSET — Magnolia excluded (aphid and scale models ONLY) ####

rd_pest <- droplevels(subset(rd, Species != "Magnolia"))


################################################################################
##  REPORTING HELPERS
################################################################################

## Residual diagnostics. DHARMa does not support clmm, so ordinal models get a
## convergence check instead.
run_diag <- function(model, label = "") {
  cat("\n========== DIAGNOSTICS:", label, "==========\n")
  if (inherits(model, "glmmTMB")) {
    res <- simulateResiduals(model, plot = TRUE)
    print(testDispersion(res))
    # High VIFs are expected in an interaction model and are not a real problem.
    print(tryCatch(check_collinearity(model),
                   error = function(e) "collinearity check unavailable"))
  } else {
    cat("clmm object: DHARMa not supported.\n")
    cat("Convergence code (0 = OK):", model$optRes$convergence, "\n")
    cat("Max |gradient| (want < 1e-3):",
        tryCatch(max(abs(model$gradient)), error = function(e) NA), "\n")
  }
}

## Type II omnibus tests, plus zone contrasts within each species and pooled.
## If DCI:Species is significant, lead with the within-species contrasts; if not,
## lead with the pooled effect.
dci_effects <- function(model, label = "", ordinal = FALSE) {
  av <- car::Anova(model, type = "II")
  cat("\n---------- OMNIBUS Anova (type II):", label, "----------\n")
  print(av)

  rn <- rownames(av)
  int_row <- rn[grepl("DCI", rn) & grepl("Species", rn)]
  pcol <- grep("^Pr", colnames(av), value = TRUE)[1]
  if (length(int_row)) cat(sprintf("\nDCI:Species interaction: p = %.4f\n", av[int_row, pcol]))

  scale_arg <- if (ordinal) "latent" else "response"
  cat("\n-- Zone contrast WITHIN each species --\n")
  print(summary(pairs(emmeans(model, ~ DCI | Species, mode = scale_arg,
                              type = scale_arg), reverse = TRUE), infer = TRUE))
  cat("\n-- Pooled DCI effect --\n")
  print(summary(pairs(emmeans(model, ~ DCI, mode = scale_arg,
                              type = scale_arg), reverse = TRUE), infer = TRUE))
  if (ordinal) cat("(ordinal: exp(estimate) = odds ratio)\n")
  invisible(av)
}

## Species-specific slopes for each local site variable, on the link scale.
## Do not back-transform these.
species_slopes <- function(model, label = "") {
  cat("\n========== SPECIES-SPECIFIC SLOPES:", label, "==========\n")
  for (v in c("Imperv25_s", "SoilMoisture_s", "SoilCompaction_s")) {
    cat("\n--- slope of", v, "within each species (link scale) ---\n")
    print(summary(emtrends(model, ~ Species, var = v), infer = TRUE))
  }
}

## Impervious slope per DCI zone — the original DCI x impervious question.
dci_imperv_slope <- function(model, label = "") {
  cat("\n---------- IMPERVIOUS SLOPE by DCI zone:", label, "----------\n")
  et <- emtrends(model, ~ DCI, var = "Imperv25_s")
  print(summary(et, infer = TRUE))
  print(pairs(et))   # difference in slopes between zones = the DCI:Imperv test
  invisible(et)
}


################################################################################
##  TREE PERFORMANCE AND CONDITION
################################################################################

#### Height — repeated measures, Gaussian on log scale ####
m_Height <- glmmTMB(
  log(Height) ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species +
    Year + (1 | TreeID),
  data = rd, family = gaussian()
)

run_diag(m_Height, "Height")
dci_effects(m_Height, "Height")
species_slopes(m_Height, "Height")
dci_imperv_slope(m_Height, "Height")
r2(m_Height)

## Effect size for the pooled height sentence ("~X% shorter in Zone 5"):
confint(pairs(emmeans(m_Height, ~ DCI, type = "response"), reverse = TRUE))


#### DBH — 2025 only, Gaussian on log scale, no random effect ####
m_DBH <- glmmTMB(
  log(DBH) ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species,
  data = rd_2025, family = gaussian()
)

run_diag(m_DBH, "DBH")
dci_effects(m_DBH, "DBH")
species_slopes(m_DBH, "DBH")
dci_imperv_slope(m_DBH, "DBH")
r2(m_DBH)


#### Canopy volume — 2025 only, Gaussian on log scale ####
m_CV <- glmmTMB(
  log(CanopyVolume) ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species,
  data = rd_2025, family = gaussian()
)

run_diag(m_CV, "Canopy volume")
dci_effects(m_CV, "Canopy volume")
species_slopes(m_CV, "Canopy volume")
dci_imperv_slope(m_CV, "Canopy volume")
r2(m_CV)


#### Condition — ordinal (1-5), repeated measures, clmm ####
# glmmTMB has no ordinal family. Watch convergence: this is the model most likely
# to struggle. If it fails, drop the soil interactions first, then the impervious
# interaction, and state the limitation.
m_Cond <- clmm(
  Condition ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species +
    Year + (1 | TreeID),
  data = rd, link = "logit", Hess = TRUE,
  control = clmm.control(maxIter = 1000, gradTol = 1e-5)
)

run_diag(m_Cond, "Condition")
dci_effects(m_Cond, "Condition", ordinal = TRUE)


################################################################################
##  INSECT PESTS
################################################################################

#### Aphid density — repeated measures, negative binomial ####
# data = rd_pest (Magnolia dropped). Species term has 3 df, not 4.
m_Aphid <- glmmTMB(
  AphidCount ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species +
    Year + (1 | TreeID),
  data = rd_pest, family = nbinom2(link = "log")
)

run_diag(m_Aphid, "Aphid density")
dci_effects(m_Aphid, "Aphid density")
species_slopes(m_Aphid, "Aphid density")
dci_imperv_slope(m_Aphid, "Aphid density")
r2(m_Aphid)
confint(pairs(emmeans(m_Aphid, ~ DCI, type = "response"), reverse = TRUE))


#### Scale density — repeated measures, negative binomial ####

#   m_Scale_noOutlier <- update(m_Scale, data = subset(rd_pest, TreeID != "26"))
m_Scale <- glmmTMB(
  ScaleCount ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species +
    Year + (1 | TreeID),
  data = rd_pest, family = nbinom2(link = "log")
)

run_diag(m_Scale, "Scale density")
dci_effects(m_Scale, "Scale density")
species_slopes(m_Scale, "Scale density")
dci_imperv_slope(m_Scale, "Scale density")
r2(m_Scale)
confint(pairs(emmeans(m_Scale, ~ DCI, type = "response"), reverse = TRUE))


#### Herbivory rate — ordinal (0-5), repeated measures, clmm ####
m_Herb <- clmm(
  HerbRate ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species +
    Year + (1 | TreeID),
  data = rd, link = "logit", Hess = TRUE,
  control = clmm.control(maxIter = 1000, gradTol = 1e-5)
)

run_diag(m_Herb, "Herbivory rate")
dci_effects(m_Herb, "Herbivory rate", ordinal = TRUE)


################################################################################
##  ARTHROPOD COMMUNITY  (2025 only, no random effect)
################################################################################

#### Arthropod family richness — negative binomial ####
m_Rich <- glmmTMB(
  SpeciesRichness ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species,
  data = rd_2025, family = nbinom2(link = "log")
)

testZeroInflation(simulateResiduals(m_Rich))   # add ziformula = ~1 if needed
run_diag(m_Rich, "Arthropod family richness")
dci_effects(m_Rich, "Arthropod family richness")
species_slopes(m_Rich, "Arthropod family richness")
r2(m_Rich)


#### Simpson diversity — beta regression ####
# beta_family needs values strictly inside (0,1). If any equal 0 or 1, squeeze
# with Smithson & Verkuilen: (y*(n-1) + 0.5)/n
m_Simp <- glmmTMB(
  Simpson_D ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species,
  data = rd_2025, family = beta_family(link = "logit")
)

run_diag(m_Simp, "Simpson diversity")
dci_effects(m_Simp, "Simpson diversity")
species_slopes(m_Simp, "Simpson diversity")
r2(m_Simp)

## The within-species contrasts printed above are the ones reported in Results
## 3.4. Take them from HERE, not from any distance-adjusted run.


################################################################################
##  ECOSYSTEM SERVICES (i-Tree Eco) — 2025 only, Gaussian on log
################################################################################
##  >>> ROW-ALIGNMENT: READ THIS <<<
##  The i-Tree exports carry their own "Tree ID" running 1-120, which is NOT the
##  field TreeID (non-contiguous, runs to 182). The files can only be joined BY
##  ROW POSITION, and the row order is Atlanta2025.csv's original file order,
##  which rd_2025 preserves.
##
##  itree_row MUST be assigned before any arrange(). Sorting first silently
##  misaligns the join across species. The check below compares the i-Tree
##  species name against the field species row for row and halts on any
##  disagreement.
################################################################################

species_order <- c("Acer", "Lagerstroemia", "Liriodendron", "Magnolia", "Quercus")

## i-Tree common name -> field genus code (used only to validate the join)
.itree_sp_map <- c(
  "Red maple"          = "Acer",
  "Common crapemyrtle" = "Lagerstroemia",
  "Tulip tree"         = "Liriodendron",
  "Sweetbay"           = "Magnolia",
  "White oak"          = "Quercus"
)

## Water Intercepted contains grouped values ("2,068.60"), which readr may guess
## as character. This is a no-op on already-numeric columns.
num_clean <- function(x) as.numeric(gsub(",", "", as.character(x)))

## Drop the trailing blank row FIRST, then number rows in native export order.
read_itree <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(`Tree ID`)) %>%
    mutate(itree_row     = row_number(),
           itree_species = unname(.itree_sp_map[trimws(`Species Name`)]))
}

carbon <- read_itree("CarbonStorage.csv") %>%
  mutate(CarbonStorage = num_clean(`Carbon Storage (kg)`)) %>%
  dplyr::select(itree_row, itree_species, CarbonStorage)

pollution <- read_itree("PollutionRemoval.csv") %>%
  mutate(PM25_removal = num_clean(`PM2.5 Removed (g/yr)`)) %>%
  dplyr::select(itree_row, itree_species, PM25_removal)

hydrology <- read_itree("HydrologyEffects.csv") %>%
  mutate(Transpiration    = num_clean(`Transpiration (l/yr)`),
         WaterIntercepted = num_clean(`Water Intercepted (l/yr)`)) %>%
  dplyr::select(itree_row, itree_species, Transpiration, WaterIntercepted)

## Base frame: number the rows BEFORE any arrange().
es_base <- rd_2025 %>%
  mutate(Species   = factor(Species, levels = species_order),
         itree_row = row_number())

check_itree_alignment <- function(es, itree_df, file_label) {
  cmp <- es %>%
    dplyr::select(itree_row, TreeID, Species) %>%
    left_join(dplyr::select(itree_df, itree_row, itree_species), by = "itree_row")

  bad <- cmp %>%
    filter(is.na(itree_species) | as.character(Species) != itree_species)

  if (nrow(bad) > 0) {
    message("\n!! MISALIGNED ROWS in ", file_label, " (showing up to 10):")
    print(as.data.frame(head(bad, 10)))
    stop(sprintf(
      paste0("JOIN MISALIGNED: %d of %d rows in %s disagree with the field ",
             "species. Do NOT interpret these models. Check that the i-Tree ",
             "export is still in Atlanta2025.csv's original row order."),
      nrow(bad), nrow(cmp), file_label), call. = FALSE)
  }

  message(sprintf("  OK  %-22s %3d/%3d rows species-matched",
                  file_label, nrow(cmp), nrow(cmp)))
  data.frame(File = file_label, Rows_checked = nrow(cmp),
             Species_matched = nrow(cmp) - nrow(bad), Status = "PASS",
             stringsAsFactors = FALSE)
}

message("\n---- i-Tree join integrity check ----")
itree_join_check <- bind_rows(
  check_itree_alignment(es_base, carbon,    "CarbonStorage.csv"),
  check_itree_alignment(es_base, pollution, "PollutionRemoval.csv"),
  check_itree_alignment(es_base, hydrology, "HydrologyEffects.csv")
)

es_data <- es_base %>%
  left_join(dplyr::select(carbon,    itree_row, CarbonStorage),    by = "itree_row") %>%
  left_join(dplyr::select(pollution, itree_row, PM25_removal),     by = "itree_row") %>%
  left_join(dplyr::select(hydrology, itree_row,
                          Transpiration, WaterIntercepted),        by = "itree_row") %>%
  arrange(Species, TreeID)     # cosmetic only; itree_row is already attached


#### PM2.5 removal — Gaussian on log, +0.05 constant for the six zeros ####
m_PM25 <- glmmTMB(
  log(PM25_removal + 0.05) ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species,
  data = es_data, family = gaussian()
)

run_diag(m_PM25, "PM2.5 removal")
dci_effects(m_PM25, "PM2.5 removal")
species_slopes(m_PM25, "PM2.5 removal")
dci_imperv_slope(m_PM25, "PM2.5 removal")
r2(m_PM25)


#### Carbon storage — Gaussian on log, no constant (min 0.1 kg) ####
m_Carbon <- glmmTMB(
  log(CarbonStorage) ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species,
  data = es_data, family = gaussian()
)

run_diag(m_Carbon, "Carbon storage")
dci_effects(m_Carbon, "Carbon storage")
species_slopes(m_Carbon, "Carbon storage")
dci_imperv_slope(m_Carbon, "Carbon storage")
r2(m_Carbon)


#### Transpiration — Gaussian on log, +0.1 constant for the one zero ####
## Because the constant sits inside the log, emmeans returns contrasts on the
## log(mu + 0.1) scale. Do NOT convert them with (exp(est)-1)*100. Use water
## interception for a back-transformable effect size.
m_Transp <- glmmTMB(
  log(Transpiration + 0.1) ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species,
  data = es_data, family = gaussian()
)

run_diag(m_Transp, "Transpiration")
dci_effects(m_Transp, "Transpiration")
species_slopes(m_Transp, "Transpiration")
dci_imperv_slope(m_Transp, "Transpiration")
r2(m_Transp)


#### Water interception — Gaussian on log, no constant (min 7.2 l/yr) ####
m_WaterInt <- glmmTMB(
  log(WaterIntercepted) ~ DCI * Imperv25_s + Species + DCI:Species + Imperv25_s:Species +
    SoilMoisture_s * Species + SoilCompaction_s * Species,
  data = es_data, family = gaussian()
)

run_diag(m_WaterInt, "Water interception")
dci_effects(m_WaterInt, "Water interception")
species_slopes(m_WaterInt, "Water interception")
dci_imperv_slope(m_WaterInt, "Water interception")
r2(m_WaterInt)

## CAVEAT FOR THE WRITE-UP: i-Tree derives transpiration AND water interception


#### SANITY CHECK AGAINST THE THESIS TABLES ####
## Height's DCI term should be chisq = 8.19, p = 0.004 (Table 2-1). If it prints
## 5.88 / 0.015, a Distance_s term has crept back into the model.
.chk <- as.data.frame(car::Anova(m_Height, type = "II"))
message(sprintf("\nCHECK  Height DCI: chisq = %.2f, p = %.3f   (expect 8.19 / 0.004)",
                .chk["DCI", grep("Chisq|Chi", names(.chk))[1]],
                .chk["DCI", grep("^Pr",       names(.chk))[1]]))

################################################################################
##  ATLANTA CHAPTER — FIGURES  (Option A: unified species-interaction models)
################################################################################

library(tidyverse)
library(patchwork)
library(scales)
library(emmeans)
library(car)

################################################################################
##  SHARED THEME, PALETTES, LABELS
################################################################################
dci_colors <- c("Zone_One" = "#D4A843", "Zone_Five" = "#4B9AC7")   # amber / blue
dci_labels <- c("Zone_One" = "Zone 1\n(Prosperous)",
                "Zone_Five" = "Zone 5\n(Distressed)")
dci_leg    <- c("Zone_One" = "Zone 1 (Prosperous)",
                "Zone_Five" = "Zone 5 (Distressed)")

# genus -> abbreviated binomial, used for italic facet strips & axis text
spp_axis <- c(Acer = "A. rubrum", Lagerstroemia = "L. indica",
              Liriodendron = "L. tulipifera", Magnolia = "M. virginiana",
              Quercus = "Q. alba")

theme_atl <- function() {
  theme_classic(base_size = 11) +
    theme(
      axis.title       = element_text(size = 10),
      axis.text        = element_text(size = 9, color = "black"),
      strip.text       = element_text(size = 8.5, face = "italic"),
      strip.background = element_blank(),
      legend.position  = "none",
      plot.title       = element_text(size = 10, face = "bold", hjust = 0),
      panel.grid       = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      plot.margin      = margin(6, 6, 4, 4)
    )
}

################################################################################
##  STAT HELPERS  — all p-values come from the fitted Option A models
################################################################################

## p-value formatter and significance star
fmt_p  <- function(p) ifelse(is.na(p), "NA",
                      ifelse(p < 0.001, "p < 0.001",
                             paste0("p = ", formatC(p, format = "f", digits = 3))))
p_star <- function(p) ifelse(is.na(p), "",
                      ifelse(p < 0.001, "***",
                      ifelse(p < 0.01,  "**",
                      ifelse(p < 0.05,  "*", "ns"))))
lab_p  <- function(p) paste0(p_star(p), " (", fmt_p(p), ")")

## first matching column name from a set of candidates
.pick <- function(df, cands) intersect(cands, names(df))[1]

## Type-II omnibus test for a term. vars = c("DCI") for a main effect,
## c("DCI","Species") for an interaction (order-insensitive).
anova_stat <- function(model, vars) {
  tryCatch({
    av  <- as.data.frame(car::Anova(model, type = "II"))
    chi <- grep("Chisq|Chi", names(av), value = TRUE)[1]
    dfc <- grep("Df",        names(av), value = TRUE)[1]
    pc  <- grep("^Pr",       names(av), value = TRUE)[1]
    rn  <- strsplit(rownames(av), ":")
    hit <- which(vapply(rn, function(z) setequal(z, vars), logical(1)))
    if (!length(hit)) return(NULL)
    list(chisq = as.numeric(av[[chi]][hit]),
         df    = as.integer(av[[dfc]][hit]),
         p     = as.numeric(av[[pc]][hit]))
  }, error = function(e) {
    # Report why a label came out blank instead of failing silently
    message("anova_stat failed (term: ", paste(vars, collapse = ":"), "): ", conditionMessage(e))
    NULL
  })
}

## two-line omnibus annotation:  "chi2(df) = x.xx\n p = ... **"
omni_txt <- function(model, vars, tag = NULL) {
  s <- anova_stat(model, vars)
  if (is.null(s)) return("")
  head <- if (is.null(tag)) "" else paste0(tag, "\n")
  paste0(head, "\u03C7\u00B2(", s$df, ") = ",
         formatC(s$chisq, format = "f", digits = 2),
         "\n", fmt_p(s$p), " ", p_star(s$p))
}

## compact "DCI x species: chi2(df) = x.xx, p = ... *" label
dcixsp_lab <- function(model) {
  s <- anova_stat(model, c("DCI", "Species"))
  if (is.null(s)) return("")
  paste0("DCI\u00D7species: \u03C7\u00B2(", s$df, ") = ",
         formatC(s$chisq, format = "f", digits = 2),
         ", ", fmt_p(s$p), " ", p_star(s$p))
}

## per-species Zone5-vs-Zone1 contrast p-values (emmeans within species)
within_spp_p <- function(model, ordinal = FALSE) {
  tryCatch({
    emm <- if (ordinal) emmeans(model, ~ DCI | Species, mode = "latent")
           else          emmeans(model, ~ DCI | Species, type = "response")
    d <- as.data.frame(summary(pairs(emm, reverse = TRUE), infer = TRUE))
    names(d)[names(d) == "p.value"] <- "p"
    data.frame(Species = as.character(d$Species), p = as.numeric(d$p),
               label = lab_p(as.numeric(d$p)), stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}

## model-estimated marginal zone means (response scale), for the diamond overlay.
## kind = "log" (clean log link/response), "logconst" (response was log(y+const)).
zone_emm <- function(model, kind = "log", const = 0) {
  d <- tryCatch(as.data.frame(emmeans(model, ~ DCI, type = "response")),
                error = function(e)
                  tryCatch(as.data.frame(emmeans(model, ~ DCI)),
                           error = function(e2) NULL))
  if (is.null(d)) return(NULL)
  est <- .pick(d, c("response", "rate", "emmean", "prob"))
  lo  <- .pick(d, c("asymp.LCL", "lower.CL"))
  hi  <- .pick(d, c("asymp.UCL", "upper.CL"))
  if (is.na(est) || is.na(lo) || is.na(hi)) return(NULL)
  out <- data.frame(DCI = d$DCI, fit = d[[est]], lo = d[[lo]], hi = d[[hi]])
  if (kind == "logconst") { out$fit <- exp(out$fit) - const
                            out$lo  <- exp(out$lo)  - const
                            out$hi  <- exp(out$hi)  - const }
  out
}

## per-species marginal zone means (response scale) for faceted overlays
zone_emm_byspp <- function(model, kind = "log", const = 0) {
  d <- tryCatch(as.data.frame(emmeans(model, ~ DCI | Species, type = "response")),
                error = function(e) NULL)
  if (is.null(d)) return(NULL)
  est <- .pick(d, c("response", "rate", "emmean", "prob"))
  lo  <- .pick(d, c("asymp.LCL", "lower.CL"))
  hi  <- .pick(d, c("asymp.UCL", "upper.CL"))
  if (is.na(est) || is.na(lo) || is.na(hi)) return(NULL)
  out <- data.frame(DCI = d$DCI, Species = d$Species,
                    fit = d[[est]], lo = d[[lo]], hi = d[[hi]])
  if (kind == "logconst") { out$fit <- exp(out$fit) - const
                            out$lo  <- exp(out$lo)  - const
                            out$hi  <- exp(out$hi)  - const }
  out
}

## species-specific slope table (emtrends) for the forest plots
slopes_df <- function(model, var) {
  d <- as.data.frame(summary(emtrends(model, ~ Species, var = var), infer = TRUE))
  names(d)[grep("\\.trend$", names(d))]             <- "slope"
  names(d)[names(d) %in% c("asymp.LCL", "lower.CL")] <- "LCL"
  names(d)[names(d) %in% c("asymp.UCL", "upper.CL")] <- "UCL"
  names(d)[names(d) == "p.value"]                    <- "p"
  d$Species <- factor(d$Species)
  d
}

################################################################################
##  REUSABLE PANEL BUILDERS
################################################################################

## Continuous response, boxplot + jitter by DCI zone, with (optional) the
## model-estimated marginal zone mean overlaid as a black diamond +/- 95% CI.
box_zone <- function(data, yvar, ylab, title,
                     model = NULL, kind = "log", const = 0,
                     ytrans = "log10", ybreaks = waiver(), yacc = 1,
                     annot = "") {
  zemm <-  NULL
  if (!is.null(zemm)) {                       # keep positive for the log axis
    zemm$fit <- pmax(zemm$fit, 1e-6)
    zemm$lo  <- pmax(zemm$lo,  1e-6)
    zemm$hi  <- pmax(zemm$hi,  1e-6)
  }
  p <- ggplot(data, aes(x = DCI, y = .data[[yvar]], fill = DCI)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.70, width = 0.5,
                 color = "black", linewidth = 0.35) +
    geom_jitter(aes(color = DCI), width = 0.15, size = 1.0,
                alpha = 0.40, shape = 16)
  if (!is.null(zemm)) {
    p <- p +
      geom_errorbar(data = zemm, aes(x = DCI, ymin = lo, ymax = hi),
                    inherit.aes = FALSE, width = 0.12, linewidth = 0.55) +
      geom_point(data = zemm, aes(x = DCI, y = fit), inherit.aes = FALSE,
                 shape = 23, size = 2.4, fill = "black", color = "black")
  }
  p +
    scale_fill_manual(values = dci_colors, labels = dci_leg, name = "DCI zone") +
    scale_color_manual(values = dci_colors, labels = dci_leg, name = "DCI zone") +
    scale_x_discrete(labels = dci_labels) +
    scale_y_continuous(trans = ytrans, breaks = ybreaks,
                       labels = label_number(accuracy = yacc)) +
    annotate("text", x = 1.5, y = Inf, vjust = 1.4, size = 2.9, label = annot) +
    labs(x = NULL, y = ylab, title = title) +
    theme_atl()
}

## Faceted (species) boxplot + jitter, per-facet contrast p annotation, and
## per-species marginal zone means overlaid.
box_zone_facet <- function(data, yvar, ylab, title, model,
                           kind = "log", const = 0,
                           ytrans = "log10", ybreaks = waiver(), yacc = 1,
                           strip = TRUE) {
  ann  <- within_spp_p(model, ordinal = FALSE)
  zemm <- NULL
  if (!is.null(zemm)) {                       # keep positive for the log axis
    zemm$fit <- pmax(zemm$fit, 1e-6)
    zemm$lo  <- pmax(zemm$lo,  1e-6)
    zemm$hi  <- pmax(zemm$hi,  1e-6)
  }
  p <- ggplot(data, aes(x = DCI, y = .data[[yvar]], fill = DCI)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.70, width = 0.55,
                 color = "black", linewidth = 0.35) +
    geom_jitter(aes(color = DCI), width = 0.15, size = 0.9,
                alpha = 0.40, shape = 16)
  if (!is.null(zemm)) {
    p <- p +
      geom_errorbar(data = zemm, aes(x = DCI, ymin = lo, ymax = hi),
                    inherit.aes = FALSE, width = 0.12, linewidth = 0.5) +
      geom_point(data = zemm, aes(x = DCI, y = fit), inherit.aes = FALSE,
                 shape = 23, size = 1.9, fill = "black", color = "black")
  }
  if (!is.null(ann)) {
    p <- p + geom_text(data = ann, aes(x = 1.5, y = Inf, label = label),
                       inherit.aes = FALSE, vjust = 1.5, size = 2.5)
  }
  p +
    facet_wrap(~ Species, nrow = 1, labeller = as_labeller(spp_axis)) +
    scale_fill_manual(values = dci_colors, labels = dci_leg, name = "DCI zone") +
    scale_color_manual(values = dci_colors, labels = dci_leg, name = "DCI zone") +
    scale_x_discrete(labels = c("Zone_One" = "Z1", "Zone_Five" = "Z5")) +
    scale_y_continuous(trans = ytrans, breaks = ybreaks,
                       labels = label_number(accuracy = yacc)) +
    labs(x = NULL, y = ylab, title = title) +
    theme_atl() +
    theme(strip.text = if (strip) element_text(face = "italic", size = 8)
                       else element_blank())
}

## Forest plot of species-specific slopes (emtrends output).
forest_slope <- function(df, title, xlab, interaction_txt) {
  df$sig  <- df$p < 0.05
  df$star <- p_star(df$p)
  ggplot(df, aes(x = slope, y = Species)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = LCL, xmax = UCL, color = sig),
                   height = 0.20, linewidth = 0.6) +
    geom_point(aes(color = sig, fill = sig), shape = 21, size = 2.6, stroke = 0.6) +
    geom_text(aes(x = UCL, label = star), hjust = -0.35, vjust = 0.75, size = 3) +
    scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "gray55")) +
    scale_fill_manual(values  = c(`TRUE` = "#B2182B", `FALSE` = "white")) +
    scale_y_discrete(labels = spp_axis) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
             size = 2.7, label = interaction_txt) +
    labs(x = xlab, y = NULL, title = title) +
    guides(color = "none", fill = "none") +   # never show the TRUE/FALSE sig key
    theme_atl() +
    theme(legend.position = "none",
          axis.text.y = element_text(face = "italic", size = 8))
}

################################################################################
##  FIGURE 2-2  —  COMMUNITY TREE SIZE & CONDITION ACROSS DCI ZONES
##  Pooled zone effect (headline where DCI x Species is NOT significant).
##  Panels: A Height  B DBH  C Canopy volume  D Condition
################################################################################
try({

  # A Height: repeated measures. DCI x Species IS significant, so the pooled
  #   view is misleading -> no diamond; annotate with the interaction & point
  #   the reader to the species panel (Fig 2-3).
  p22a <- box_zone(rd, "Height", "Height (m)", "A  Height",
                   model = NULL, ytrans = "log10",
                   ybreaks = c(1, 2, 5, 10, 20), yacc = 1,
                   annot = paste0(omni_txt(m_Height, c("DCI", "Species"),
                                           tag = "Zone \u00D7 species")
                                  ))

  # B DBH (2025 only): DCI x Species ns -> pooled headline (here ns)
  p22b <- box_zone(rd_2025, "DBH", "DBH (cm)", "B  DBH",
                   model = m_DBH, ytrans = "log10",
                   ybreaks = c(1, 2, 5, 10, 20), yacc = 1,
                   annot = omni_txt(m_DBH, c("DCI"), tag = "Zone effect"))

  # C Canopy volume (2025 only): pooled headline, ~37% lower in Zone 5
  p22c <- box_zone(rd_2025, "CanopyVolume",
                   expression("Canopy volume (m"^3*")"), "C  Canopy volume",
                   model = m_CV, ytrans = "log1p",
                   ybreaks = c(0, 5, 20, 100, 500), yacc = 1,
                   annot = omni_txt(m_CV, c("DCI"), tag = "Zone effect"))

  # D Condition (ordinal): proportional stacked bar + omnibus DCI test
  cond_sum <- rd %>% filter(!is.na(Condition)) %>%
    count(DCI, Condition) %>% group_by(DCI) %>%
    mutate(prop = n / sum(n)) %>% ungroup() %>%
    mutate(Condition = factor(Condition, levels = 5:1),
           DCI = factor(DCI, levels = c("Zone_One", "Zone_Five")))
  cond_cols <- c("5" = "#1D6B3C", "4" = "#6AB187", "3" = "#D4C44A",
                 "2" = "#E08B35", "1" = "#C0392B")
  p22d <- ggplot(cond_sum, aes(x = DCI, y = prop, fill = Condition)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = cond_cols, name = "Condition\nscore",
                      breaks = c("5", "4", "3", "2", "1"),
                      labels = c("5 (Excellent)", "4", "3", "2", "1 (Poor)")) +
    scale_x_discrete(labels = dci_labels) +
    # Headroom above 100% so the zone-effect label sits clear of the bars
    scale_y_continuous(labels = label_percent(), breaks = seq(0, 1, 0.25),
                       limits = c(0, 1.2), expand = expansion(mult = c(0, 0))) +
    annotate("text", x = 1.5, y = 1.03, vjust = 0, size = 2.9,
             label = omni_txt(m_Cond, c("DCI"), tag = "Zone effect")) +
    labs(x = NULL, y = "Proportion of trees", title = "D  Condition") +
    theme_atl() +
    theme(legend.position = "right", legend.key.size = unit(0.32, "cm"),
          legend.text = element_text(size = 7.5),
          legend.title = element_text(size = 7.5))

  fig22 <- (p22a | p22b | p22c | p22d) +
    plot_layout(guides = "collect")

  ggsave("Figure2-2_TreeSize_Community.pdf", fig22,
         device = cairo_pdf,
         width = 10.5, height = 7.6, units = "in")
  # TIFF copy for embedding in the Word manuscript (300 dpi, LZW-compressed)
  ggsave("Figure2-2_TreeSize_Community.tiff", fig22, device = "tiff", type = "cairo",
         width = 10.5, height = 7.6, units = "in",
         dpi = 300, compression = "lzw", bg = "white")
  message("Figure 2-2 saved.")
}, silent = FALSE)

################################################################################
##  FIGURE 2-3  —  SPECIES-SPECIFIC TREE SIZE WITHIN EACH DCI ZONE
##  Per-species Zone5-vs-Zone1 contrast (emmeans ~ DCI | Species).
##  Rows: A Height  B DBH  C Canopy volume
################################################################################
try({

  p23a <- box_zone_facet(rd, "Height", "Height (m)",
                         paste0("A  Height   |   ", dcixsp_lab(m_Height)),
                         model = m_Height, ytrans = "log10",
                         ybreaks = c(1, 2, 5, 10, 20), yacc = 1, strip = TRUE)

  p23b <- box_zone_facet(rd_2025, "DBH", "DBH (cm)",
                         paste0("B  DBH   |   ", dcixsp_lab(m_DBH)),
                         model = m_DBH, ytrans = "log10",
                         ybreaks = c(1, 2, 5, 10, 20), yacc = 1, strip = FALSE)

  p23c <- box_zone_facet(rd_2025, "CanopyVolume",
                         expression("Canopy volume (m"^3*")"),
                         paste0("C  Canopy volume   |   ", dcixsp_lab(m_CV)),
                         model = m_CV, ytrans = "log1p",
                         ybreaks = c(0, 5, 20, 100, 500), yacc = 1, strip = FALSE)

  fig23 <- (p23a / p23b / p23c) +
    plot_layout(guides = "collect")

  ggsave("Figure2-3_TreeSize_Species.pdf", fig23,
         device = cairo_pdf, width = 9.5, height = 7.2, units = "in")
  # TIFF copy for embedding in the Word manuscript (300 dpi, LZW-compressed)
  ggsave("Figure2-3_TreeSize_Species.tiff", fig23, device = "tiff", type = "cairo",
         width = 9.5, height = 7.2, units = "in",
         dpi = 300, compression = "lzw", bg = "white")
  message("Figure 2-3 saved.")
}, silent = FALSE)

################################################################################
##  FIGURE 2-4  —  COMMUNITY INSECT-PEST METRICS ACROSS DCI ZONES
##  A Aphid density   B Scale density   C Foliar herbivory (ordinal)
##  Aphid & scale use rd_pest (Magnolia excluded) to match the fitted models.
################################################################################
try({

  # A Aphids: strong pooled zone effect (~4.7x higher in Zone 5)
  aph_max <- max(rd_pest$AphidCount, na.rm = TRUE)
  p24a <- box_zone(rd_pest, "AphidCount", "Aphids per sample", "A  Aphid density",
                   model = m_Aphid, kind = "log",
                   ytrans = "log1p", ybreaks = c(0, 1, 5, 20, 100), yacc = 1,
                   annot = paste0(omni_txt(m_Aphid, c("DCI"), tag = "Zone effect"),
                                  "\n(Magnolia excluded)"))

  # B Scale: NO pooled zone effect; the signal is DCI x impervious (-> Fig 2-5)
  p24b <- box_zone(rd_pest, "ScaleCount", "Scale insects per sample",
                   "B  Scale density",
                   model = m_Scale, kind = "log",
                   ytrans = "log1p", ybreaks = c(0, 1, 5, 20, 100), yacc = 1,
                   annot = paste0(omni_txt(m_Scale, c("DCI"), tag = "Zone effect")))

  # C Herbivory (ordinal 0-5): proportional stacked bar + omnibus DCI test
  herb_lv  <- sort(unique(as.integer(as.character(rd$HerbRate[!is.na(rd$HerbRate)]))))
  herb_sum <- rd %>% filter(!is.na(HerbRate)) %>%
    count(DCI, HerbRate) %>% group_by(DCI) %>%
    mutate(prop = n / sum(n)) %>% ungroup() %>%
    mutate(HerbRate = factor(HerbRate, levels = rev(herb_lv)),
           DCI = factor(DCI, levels = c("Zone_One", "Zone_Five")))
  herb_cols <- setNames(
    colorRampPalette(c("#E9EDC9", "#DDA15E", "#9C2A1E"))(length(herb_lv)),
    as.character(herb_lv))
  p24c <- ggplot(herb_sum, aes(x = DCI, y = prop, fill = HerbRate)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = herb_cols, name = "Herbivory\nscore",
                      breaks = as.character(rev(herb_lv))) +
    scale_x_discrete(labels = dci_labels) +
    # Headroom above 100% so the zone-effect label sits clear of the bars
    scale_y_continuous(labels = label_percent(), breaks = seq(0, 1, 0.25),
                       limits = c(0, 1.2), expand = expansion(mult = c(0, 0))) +
    annotate("text", x = 1.5, y = 1.03, vjust = 0, size = 2.9,
             label = omni_txt(m_Herb, c("DCI"), tag = "Zone effect")) +
    labs(x = NULL, y = "Proportion of trees", title = "C  Foliar herbivory") +
    theme_atl() +
    theme(legend.position = "right", legend.key.size = unit(0.32, "cm"),
          legend.text = element_text(size = 7.5),
          legend.title = element_text(size = 7.5))

  fig24 <- (p24a | p24b | p24c) +
    plot_layout(guides = "collect")

  ggsave("Figure2-4_PestMetrics_Community.pdf", fig24,
         device = cairo_pdf, width = 9.5, height = 5.7, units = "in")
  # TIFF copy for embedding in the Word manuscript (300 dpi, LZW-compressed)
  ggsave("Figure2-4_PestMetrics_Community.tiff", fig24, device = "tiff", type = "cairo",
         width = 9.5, height = 5.7, units = "in",
         dpi = 300, compression = "lzw", bg = "white")
  
  message("Figure 2-4 saved.")
}, silent = FALSE)

################################################################################
##  FIGURE 2-5  —  LOCAL SITE-CONDITION EFFECTS THAT VARY BY ZONE OR SPECIES
##  (the Option A payoff)
##  A  DCI x impervious for scale insects  (the one significant zone interaction)
##  B  Impervious slope by species — DBH        (Imperv x Species significant)
##  C  Impervious slope by species — Condition  (Imperv x Species significant)
##  D  Soil-moisture slope by species — Aphids  (Species x SoilMoisture sig.)
################################################################################
try({

  ## ---- Panel A: predicted scale density vs impervious, by zone ----
  ## Map the standardized predictor back to raw % using the std2 constants
  ## (std2 used mean & 2*SD of rd$Impervious25).
  mu_imp <- mean(rd$Impervious25, na.rm = TRUE)
  sd_imp <- sd(rd$Impervious25,   na.rm = TRUE)
  seq_s  <- seq(min(rd_pest$Imperv25_s, na.rm = TRUE),
                max(rd_pest$Imperv25_s, na.rm = TRUE), length.out = 60)

  emm_sc <- as.data.frame(emmeans(m_Scale, ~ DCI, by = "Imperv25_s",
                                  at = list(Imperv25_s = seq_s),
                                  type = "response"))
  names(emm_sc)[names(emm_sc) %in% c("response", "rate", "emmean")] <- "fit"
  names(emm_sc)[names(emm_sc) %in% c("asymp.LCL", "lower.CL")]      <- "lo"
  names(emm_sc)[names(emm_sc) %in% c("asymp.UCL", "upper.CL")]      <- "hi"
  emm_sc$Impervious25 <- emm_sc$Imperv25_s * 2 * sd_imp + mu_imp

  # per-zone slope tests + the DCI x impervious interaction test (for the label)
  sc_zone <- as.data.frame(summary(
    emtrends(m_Scale, ~ DCI, var = "Imperv25_s"), infer = TRUE))
  names(sc_zone)[names(sc_zone) == "p.value"] <- "p"
  pz1 <- sc_zone$p[sc_zone$DCI == "Zone_One"]
  pz5 <- sc_zone$p[sc_zone$DCI == "Zone_Five"]
  sci <- anova_stat(m_Scale, c("DCI", "Imperv25_s"))
  a_lab <- paste0("DCI\u00D7impervious: \u03C7\u00B2(", sci$df, ") = ",
                  formatC(sci$chisq, format = "f", digits = 2), ", ",
                  fmt_p(sci$p), " ", p_star(sci$p),
                  "\nZone 1 slope ", fmt_p(pz1), " ", p_star(pz1),
                  "  |  Zone 5 slope ", fmt_p(pz5), " ", p_star(pz5))

  p25a <- ggplot() +
    geom_ribbon(data = emm_sc, aes(x = Impervious25, ymin = lo, ymax = hi,
                                   fill = DCI), alpha = 0.15) +
    geom_line(data = emm_sc, aes(x = Impervious25, y = fit, color = DCI),
              linewidth = 1.0) +
    geom_point(data = rd_pest, aes(x = Impervious25, y = ScaleCount, color = DCI),
               size = 1.2, alpha = 0.40, shape = 16) +
    scale_color_manual(values = dci_colors, labels = dci_leg, name = "DCI zone") +
    scale_fill_manual(values  = dci_colors, labels = dci_leg, name = "DCI zone") +
    scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 20, 100)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 1.3, size = 2.6,
             label = a_lab) +
    labs(x = "Local impervious cover (%)", y = "Scale insects per sample",
         title = "A  Scale insects: DCI \u00D7 impervious") +
    theme_atl() +
    theme(legend.position      = c(0.98, 0.78),
          legend.justification = c(1, 1),
          legend.title         = element_text(size = 8),
          legend.text          = element_text(size = 7.5),
          legend.key.size      = unit(0.4, "cm"),
          legend.background    = element_rect(fill = alpha("white", 0.75),
                                              color = NA),
          legend.margin        = margin(2, 4, 2, 4))

  ## ---- Panels B-D: species-specific slope forest plots (emtrends) ----
  int_lab <- function(model, vars, name) {
    s <- anova_stat(model, vars)
    paste0(name, ": \u03C7\u00B2(", s$df, ") = ",
           formatC(s$chisq, format = "f", digits = 2), ", ",
           fmt_p(s$p), " ", p_star(s$p))
  }

  p25b <- forest_slope(
    slopes_df(m_DBH, "Imperv25_s"), "B  DBH \u2013 impervious slope",
    "Slope (log DBH per 2 SD impervious)",
    int_lab(m_DBH, c("Imperv25_s", "Species"), "Impervious\u00D7species"))

  p25c <- forest_slope(
    slopes_df(m_Cond, "Imperv25_s"), "C  Condition \u2013 impervious slope",
    "Slope (latent logit per 2 SD impervious)",
    int_lab(m_Cond, c("Imperv25_s", "Species"), "Impervious\u00D7species"))

  p25d <- forest_slope(
    slopes_df(m_Aphid, "SoilMoisture_s"), "D  Aphids \u2013 soil-moisture slope",
    "Slope (log aphids per 2 SD moisture)",
    int_lab(m_Aphid, c("SoilMoisture_s", "Species"), "Moisture\u00D7species"))

  fig25 <- (p25a | p25b) / (p25c | p25d)

  ggsave("Figure2-5_SiteInteractions.pdf", fig25,
         device = cairo_pdf, width = 9.0, height = 7.0, units = "in")
  # TIFF copy for embedding in the Word manuscript (300 dpi, LZW-compressed)
  ggsave("Figure2-5_SiteInteractions.tiff", fig25, device = "tiff", type = "cairo",
         width = 9.0, height = 7.0, units = "in",
         dpi = 300, compression = "lzw", bg = "white")
  message("Figure 2-5 saved.")
}, silent = FALSE)




################################################################################
##  FIGURE 2-6  —  ECOSYSTEM SERVICES ACROSS DCI ZONES  (REDESIGNED)




box_zone_es <- function(data, yvar, ylab, title, model,
                        kind = "log", const = 0,
                        ytrans = "log10", ybreaks = waiver(), yacc = 1,
                        annot = "") {
  zemm <- zone_emm(model, kind = kind, const = const)
  if (!is.null(zemm)) {                        # keep positive for the log axis
    zemm$fit <- pmax(zemm$fit, 1e-6)
    zemm$lo  <- pmax(zemm$lo,  1e-6)
    zemm$hi  <- pmax(zemm$hi,  1e-6)
  }
  p <- ggplot(data, aes(x = DCI, y = .data[[yvar]], fill = DCI)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.70, width = 0.5,
                 color = "black", linewidth = 0.35) +
    geom_jitter(aes(color = DCI), width = 0.15, size = 1.0,
                alpha = 0.40, shape = 16)
  if (!is.null(zemm)) {
    p <- p +
      geom_errorbar(data = zemm, aes(x = DCI, ymin = lo, ymax = hi),
                    inherit.aes = FALSE, width = 0.12, linewidth = 0.55) +
      geom_point(data = zemm, aes(x = DCI, y = fit), inherit.aes = FALSE,
                 shape = 23, size = 2.4, fill = "black", color = "black")
  }
  p +
    scale_fill_manual(values = dci_colors, labels = dci_leg, name = "DCI zone") +
    scale_color_manual(values = dci_colors, labels = dci_leg, name = "DCI zone") +
    scale_x_discrete(labels = dci_labels) +
    scale_y_continuous(trans = ytrans, breaks = ybreaks,
                       labels = label_number(accuracy = yacc)) +
    annotate("text", x = 1.5, y = Inf, vjust = 1.3, size = 2.7, label = annot) +
    labs(x = NULL, y = ylab, title = title) +
    theme_atl()
}

## pooled DCI annotation, plus the (ns) DCI x species status that justifies
## pooling in the first place — keeps the reader from asking "why not by species?"
pooled_lab <- function(model) {
  s  <- anova_stat(model, c("DCI"))
  si <- anova_stat(model, c("DCI", "Species"))
  if (is.null(s)) return("")
  paste0("\u03C7\u00B2(", s$df, ") = ", formatC(s$chisq, format = "f", digits = 2),
         "\n", fmt_p(s$p), " ", p_star(s$p),
         if (!is.null(si)) paste0("\nDCI\u00D7species ", p_star(si$p)) else "")
}

## impervious x species annotation for the forest panels
impxsp_lab <- function(model) {
  s <- anova_stat(model, c("Imperv25_s", "Species"))
  if (is.null(s)) return("")
  paste0("Impervious\u00D7species: \u03C7\u00B2(", s$df, ") = ",
         formatC(s$chisq, format = "f", digits = 2),
         ", ", fmt_p(s$p), " ", p_star(s$p))
}


try({
  
  ## Display-only offsets so zeros survive the log10 axis (PM2.5 has 6 zeros,
  ## transpiration has 1). The diamond overlay is separate: it is back-
  ## transformed on the true response scale using each model's own constant.
  es_fig <- es_data %>%
    mutate(PM25_plot   = PM25_removal  + 0.001,
           Transp_plot = Transpiration + 0.01)
  
  ## ---- Row 1: pooled community zone effect, all four services --------------
  p26a <- box_zone_es(es_fig, "PM25_plot",
                      expression("Annual PM"[2.5]*" removal (g yr"^{-1}*")"),
                      "A  PM2.5 removal",
                      model = m_PM25, kind = "logconst", const = 0.05,
                      ytrans = "log10", yacc = 0.01,
                      annot = pooled_lab(m_PM25))
  
  p26b <- box_zone_es(es_fig, "CarbonStorage",
                      "Carbon storage (kg)",
                      "B  Carbon storage",
                      model = m_Carbon, kind = "log",
                      ytrans = "log10", yacc = 0.1,
                      annot = pooled_lab(m_Carbon))
  
  p26c <- box_zone_es(es_fig, "Transp_plot",
                      expression("Transpiration (l yr"^{-1}*")"),
                      "C  Transpiration",
                      model = m_Transp, kind = "logconst", const = 0.1,
                      ytrans = "log10", yacc = 0.01,
                      annot = pooled_lab(m_Transp))
  
  p26d <- box_zone_es(es_fig, "WaterIntercepted",
                      expression("Water interception (l yr"^{-1}*")"),
                      "D  Water interception",
                      model = m_WaterInt, kind = "log",
                      ytrans = "log10", yacc = 1,
                      annot = pooled_lab(m_WaterInt))
  
  ## ---- Row 2: ONLY the services with a significant impervious x species ----
  p26e <- forest_slope(
    slopes_df(m_PM25, "Imperv25_s"),
    "E  PM2.5 \u2013 impervious slope by species",
    "Slope (log PM2.5 per 2 SD impervious)",
    impxsp_lab(m_PM25))
  
  p26f <- forest_slope(
    slopes_df(m_Carbon, "Imperv25_s"),
    "F  Carbon \u2013 impervious slope by species",
    "Slope (log carbon per 2 SD impervious)",
    impxsp_lab(m_Carbon))
  
  fig26 <- ((p26a | p26b | p26c | p26d) / (p26e | p26f)) +
    plot_layout(heights = c(1, 0.85), guides = "collect")
  
  ggsave("Figure2-6_EcosystemServices.pdf", fig26,
         device = cairo_pdf, width = 11, height = 7.2, units = "in")
  # TIFF copy for embedding in the Word manuscript (300 dpi, LZW-compressed)
  ggsave("Figure2-6_EcosystemServices.tiff", fig26, device = "tiff", type = "cairo",
         width = 11, height = 7.2, units = "in",
         dpi = 300, compression = "lzw", bg = "white")
  message("Figure 2-6 saved.")
  
}, silent = FALSE)
