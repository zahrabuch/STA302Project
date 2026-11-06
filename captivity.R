remotes::install_github("datawookie/hagr")
library(hagr)
library(dplyr)
library(ggplot2)

max_missing <- 2500
min_n       <- 10      # smallest class that gets a panel
use_log     <- FALSE   # set TRUE only if transformed axes are allowed

dat <- age
cat("Rows loaded:", nrow(dat), "\n")

#keep columns with at most 2500 missing values
resp <- "maximum_longevity_yr"
keep <- colSums(is.na(dat)) <= max_missing
cat("Dropped for missingness:", paste(names(dat)[!keep], collapse = ", "), "\n")
dat <- dat[, keep]
stopifnot(resp %in% names(dat), "adult_mass_g" %in% names(dat), "class" %in% names(dat))

# remove the listed columns 
drop_cols <- c("hagrid", "common_name", "data_quality", "source",
               "sample_size", "domain", "kingdom", "phylum")
dat <- dat %>% select(-any_of(drop_cols))

# summary of remaining predictors (BEFORE dropping missing rows)
skew <- function(x) { x <- x[!is.na(x)]; mean((x - mean(x))^3) / sd(x)^3 }

preds   <- setdiff(names(dat), resp)
num_var <- preds[sapply(dat[preds], is.numeric)]
cat_var <- setdiff(preds, num_var)

num_summary <- data.frame(
  Variable      = num_var,
  N_missing     = sapply(dat[num_var], function(x) sum(is.na(x))),
  Pct_missing   = sapply(dat[num_var], function(x) round(100 * mean(is.na(x)), 1)),
  Mean          = sapply(dat[num_var], mean,   na.rm = TRUE),
  Median        = sapply(dat[num_var], median, na.rm = TRUE),
  Variance      = sapply(dat[num_var], var,    na.rm = TRUE),
  Skew          = sapply(dat[num_var], skew),
  row.names     = NULL
)
num_summary[, c("Mean", "Median", "Variance", "Skew")] <-
  lapply(num_summary[, c("Mean", "Median", "Variance", "Skew")], signif, 5)
print(num_summary)

# levels and counts for each categorical variable (before dropping missing rows) 
cat_var <- setdiff(names(dat)[!sapply(dat, is.numeric)], resp)

level_counts <- do.call(rbind, lapply(cat_var, function(v) {
  x  <- na_if(trimws(as.character(dat[[v]])), "")      # empty text counts as missing
  tb <- table(x, useNA = "ifany")
  data.frame(Variable = v,
             Level    = ifelse(is.na(names(tb)), "(missing)", names(tb)),
             Count    = as.integer(tb),
             Percent  = round(100 * as.integer(tb) / length(x), 1),
             row.names = NULL)
}))
level_counts <- level_counts[order(level_counts$Variable, -level_counts$Count), ]
print(level_counts, row.names = FALSE)
write.csv(level_counts, "categorical_level_counts.csv", row.names = FALSE)

# histogram of the response BEFORE dropping missing rows 
cat("Species with a recorded longevity:", sum(!is.na(dat[[resp]])), "\n")

p_hist_before <- ggplot(dat, aes(maximum_longevity_yr)) +
  geom_histogram(bins = 40, na.rm = TRUE, fill = "steelblue", colour = "white") +
  labs(x = "Maximum longevity (years)", y = "Number of species",
       title = "Maximum longevity (before dropping missing rows)")
print(p_hist_before)
ggsave("hist_longevity_before_dropping.png", p_hist_before,
       width = 7, height = 4.5, dpi = 300)

write.csv(num_summary, "predictor_summary_numeric.csv", row.names = FALSE)

# drop every row with a missing value
# treat empty text as missing too
dat <- dat %>% mutate(across(where(is.character), ~ na_if(trimws(.x), "")))
dat <- na.omit(dat)
cat("Complete rows:", nrow(dat), "\n")
cat("Any NA left:", anyNA(dat), "\n")      # should print FALSE

#convert specimen_origin to a categorical (factor)
dat$specimen_origin <- factor(dat$specimen_origin)

levels(dat$specimen_origin)        # see the categories
table(dat$specimen_origin)         # rows per category
class(dat$specimen_origin)         # should print "factor"

dat$class <- factor(dat$class)
dat <- droplevels(dat)
cat("Remaining columns:", paste(names(dat), collapse = ", "), "\n")

# histogram of the response 
p_hist <- ggplot(dat, aes(maximum_longevity_yr)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  labs(x = "Maximum longevity (years)", y = "Frequency",
       title = "Histogram Showing Distribution of Maximum Longevity in Years")
print(p_hist)
ggsave("hist_longevity.png", p_hist, width = 7, height = 4.5, dpi = 300)

# adult mass vs longevity, one panel per class 
counts <- table(dat$class)
print(sort(counts, decreasing = TRUE))
keep_cls <- names(counts)[counts >= min_n]

plot_dat <- dat %>%
  filter(class %in% keep_cls) %>%
  droplevels() %>%
  mutate(class_lab = factor(paste0(class, " (n = ", counts[as.character(class)], ")")))

p_cls <- ggplot(plot_dat, aes(adult_mass_g, maximum_longevity_yr)) +
  geom_point(alpha = 0.4, colour = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, colour = "firebrick", linewidth = 0.8) +
  facet_wrap(~ class_lab, scales = "free") +
  labs(x = "Adult mass (g)", y = "Maximum longevity (years)",
       title = "Adult mass vs maximum longevity by class") +
  theme_bw()

if (use_log) p_cls <- p_cls + scale_x_log10() + scale_y_log10()

print(p_cls)
n_panels <- length(unique(plot_dat$class))
ggsave("mass_vs_longevity_by_class.png", p_cls,
       width = 10, height = max(4, 2.6 * ceiling(n_panels / 3)), dpi = 300)

# ---- prepare the model data ----
library(forcats)
mod_dat <- dat

# discrete (count) variables, computed on all 4,219 species in `age`
count_of <- function(col) {
  tb <- table(age[[col]])
  as.integer(tb[as.character(mod_dat[[col]])])
}
mod_dat$species_count <- count_of("species")   # species sharing this epithet
mod_dat$genus_size    <- count_of("genus")     # species in the same genus
mod_dat$family_size   <- count_of("family")    # species in the same family

# remove the original text columns they replace
mod_dat <- mod_dat %>% select(-any_of(c("species", "genus", "family")))

# class and order: top 5 levels + "Other"
mod_dat$class <- fct_lump_n(factor(mod_dat$class), n = 5, other_level = "Other")
mod_dat$order <- fct_lump_n(factor(mod_dat$order), n = 5, other_level = "Other")

# any remaining text columns (e.g. specimen_origin) become factors
chr <- names(mod_dat)[sapply(mod_dat, is.character)]
mod_dat[chr] <- lapply(mod_dat[chr], factor)

for (v in c("class", "order")) { cat("\n", v, "\n"); print(table(mod_dat[[v]])) }
summary(mod_dat[c("species_count", "genus_size", "family_size")])

# use ALL leftover columns as predictors
preds <- setdiff(names(mod_dat), resp)
cat("\nPredictors in the model:", paste(preds, collapse = ", "), "\n")
cat("Number of predictors:", length(preds), "\n")

# check the model is estimable BEFORE fitting 
form <- reformulate(c(preds, "adult_mass_g:class"), response = resp)
X <- model.matrix(form, data = mod_dat)
cat("Observations:", nrow(mod_dat), " Model coefficients:", ncol(X),
    " Rank:", qr(X)$rank, "\n")

#fit the preliminary model: all predictors + adult_mass_g x class interaction
form <- reformulate(c(preds, "adult_mass_g:class"), response = resp)
fit  <- lm(form, data = mod_dat)

summary(fit)
anova(fit)
cat("Observations used:", nobs(fit), "\n")

# multicollinearity check (if the car package is installed)
if (requireNamespace("car", quietly = TRUE)) print(car::vif(fit))

#residual diagnostics
res_dat <- data.frame(
  fitted = fitted(fit),
  resid  = resid(fit),
  std    = rstandard(fit),
  class  = mod_dat$class[as.integer(names(resid(fit)))]
)

png("residual_diagnostics.png", width = 2000, height = 2000, res = 250)
par(mfrow = c(2, 2))
plot(fit)            # residuals vs fitted, Q-Q, scale-location, residuals vs leverage
dev.off()
par(mfrow = c(2, 2)); plot(fit); par(mfrow = c(1, 1))

# residuals vs fitted, coloured by class
p_rf <- ggplot(res_dat, aes(fitted, resid, colour = class)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values", y = "Residuals", colour = "Class",
       title = "Residuals vs fitted")
print(p_rf)
ggsave("resid_vs_fitted.png", p_rf, width = 7, height = 5, dpi = 300)

# histogram of residuals
p_rh <- ggplot(res_dat, aes(resid)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  labs(x = "Residuals", y = "Count", title = "Distribution of residuals")
print(p_rh)
ggsave("resid_hist.png", p_rh, width = 7, height = 4.5, dpi = 300)

# residuals vs each numeric predictor
library(tidyr)
num_preds <- preds[sapply(mod_dat[preds], is.numeric)]
long <- cbind(mod_dat[as.integer(names(resid(fit))), num_preds, drop = FALSE],
              resid = resid(fit)) |>
  pivot_longer(-resid, names_to = "predictor", values_to = "value")

p_rp <- ggplot(long, aes(value, resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ predictor, scales = "free_x") +
  labs(x = "Predictor value", y = "Residuals", title = "Residuals vs each numeric predictor")
print(p_rp)
ggsave("resid_vs_predictors.png", p_rp, width = 9, height = 6, dpi = 300)

# most influential species (Cook's distance)
cd <- cooks.distance(fit)
head(sort(cd, decreasing = TRUE), 5)