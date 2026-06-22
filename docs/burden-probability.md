# Burden tail-probability model

This document gives the exact computation used by [`scripts/QTLBurden.R`](../scripts/QTLBurden.R) to produce the burden tail-probability columns.

## Per-individual, per-gene model

For one gene and one individual:

- Let `β_i` be the aFC effect (log2 scale) for variant `i`.
- Let `SE_i` be its standard error.
- Let `d_i` be genotype dosage (`0`, `1`, `2`; missing values are treated as `0` for the burden score and tail-probability calculations).

The model computes:

```text
L = Σ_i d_i * β_i
```

and assumes variant effects are independent normals:

```text
β_i ~ Normal(β_i, SE_i^2)
```

Then:

```text
L ~ Normal(μ_L, σ_L^2)
μ_L = Σ_i d_i * β_i
σ_L^2 = Σ_i d_i^2 * SE_i^2
```

### Intuition

Think of each individual’s genotype dosages as a fixed **dosage-like perturbation profile** (the weights `d_i`) applied to the variant effect sizes.

If you could repeatedly draw effect sizes from their uncertainty distributions (`β_i ~ Normal(β_i, SE_i^2)`), the gene burden `L` varies around `μ_L`.

Then:

- `burden_probability_loss50` is the probability this perturbation shifts expression at least as low as the loss cutoff.
- `burden_probability_gain50` is the probability this perturbation shifts expression at least as high as the gain cutoff.

If no standard error is available for a variant (`SE_i` missing), `SE_i` is set to `0` for that variant.

`QTLBurden.R` also accepts 95% CI bounds directly:

```text
SE_i = (log2_aFC_upper - log2_aFC_lower) / (2 * 1.96)
```

This is applied in `auto` mode when `VariantEffectSEColumn` is not explicitly set and no direct SE column is present.

## One-sided tail probabilities

The log2 thresholds are:

- `loss_threshold = log2(0.5)` by default (`-0.5849625007`)
- `gain_threshold = log2(1.5)` by default (`0.5849625007`)

Then:

```text
P(loss) = P(L <= loss_threshold)
       = Φ((loss_threshold - μ_L) / σ_L)

P(gain) = P(L >= gain_threshold)
       = 1 - Φ((gain_threshold - μ_L) / σ_L)
```

where `Φ` is the standard normal CDF.

If `σ_L` is effectively zero, the code returns a hard 0/1 tail score:

- `P(loss) = I(μ_L <= loss_threshold)`
- `P(gain) = I(μ_L >= gain_threshold)`

## Output fields

From one burden run, for each row (sample, gene):

- `burden_probability_loss50` is `P(loss)`.
- `burden_probability_gain50` is `P(gain)`.
- `burden_probability` is direction-specific:
  - if `BurdenDirection = loss` / `deletion`, it equals the loss tail probability `P(loss)`
  - if `BurdenDirection = gain` / `duplication`, it equals the gain tail probability `P(gain)`
  - if `BurdenDirection = both`, it is currently set to the loss tail probability (`P(loss)`)

In the cleaning step, quality flags and downstream weighted summaries use:

- `is_noisy_extreme_call` with `BurdenTailProbability` cutoff.
- `burden_tail_weight` for weighted medians (loss-tail probability for negative-extreme bins, gain-tail probability for positive-extreme bins, and `1` otherwise).

## Interpretation

These are analytic tail probabilities under an approximate normal model, not empirical permutation frequencies or bootstrap estimates.
They can be read as:

```text
P(L ≥ gain_threshold | {β_i, SE_i}, dosage profile)  [for gain-like direction]
P(L ≤ loss_threshold | {β_i, SE_i}, dosage profile)  [for loss-like direction]
```

In words: the chance of observing an extreme expression effect (at or beyond the configured threshold) for that sample-gene pair, under uncertainty in effect sizes but with the observed dosage profile fixed.

`BurdenProbability` values can be interpreted as tail-based confidence for extreme expression change, and are suitable for continuous/weighted summaries, flagging, or filtering.
