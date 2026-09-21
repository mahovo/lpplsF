# CRAN-readiness review — lpplsF 0.1.0

Reviewed 2026-08-24 against `dev` @ `2124ca6`, in worktree
`lpplsf-cran-submission-0fa4c9`.

**Status: findings only, except B2, which has been implemented** (see §2 and §3).
Nothing has been committed — `dev` is untouched and every change below sits in the
worktree for review.

Baseline established fresh, not taken on trust: `R CMD build` + `R CMD check
--as-cran` on R 4.5.2, aarch64-apple-darwin20. Tarball **168 KB**, check **45 s**,
tests pass, **2 NOTEs**.

The second NOTE is *local toolchain only* — "'tidy' doesn't look like recent
enough HTML Tidy" and "package 'V8' unavailable". It will not appear on CRAN's
machines. So the handoff's "one NOTE" is right in substance.

One new signal the handoff didn't record, from `checking package dependencies`:

```
INFO: Vignette dependency required without any vignettes: 'knitr'
```

---

## 1. The vignette decision — both stated constraints are false

The decision was being blocked by two numbers. Both were measured, and both
turn out not to say what they appeared to say.

### "The cache is 7.4 MB, against CRAN's 5 MB guidance"

The 7.4 MB is **not data**. It is serialised `ggplot` objects. Measured, on disk,
across all 11 `.rds`:

| cache contents | gzip | bzip2 | xz |
|---|---|---|---|
| as-is (today) | **7.40 MB** | 4.88 MB | **1.55 MB** |
| drop `$fit_plot` only | 6.97 MB | 4.65 MB | 1.51 MB |
| drop all plot slots | **0.27 MB** | 0.23 MB | **0.17 MB** |

Two things follow.

**Re-saving with `compress = "xz"` and changing nothing else takes the cache from
7.40 MB to 1.55 MB — 4.8x, zero behavioural change, zero code change.** That
alone settles the size question.

**Dropping the stored plots takes it to 0.17 MB — 44x.** The numerical content
(`fit`, `mpl_output`, `log_price`, `fit_args`) is only ~700 KB in memory; the
rest is plot objects carrying their enclosing environments.

Caveat on that second row, found by checking rather than assuming: the document
*uses* the cached plot slots directly — `$fit_plot`, `$param_plot`,
`$trace_plot_mo/bm/bo`, `$mpl_plot` — and `plot.lppls_fit()` only *retrieves*
stored slots (`R/lppls_fitter.R:1384-1394`); it has no rebuild path. Only
`$fit_plot` can be rebuilt today, via the exported `fit_plot()`. So the 0.17 MB
row costs new rebuild paths for the param/trace/mpl plots. The 1.55 MB xz row
costs nothing.

Note also that `object.size()` is actively misleading here: it reports
`$fit_plot` as 15.9 MB of each fit's 18.7 MB, but dropping it saves only 0.43 MB
on disk, because it re-counts environments shared with the other plot slots. The
on-disk table above is the one to trust.

### "A full build takes ~20 minutes"

That is the `recompute = TRUE` path. The shipped path is `recompute = FALSE`,
which loads the cache. Measured end-to-end:

| render | time | output |
|---|---|---|
| as-is (`html_document` + plotly) | **11.5 s** | **14.2 MB** |
| as `rmarkdown::html_vignette`, static surfaces, dpi 96 | **12.5 s** | **2.72 MB** |

Both rendered clean — 41 figures, no chunk errors.

### The constraint that is actually real

Not the cache and not the time: the **rendered output**. 14.2 MB as it stands,
and that HTML is what `R CMD build` puts in the tarball.

Composition of the 14.2 MB: 7.66 MB of `<script>` (plotly.js, jQuery, Bootstrap)
+ 6.03 MB of base64 PNG. Switching to `html_vignette` and taking the static
`persp` branch the document *already carries* for LaTeX output
(`show_surface()`, `misc/reproduce_thesis_plots.Rmd:84-97`) removes the
JavaScript entirely: 2.72 MB, of which **95% is base64 PNG** — so the remainder
is tunable by `dpi`/`fig.width` alone.

Realistic tarball if the guided tour ships as a vignette:
`0.17 (pkg) + 1.55 (xz cache) + 2.72 (built HTML) ≈ 4.4 MB` — under the 5 MB
guidance, but not comfortably. With the plot slots dropped instead
(`0.17 + 0.17 + 2.72 ≈ 3.1 MB`), or with dpi tuning, comfortably under.

### Work the conversion needs, whichever shape is chosen

Found by reading the document, not estimated:

1. `pkgload::load_all(pkg_root)` (line 46) → `library(lpplsF)`. A vignette cannot
   `load_all`.
2. `lpplsF:::overlay_legend_theme(...)` (line 997) — a `:::` call into the
   package's own internals. Export it or inline it.
3. Undeclared packages used by the document: **`viridisLite`**, **`pkgload`**,
   **`rugarch`**. Any that survive into a vignette must go in `Suggests`, or the
   check raises "unstated dependencies in vignettes".
4. Add the `%\VignetteIndexEntry{}` / `%\VignetteEngine{}` block; drop the
   `pdf_document`/`html_document` output pair.

---

## 2. Findings, ranked

### Blockers

**B1 — `URL`/`BugReports` 404s. PARTIALLY FIXED ON `dev` — one still live.**

Originally both `DESCRIPTION:39-40` fields pointed at `github.com/yourusername/lpplsF`,
propagating into `man/lpplsF-package.Rd`, plus the README's BibTeX block.

`dev` commit **`d03fddb` "Fix url in DESCRIPTION and README"** (made after this
review's `2124ca6` baseline) fixed two of the three:

| location | state on `dev` |
|---|---|
| `DESCRIPTION` → `BugReports` | fixed → `mahovo/lpplsF/issues` |
| `README.md` BibTeX `url` | fixed → `mahovo/lpplsF` |
| `DESCRIPTION` → `URL` | **still `yourusername`** |

```
URL: https://github.com/yourusername/lpplsF          <- still 404
BugReports: https://github.com/mahovo/lpplsF/issues  <- fixed
```

So `R CMD check --as-cran` will still report one invalid URL, from `DESCRIPTION`
and `man/lpplsF-package.Rd`. One line left to change.

**B2 — `inst/extdata/SPX.csv` is CRSP data and cannot ship. CONFIRMED.**
14,727 daily observations, 1962-07-02 to 2020-12-31, column `spindx`. The user
confirmed the source from the thesis's "Data sources" section: **Wharton Research
Data Services, Center for Research in Security Prices, LLC (CRSP)**, CRSP Stock
Market Indexes, library `crspa`, variable `spindx` ("Level on S&P Composite
Index"), range 1925-12-31 to 2020-12-31, retrieved 2021-02-03.

Redistributing it in a public package would breach the WRDS/CRSP subscription
terms — an obligation independent of CRAN, and the one unambiguous, personally
binding exposure here. Thesis use is unaffected: that is exactly what the
academic licence covers, and the thesis documents it properly.

Three facts that shape the fix:

1. **`SPX.csv` is referenced nowhere in `R/`, `man/` or `tests/`** — only by
   `misc/` documents. No package code depends on it.
2. **Only rows 11800:14727 are ever used** — 2,928 daily closes, 2009-05-18 to
   2020-12-31. The other 11,799 rows are dead weight.
3. **The series is in the package twice.** `spx_fit.rds$log_price` embeds the same
   2,928 values; `exp(log_price)` reproduces `spindx` exactly (909.71, 908.13,
   903.47, …). Deleting `SPX.csv` alone would only make the CRSP data less
   obvious, not remove it. `fit_e008.rds`/`fit_e009.rds` also carry `log_price`,
   but those series are synthetic — no issue.

**Rebasing is lossless — verified.** The LPPLS design matrix carries an intercept
column, so an additive shift in log-price moves only `A`. Confirmed numerically:
`tc`, `m`, `omega` bit-identical; `A` absorbs the shift exactly; `B`, `C1`, `C2`
agree to ~1e-18; MPLE and all three likelihood intervals identical. So the package
*could* ship the series **rebased** (first close = 100) at no scientific cost, and
the shipped object would then demonstrably not be the S&P 500 index level.

**Not adopted.** The standard chosen for this package is CRAN precedent, which
`gt::sp500` and `rugarch::sp500ret` both meet with unmodified vendor levels plus a
source citation; rebasing exceeds that standard at the cost of shifting the SPX
figures' y-axis away from recognisable index levels, weakening the reproduction.
The shipped series is therefore the **raw Yahoo closes**. Rebasing remains a
one-line change if that trade is ever re-judged, with the guarantee above that no
reported quantity moves.

**THE SWAP IS VERIFIED SAFE.** Yahoo `^GSPC` was fetched via
`quantmod::getSymbols` for 2009-05-01 → 2020-12-31 and compared against the CRSP
window. The trading calendars align **exactly** — 2,928 dates in both, zero on
either side only. Levels agree to 1.6e-7 on average; the worst single day is
0.39 index points (2019-08-12, 2883.09 vs 2882.70), i.e. 1.35e-4 relative. In log
space, max diff 1.35e-4, sd 3.6e-6.

Re-running the vignette's own calibration (`fh = 500`, `hold_out = 500`,
`num_searches = 1`, MPL, cutoffs 0.05/0.10/0.50) on all three series:

| series | tc | m | omega | MPLE | critical date |
|---|---|---|---|---|---|
| CRSP (thesis) | 2717.00880442 | 0.9 | 6.15965123098 | **2717** | **2020-03-03** |
| Yahoo `^GSPC` | 2717.00874081 | 0.9 | 6.15964784842 | **2717** | **2020-03-03** |
| Yahoo rebased to 100 | 2717.00874083 | 0.9 | 6.15964784842 | **2717** | **2020-03-03** |

All three likelihood intervals are **bit-identical** across the three series:
LI(0.05) 2696–2757, LI(0.10) 2696–2757, LI(0.50) 2696–2744.

So every quantity the thesis and the vignette report survives the swap unchanged.
Only the Brent-refined `tc` moves, in its 7th decimal (2.3e-8 relative) — the same
character as the `factr` and objective-fusion changes already accepted.

**Noticed while verifying, unrelated to the swap:** this calibration is a
**boundary solution** — `m` sits exactly on its upper bound of 0.9, and the
package's own diagnostic warns that "the modified profile likelihood and its
likelihood intervals are unreliable for this fit". That warning fires on the CRSP
baseline too, so it is a property of the vignette's configuration, not of the data
source. Worth a decision separately from B2.

**Source landscape** (searched 2026-08-24): no authoritative permissively-licensed
*daily* S&P 500 series exists. FRED's `SP500` carries an explicit S&P Dow Jones
Indices notice — "reproduction … prohibited". datahub/`s-and-p-500` is ODC-PDDL
but **monthly** (Shiller-derived), which cannot resolve a 2,928-day window or the
March 2020 crash. Kaggle CC0 tags are uploader-declared and cannot convey rights
the uploader lacks. Yahoo `^GSPC` and Stooq are silent/ambiguous but are the de
facto norm, with heavy CRAN precedent: `gt::sp500` (Posit, daily 1950-2015, no
source stated), `rugarch::sp500ret` (Yahoo), `sp500SlidingWindow::SP500`,
`Lock5Data::SandP500`. Stooq now serves a JavaScript bot-detection challenge, so
it is a by-hand download, not a scriptable one.

#### Which CRAN precedent? `rugarch` and `tsgarch` are not the same standard

lpplsF already lists `tsgarch` in `Suggests`, so it is worth asking what that
package does. It ships two datasets — `dmbp` and `nikkei` — and **no S&P 500 data
at all**. Both come from **journal replication archives**, cited to the archive:

| package | dataset | stated source |
|---|---|---|
| `rugarch` | `sp500ret` | "Yahoo Finance" — the entire citation |
| `tsgarch` | `dmbp` | "Journal of Business & Economic Statistics Data Archive" |
| `tsgarch` | `nikkei` | "Journal of Applied Econometrics Data Archive, 2003-v18.6/giot-laurent", plus the full Giot & Laurent (2003), *JAE* 18(6), 641–664 citation |

Four things characterise the `tsgarch` standard beyond the source itself:

1. **Returns rather than levels** — both datasets are daily log returns in percent.
2. **A stated reason for inclusion** — "included as it is used for the GARCH
   benchmark" / "the APARCH benchmark". The dataset earns its place.
3. **Transformations disclosed in full.** The `nikkei` page spends a paragraph on a
   duplicate date at 2000-08-31, exactly how the rows were shifted, and why it does
   not affect estimation.
4. **Honest about limits** — `dmbp` records that "exact dates are not known as this
   dataset did not provide an index".

Both packages are by the same author, Alexios Galanos. The practice **tightened**
between them: the older `rugarch` cites a vendor in two words, the newer `tsgarch`
ships only journal-archive benchmark data with full provenance. So "CRAN precedent"
is not one standard but two, and the author chose the stricter one when he rebuilt.

lpplsF cannot match the `tsgarch` *source* — journal archives are tied to specific
published papers, and none carries daily SPX for 2009–2020. It can and does match
the *discipline*: precise source and retrieval date, a stated reason for inclusion,
and disclosure of what was done to the series. That lands above the `rugarch`/`gt`
precedent and as close to `tsgarch` as this window allows.

Note on form: point 1 does **not** transfer. GARCH models returns; LPPLS models the
log-price path, and `fit_lppls()` takes `log(close)` directly. The series therefore
ships as **levels**. Shipping returns was tested and is equally lossless — the
reconstruction `c(0, cumsum(diff(log_price)))` is exact to 0.000e+00 and gives the
same MPLE, LIs and critical date — but it would only add a reconstruction step.

#### B2 — RESOLVED, implemented in the worktree

| | before | after |
|---|---|---|
| dataset | `inst/extdata/SPX.csv`, CRSP, 232 KB, 14,727 rows | `data/spx.rda`, Yahoo `^GSPC`, **8.7 KB**, 2,928 rows |
| second copy | `spx_fit.rds$log_price`, CRSP, 2,274 KB | regenerated from Yahoo, xz, **455 KB** |
| documentation | none | `R/data.R` → `man/spx.Rd`, `tsgarch` style |

* `inst/extdata/SPX.csv` deleted; `inst/` removed (it held nothing else).
* `data/spx.rda` added — `date` (Date) and `close` (numeric), 2009-05-18 to
  2020-12-31, exactly the calibrated window. This also makes the existing
  `LazyData: true` correct, closing **M2**.
* `R/data.R` documents it in the `tsgarch` manner: source and retrieval date,
  reason for inclusion, the CRSP substitution and its measured effect.
* `spx_fit.rds` regenerated with the vignette's `spx-compute` config verbatim.
  Verified: **MPLE 2717 → 2020-03-03, m = 0.9, all three LIs bit-identical** to the
  CRSP fit; only Brent-refined `tc` differs (7th decimal).
* Load sites updated in `misc/reproduce_thesis_plots.Rmd`, `misc/performance.Rmd`
  and `misc/lagrange_inference_report.Rmd`; a provenance note added to the vignette.
* Vignette re-rendered on the cached path: 11.2 s, exit 0, **0 chunk errors**, 35
  figures, `2020-03-03` present.

**Not resolved by this, and worth deciding:** the CRSP file remains in the repo's
**git history**. Removing it from `HEAD` does not remove it from earlier commits,
and if `mahovo/lpplsF` is or becomes public, the data stays reachable there. Only
a history rewrite removes it — a call only the repository owner can make.

**B3 — `VignetteBuilder: knitr` with no `vignettes/` directory.** Produces both
the NOTE line and the INFO line. Resolved either way by §1's decision: add
vignettes, or drop the field *and* drop `knitr`/`rmarkdown` from `Suggests`
(they are otherwise unused — 0 uses in `R/` and `tests/`).

### High — a reviewer will almost certainly raise these

**H1 — `\dontrun{}` on six help pages, including the flagship `fit_lppls()`.**
`fit_lppls`, `lppls_rolling`, `lppls_lagrange`, `lr_trend_test`,
`contour_plot_sse`, `fit_plot`. CRAN's standard reviewer response is that
`\dontrun{}` is for examples that genuinely *cannot* run, and to replace it with
`\donttest{}` or unwrap it.

The `fit_lppls()` example, run exactly as written: **1.18 s**. It is well inside
CRAN's limits and has no reason to be unexecuted. Consequence today: every
example timing in the check is ~0, i.e. **no exported function's example is
actually being exercised by `R CMD check`.**

**H2 — `fit_lppls()` changes the user's random number stream.**
`R/lppls_fitter.R:625` (`set.seed(i)` per search) and `:917` (`set.seed(fit2_best$ID)`
for the trace replay). Verified empirically: `.Random.seed` differs across a
call, and the user's next `runif(1)` returns 0.3890714 instead of 0.3849424.
CRAN policy is that packages must not alter the user's session state.

The fix that keeps behaviour bit-identical is to save and restore around the
function, leaving the internal `set.seed(i)` calls exactly as they are:

```r
if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
  .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
} else {
  on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
}
```

Same applies to `.sim_ar_garch()` (`R/lr_trend_test.R:16`), though there the seed
is user-supplied and opt-in.

Worth documenting either way: the hard-coded `set.seed(i)` means the "random"
restarts are a fixed set, so `num_searches` is deterministic across calls and
does not respond to the user's own `set.seed()`. That is a deliberate
reproducibility choice, but it is currently undocumented and will read as a bug.

**H3 — the two method references carry no DOIs.** CRAN asks for
`Authors (year) <doi:...>` in `Description`. Verified by resolving each against
doi.org:

- `<doi:10.1016/j.physa.2013.04.012>` — Filimonov & Sornette, Physica A 392(17), 3698-3707
- `<doi:10.1080/14697688.2016.1276298>` — Filimonov, Demos & Sornette, Quant. Finance 17(8), 1167-1186

(Beware: a web search returned `...04.013` for the first. That DOI is a different
paper — "Geometry of complex networks and topological centrality". Resolved both
before recording them.)

**H4 — `Depends: R (>= 4.0.0)` is wrong.** Measured against current CRAN:
`symengine` needs R >= 4.2; `ggplot2` >= 4.1; `dplyr` >= 4.1.0; `tidyr` >= 4.1.0.
The true floor is **R >= 4.2**. As declared, the package cannot install on the R
versions it claims to support.

### Medium

**M1 — three examples are not self-contained.** `man/fit_plot.Rd`,
`man/contour_plot_sse.Rd` and `man/lr_trend_test.Rd` use `log_p` / `t` without
ever defining them, so they would error if run. This is *why* they are in
`\dontrun{}`; unwrapping per H1 requires fixing them first.
`man/contour_plot_sse.Rd` also names a local `par`, shadowing `graphics::par`.

**M2 — `LazyData: true` with no `data/` directory. CLOSED by the B2 work.**
`R CMD build` was silently dropping it ("Omitted 'LazyData' from DESCRIPTION"),
which is why no NOTE appeared. `data/spx.rda` now exists, so the field is correct
as it stands — and this is also the first piece of the `tsbs` layout the package
was always going to need.

**M3 — `README.md:105` links to `DEVELOPMENT.md`, which does not exist.** The
README ships in the tarball, so this is a dead link in the published package.

**M4 — no `inst/CITATION`.** The README instead carries a hand-written BibTeX
block with the 404 URL. A `CITATION` file would make `citation("lpplsF")` correct
and is conventional for a package backing a thesis.

**M5 — the flagship example produces a degenerate fit.** Running
`man/fit_lppls.Rd`'s example emits `1000 of 1000 estimates had D out of range`.
The synthetic parameters put every estimate outside the damping filter. It works,
but it advertises the package badly, and it will be the first thing a reviewer
runs once H1 is fixed.

**M6 — citation year wrong.** `man/lppls_lagrange.Rd` gives "Demos, G., &
Sornette, D. (2019) ... Quantitative Finance, 17(5), 657-675". Volume 17(5) is
2017 (`doi:10.1080/14697688.2016.1231417`, resolved).

**M7 — no `cran-comments.md`, no `NEWS.md`.** `cran-comments.md` is effectively
expected at submission; `NEWS.md` is conventional from 0.1.0 onward.

**M8 — `.Rbuildignore:11-12` has unescaped dots**: `^.git$`, `^.gitignore$`.
Harmless in practice (both are excluded by default anyway) but the patterns are
regexes and do not mean what they appear to.

### Low

- **L1** — `importFrom(dplyr, filter)` is unused in `R/lppls_fitter.R` after the
  performance work. `dplyr::filter` is still used qualified in `R/plotting.R`, so
  the dependency itself stays. (Carried over from the handoff; confirmed.)
- **L2** — the HTML-manual NOTE is local only (old HTML Tidy, no V8). Nothing to fix.
- **L3** — `misc/AI_bubble.Rmd` still passes `factr = 1e-08`. `misc/` is
  `.Rbuildignore`d so it never reaches CRAN. Per the handoff, leave
  `misc/lpplsF_equivalence_checks.Rmd:197,206` alone.

---

## 3. What is already in good shape

Worth recording, since "assume nothing has been checked" cuts both ways:

- All 16 declared dependencies are current on CRAN; none archived.
- `\value` is present and meaningful on all 18 function help pages.
- No writing to the user's filespace; no `setwd`, `sink`, `library()`/`require()`
  in package code; no `:::` in `R/`; no `T`/`F` for `TRUE`/`FALSE`.
- `graphics::par()` is correctly saved and restored under `on.exit`
  (`R/sse_surface_plot.R:25-26`).
- All `cat()` calls are inside `print`/`summary` methods, where they belong.
- `Authors@R` is well-formed with an ORCID; `MIT + file LICENSE` is correct.
- `.DS_Store` in `tests/testthat/` does *not* reach the tarball — `R CMD build`
  excludes it.
- Tests pass in the built tarball.
