# synon

**Taxonomic translations to a custom checklist**

`synon` is an R package for standardizing vascular plant scientific names in large, messy datasets—especially herbarium specimen data—by translating synonyms to accepted names using a prioritized, transparent workflow.

It is designed for real-world taxonomic chaos: outdated names, inconsistent infraspecific ranks, regional checklists, and spelling variation.

---

## ✨ Key Features

- Works with **large datasets** (tens to hundreds of thousands of names)
- Supports **custom target checklists**
- Uses a **hierarchical translation strategy**:
  1. Direct matches to the checklist
  2. User-supplied synonym tables
  3. Built-in synonym tables
  4. Optional exact and fuzzy matching via **WCVP**
- Preserves **provenance** of every translation
- Never overwrites an already-resolved name
- Optional handling of **infraspecific names** via binomial fallback
- Caches expensive fuzzy matching results for reuse

---

## 📦 Built-in Synonym Sources

`synon` ships with curated synonym lookup tables from:

- **NatureServe Network (USA/Canada)**
- **SEINet**
- **USDA PLANTS**
- **World Checklist of Vascular Plants (WCVP)**

These tables are filtered dynamically to your target checklist to minimize spurious translations.

---

## 🔧 Installation

```r
# install.packages("devtools")
devtools::install_github("/synon")
