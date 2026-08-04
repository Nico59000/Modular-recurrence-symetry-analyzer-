# Modular Recurrence Analyzer

Analyseur autonome en Python des périodes et symétries de suites modulo
un entier, avec spécialisation des diagnostics à \(\mathbb F_3\).

## Prérequis

- Python 3.10 ou plus récent ;
- aucune dépendance externe.

## Conventions

Les coefficients sont fournis du terme le plus ancien au plus récent :

```text
u[n+k] = c[0]u[n] + ... + c[k-1]u[n+k-1]  (mod m)
```

## Exemples

### Fibonacci modulo 3

```bash
python modular_recurrence_analyzer.py recurrence \
  --coefficients 1,1 \
  --seed 0,1 \
  --modulus 3 \
  --name fibonacci
```

Le résultat contient :

```text
period_word = 01120221
global_antiperiod_h = 4
half_antiperiodic = true
reversed_second_half = 1220
```

### Analyser un mot périodique arbitraire

```bash
python modular_recurrence_analyzer.py word \
  --word 01120221 \
  --modulus 3
```

### Énumérer tous les cycles d'une récurrence

```bash
python modular_recurrence_analyzer.py cycles \
  --coefficients 1,1,1,0 \
  --modulus 3 \
  --name tritetranacci \
  --json-out tritetranacci.json
```

### Reproduire les quinze familles du dépôt étudié

```bash
python modular_recurrence_analyzer.py legacy \
  --modulus 3 \
  --json-out legacy.json \
  --csv-out legacy.csv
```

### Balayer tous les vecteurs de coefficients binaires d'ordres 2 à 6

```bash
python modular_recurrence_analyzer.py scan \
  --min-order 2 \
  --max-order 6 \
  --modulus 3 \
  --workers 4 \
  --json-out binary_families.json
```

## Sous-commandes

- `word` : toute période explicite, sans hypothèse de récurrence ;
- `recurrence` : une récurrence et une graine, via Brent ;
- `cycles` : tous les cycles d'une récurrence ;
- `legacy` : les quinze familles nommées dans le dépôt ;
- `scan` : tous les masques binaires non nuls dans une plage d'ordres.

## Diagnostics produits

- prépériode et période d'état ;
- mot périodique primitif ;
- opposition modulo \(m\) ;
- antipériodicité de demi-période ;
- paires de cycles opposés ;
- renversement de la seconde moitié ;
- symétries affines et diédriques ;
- test spectral pair/impair modulo 3 ;
- recherche de \(M^h=-I\).

## Tests

```bash
python -m unittest -v test_modular_recurrence_analyzer.py
```

## Limites

La commande `cycles` explore exactement \(m^k\) états et impose donc une
limite configurable. Pour les ordres beaucoup plus grands, une méthode
par factorisation du polynôme caractéristique sur un corps fini est plus
adaptée que l'énumération exhaustive.

# Candidate extraction and OEIS comparison

Two additional scripts extend the modular-period analyzer.

## 1. Generate integer candidates from period windows

The construction used by OEIS A276275 is reproduced by:

```bash
python pisano_seed_sequence_generator.py single \
  --period 1112201210010 \
  --coefficients 1,1,0 \
  --modulus 3 \
  --family padovan \
  --start 4 \
  --one-based \
  --terms 48 \
  --json-out A276275_reproduction.json
```

The selected window is

```text
period: 1112201210010
start:  4
order:  3
seed:   220
```

and the lifted sequence begins

```text
2, 2, 0, 4, 2, 4, 6, 6, 10, 12, 16, 22, ...
```

To generate one candidate for every cyclic order-length window:

```bash
python pisano_seed_sequence_generator.py all \
  --period 1112201210010 \
  --coefficients 1,1,0 \
  --modulus 3 \
  --family padovan \
  --terms 50 \
  --deduplicate seed \
  --json-out padovan_period_candidates.json \
  --csv-out padovan_period_candidates.csv
```

The coefficient convention is oldest-to-newest.  Therefore Padovan uses
`1,1,0`, while the OEIS/Mathematica `LinearRecurrence` signature is
reported as `0,1,1`.

## 2. Search and compare with OEIS

The OEIS provides read-only JSON search results, but publication remains a
manual authenticated workflow.

Online comparison:

```bash
python oeis_sequence_compare.py \
  --candidates padovan_period_candidates.json \
  --mode online \
  --query-terms 10 \
  --delay 2.0 \
  --cache-dir .oeis_cache \
  --json-out oeis_matches.json \
  --csv-out oeis_matches.csv \
  --draft-dir review_drafts
```

Offline comparison against the official `stripped.gz` dataset:

```bash
python oeis_sequence_compare.py \
  --candidates padovan_period_candidates.json \
  --mode offline \
  --stripped stripped.gz \
  --min-match 8 \
  --json-out local_matches.json
```

The comparator is deliberately read-only.  It does not log in, edit, or
submit entries.  Draft files are review aids for manual submission only.
