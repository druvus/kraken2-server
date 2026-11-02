# Multi-database classification

Status: implemented in v0.2.0 (`--db` repeated, `--strict-merge`); see
"Verification" for how it is tested. The rest of this document is the
design as written before the implementation, kept because it records why
the server behaves as it does.

Design for classifying reads against several Kraken 2 databases in one server,
following what `k2 classify --db a,b` does in Kraken 2.17. This document records
the upstream behaviour as read from the submodule at commit 01fb1d9 (`scripts/k2`,
`src/merge.cc`, `src/libtax.cc`) and the choices made for the server.

## What upstream does

`k2 classify --db a,b,c reads.fq` runs four steps.

1. **Sanity check.** Every database directory must contain `nodes.dmp` (or
   `taxonomy/nodes.dmp`). A taxid that appears in several files must have the
   same parent everywhere, otherwise `k2` exits.
2. **Merged taxonomy.** The `seqid2taxid.map` files of all databases are
   concatenated. `libtax.generate_taxonomy` loads `names.dmp` and `nodes.dmp`
   from the *first* database, marks every taxid from the concatenated map, and
   writes a `taxo.k2d` containing the marked nodes and their ancestors. This
   is the same code path `kraken2-build` uses, so the merged taxonomy is a
   normal Kraken taxonomy whose internal ids differ from every input
   database's ids; external (NCBI) ids are the common currency.
3. **Per-database classification.** The reads are classified once per
   database with the ordinary classifier and the user's options (confidence,
   minimum hit groups, quality masking). Each run writes a standard kraken2
   output file.
4. **Pairwise merge.** The `merge` program combines two output files into one,
   repeatedly, until one remains. For each read, with the two lines aligned:
   - Both hit lists are parsed into runs of `taxid:count`. Markers `A`, `|:|`
     and `-:-` are kept as special taxids.
   - The runs are walked in lockstep by minimizer count. For each aligned
     segment the merged taxid is the LCA of the two external taxids in the
     merged taxonomy, where LCA(x, 0) = x. A marker segment keeps its marker.
   - Merged hit counts accumulate per merged taxid over non-marker segments.
     `total_minimizers` counts every segment including markers.
   - The call is recomputed with `resolve_tree` over the merged hit counts,
     with **confidence 0** (`k2` does not pass `-i`) and **no minimum hit
     groups**. The per-database C/U calls are never consulted.
   - The merged hit list is written with adjacent equal taxids collapsed.
   - Read counts for the report come from the merged calls. Distinct k-mer
     counts are only available if the per-database runs dumped counters.

Consequences worth stating:

- A read unclassified by every database on its own can be classified after
  the merge if its hit lists carry hits, because the hit-group filter is not
  reapplied. A read classified by one database can lose its call if hits in
  another database pull the LCA to root and the (zero) confidence is met by
  root, though that requires root itself to carry hits.
- The documented `--merge-policy first-classified|last-classified|resolve-lca`
  option is parsed by `k2` but not passed to `merge` at this commit; only
  the LCA behaviour exists.
- `merge` is a separate binary. The bioconda `kraken2` 2.17.2 package
  ships it under `share/kraken2-2.17.2/libexec`, where `k2` finds it.

## Server design

The server holds all databases in one process and has each read's minimizer
token stream in memory, so the file round trips disappear and the merge
becomes an element-wise operation.

### Loading

- `--db` may be repeated. Each database becomes an `Index` with its
  `IndexOptions`, `Taxonomy`, hash table and name (its directory basename).
- All databases must agree on `k`, `l`, `spaced_seed_mask`, `toggle_mask`,
  `dna_db` and `revcom_version`, since one token stream is shared. Otherwise
  the server logs the mismatch and marks the index broken.
  `minimum_acceptable_hash_value` may differ: it is applied per database at
  lookup time (a key below a database's threshold counts as no hit there).
- The indexes load in sequence in the loader thread; `ServerReady` stays
  false until all are loaded.
- With one `--db` nothing changes: the single-database path is a vector of
  one and the merged taxonomy *is* that database's taxonomy.

### Merged taxonomy

Instead of `names.dmp`, `nodes.dmp` and `seqid2taxid.map`, which many
databases distributed as `*.k2d` only do not carry, the merged taxonomy is
built from the input `taxo.k2d` taxonomies themselves:

- Collect every node of every database as (external id, parent external id,
  name, rank). A Kraken taxonomy already contains exactly the marked nodes
  and their ancestors, so the union equals what `generate_taxonomy` would
  produce from the concatenated maps.
- Consistency check, mirroring `k2`: an external id present in several
  databases must have the same parent external id in all of them. A conflict
  is an error at load time.
- Lay the union out as a Kraken taxonomy (children of a node in one
  contiguous block, parents before children) and write it to a temporary
  `taxo.k2d`, then load it as a `Taxonomy` and call
  `GenerateExternalToInternalIDMap()`. Reusing the file format keeps the
  report writer and the existing taxonomy code unchanged.
- For each database, precompute `to_merged[internal id] -> merged internal
  id` once, so the per-read path is an array lookup.

### Per read

1. Tokenise minimizers once (the existing pass 1).
2. For each database: one `GetBatch`, then a replay that yields a taxon
   vector in that database's internal ids, mapped through `to_merged`, with
   this database's `minimum_acceptable_hash_value` applied. Markers are
   pushed unchanged.
3. Merge element-wise across databases: for position *i*, the merged taxon
   is the LCA over databases of the merged-internal taxa (LCA(x, 0) = x);
   markers pass through. This is what `merge` computes with run-length
   alignment, without the alignment, because every database saw the same
   minimizers at the same positions.
4. Hit counts over merged taxa, `total_minimizers` = number of positions
   (markers included, as `merge` does), then `ResolveTree` with
   **confidence 0** and no hit-group filter, to match `k2` output. The
   user's `--confidence` and `--hit-groups` are applied only when a single
   database is configured, exactly as today.
5. Result: `tax_id` is the merged external id; the hit list is the merged
   taxon vector formatted by `AddHitlistString`, which already collapses
   runs; `name` stays unsent.
6. Report counters are kept per merged taxon, the same as the
   single-database path: with `--report-kmer`, `add_kmer` is called with the
   minimizer on a new lookup and `increaseKmerCount(1)` on a repeat, so the
   report gives real distinct counts rather than the sums `merge` produces
   without dumps. (Since v0.2.2 the distinct-count sketch is only kept with
   `--report-kmer`, which is the only report that shows it.)

### Options and compatibility

- `--db` repeated is the only new option. The protocol does not change; a
  client cannot tell how many databases a server holds. `GetSummary`
  reports over the merged taxonomy.
- `--confidence` and `--hit-groups` given with several databases print a
  note that they apply per database in `k2` but have no effect on the
  merged call, since `merge` discards the per-database calls and the server
  does not compute them.
- `--strict-merge` applies both thresholds to the merged call: the user's
  confidence in `ResolveTree` and the minimum hit-group count over merged
  minimizers. This filters the way a single-database run does and can
  differ from `k2`. It is off by default so that output stays comparable
  to `k2 classify --db a,b`.
- Memory: hash tables are additive; `--memory-mapping` applies to all.

## Verification

- Unit tests on the synthetic taxonomy: merged taxonomy from two databases
  with overlapping and disjoint leaves; conflict detection; element-wise
  LCA merge on hand-built taxon vectors including markers; call
  recomputation with hits spread across databases.
- Integration (`testing/parity_test.sh`, multi-database section): two tiny
  nucleotide databases sharing a species, built by `testing/make_multidb.sh`,
  are classified by the server and by `k2 classify --db a,b` from bioconda,
  single-end and paired. The per-read output and the report are identical.
  The section also checks that `--hit-groups` is ignored by default, as in
  `k2`, and applied with `--strict-merge`, that a single small database
  still matches `kraken2`, and that databases built with different settings
  are refused at load time.
- Single-database output remains byte-identical to `kraken2` (the rest of
  the parity script).
