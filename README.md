# Standalone FLASHTnT port candidate

This experimental package ports FLASHTnT from the BSD-3-Clause source at
`t0mdavid-m/OpenMS@3f508829ad81c91354d397966f28d428e29e5329` to installed,
exactly pinned OpenMS Core, CLI and FLASH SDKs. Native qualification is pending.
The source provenance records the original Git objects, including the extracted
Tag/DAG helper definitions. No Core, FLASH, third-party or GUI source is built.

Configure with an installed matching SDK chain in `CMAKE_PREFIX_PATH`, then build
and run CTest. Source archives must set `OPENMS4_SOURCE_REVISION` and
`OPENMS4_SOURCE_DIRTY` explicitly; ordinary Git checkouts derive their identity.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/sdk
cmake --build build --parallel 16
ctest --test-dir build --output-on-failure
cmake --install build --prefix /path/to/install
```

The public executable interface remains `FLASHTnT -in out_deconv.mzML -fasta
database.fasta -out_tag tags.tsv -out_pro protein.tsv -out_prsm prsms.tsv`.
`-write_ini` and `-write_ctd` expose `tag:*`, `ex:*`, ion types and FDR parameters
used by FLASHApp. Input mzML must carry FLASHDeconv's `DeconvMassInfo` metadata;
ordinary raw/profile mzML is not tagging input.

Algorithms and Tag/DAG records are private to this executable. The local
FLASHTnTFile writer renders this new upstream tool's result records; the Core
SDK retains all existing readers/writers and supplies mzML, FASTA and ProForma.
The writer now uses Core's ProForma AST instead of the incompatible legacy
ProForma class. This corrects the upstream writer's nested brackets for range
modifications and validates modification positions. It requires native contract
tests and representative scientific output checks before release.

The port also replaces removed OpenMS String operations with standard strings,
removes an unused Qt include, gives the private classes ordinary static linkage,
and stores uppercase tag sequences without mutating their original case.
It rejects missing/truncated deconvolution metadata before indexing peak arrays
and uses a sequential scan number when a spectrum has no native ID.
No compatible FLASHTnT artifact or scientific acceptance is claimed yet.
