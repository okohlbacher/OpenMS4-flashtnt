# Standalone FLASHTnT port candidate

This experimental package ports FLASHTnT from the BSD-3-Clause source at
`t0mdavid-m/OpenMS@3f508829ad81c91354d397966f28d428e29e5329` to installed,
exactly pinned OpenMS Core, CLI and FLASH SDKs. Native Release compilation,
metadata generation and port contracts passed against the installed ci.2 SDKs
on IBMI dax (Linux x64, GCC 14.4). Five-platform CI exercises the same package.
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
modifications and validates modification positions. Singleton modification
candidate lists become named modifications; ambiguous lists retain the observed
mass delta instead of inventing a single identification.

The port also replaces removed OpenMS String operations with standard strings,
removes an unused Qt include, gives the private classes ordinary static linkage,
and preserves the upstream normalization of tag sequences to uppercase.
It rejects missing/truncated deconvolution metadata before indexing peak arrays
and uses a sequential scan number when a spectrum has no native ID.
Open-ended extension uses an endpoint of -1. The port fixes an inherited
out-of-bounds read by inspecting the endpoint spectrum entry only when the
endpoint is specified. The AQPZ test exercises this path and reproduced the
invalid read under AddressSanitizer before the fix. Precursor estimation also
skips an exhausted complementary-path iterator before dereferencing it; the
same fixture exposed this second inherited out-of-bounds read. To check with GCC or
Clang, use a separate build with
`-DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer -g"`
and run the same CTest suite. Dependencies must use a compatible compiler/runtime.
The Linux x64 CI driver runs this additional instrumented build with address,
undefined-behavior and leak detection before it publishes the ordinary Release
payload. It retains the sanitizer logs and scientific comparison output.
`-DOPENMS4_WARNINGS_AS_ERRORS=ON` promotes package compiler warnings to errors.

Tag overlap checks explicitly use floating-point `std::abs`. In this port, the
unqualified call inherited from upstream selected integer `abs` with GCC 14.4 but floating-point `abs` with
Apple Clang 21, truncating sub-Dalton differences on Linux. The port contract
checks two peptide mass ladders separated by 0.5 Da to prevent this regression.

CTest also runs the bundled AQPZ example with FLASHApp's saved parameters. It
requires nonempty result tables, the expected strongest AQPZ identification and
full database/matched sequences, valid masses/positions and reconstructed AQPZ
tags. This is a biological identity and executable-interface check. It does
not establish complete numerical equivalence to the retained May 2025 cache:
that cache has no recorded binary revision, and the pinned upstream algorithms
changed subsequently, including the December 2025 inverse-tag-direction fix
(`8677e5015860a6ecce5e7bcc35295978040291b1`).

Before the memory and floating-point fixes, the four Linux tests passed in
10.47 seconds (AQPZ: 10.13 seconds). That revision produced 622 tags and
10 protein/PrSM rows; the historical cache has
2,968 tags and 17 rows. Its best AQPZ score/matching-fragment count was 559/79,
compared with 505/69 historically. Database and matched sequences and positions
1–240 agree, while inferred mass and coverage differ. Restoring the upstream
mutating uppercase operation removed an unintended port behavior change; the
remaining historical differences have not been causally attributed or approved
as scientifically equivalent.

At source revision `250debb`, three Mac ARM runs produced byte-identical result
tables with 698 tags (114 AQPZ tags), while repeated Linux x64 runs produced
622 tags (96 AQPZ tags). The strongest AQPZ sequence, score, mass, fragment count
and coverage agree across these platforms; tag counts do not. These observations
predate the memory and floating-point fixes.

At `93eb47c`, all four tests pass on Linux x64 (10.10 seconds), macOS ARM
(5.69 seconds), and Linux with ASan/UBSan and leak detection (36.97 seconds).
Linux now produces the same 698 tags (114 AQPZ tags), 10 protein/PrSM rows and
strongest identification fields as macOS. Complete historical numerical
equivalence remains unqualified, and the historical reference is unchanged.

The untouched historical counts, tag sequence multiplicities, mass, score and
fragment fields remain in `tests/data/aqpz/reference.json`. Every run compares
them and writes `aqpz-results/acceptance.json`; differences are not replaced with
new golden numbers. Require those additional historical comparisons explicitly:

```sh
OPENMS_TOOL_PREFIX_PATH=build python3 tests/test_aqpz.py \
  build/bin/FLASHTnT build/aqpz-reference --strict-reference
```

The input mzML, complete FASTA database and parameter file are copied unchanged
from public FLASHApp revision `f8e9eba435ea0843660c63fe86c58c64798e7f71`.
`source-provenance.json` records their original paths, Git objects and SHA-256
digests, plus the source result tables used for the compact reference. These
experimental checks do not yet qualify FLASHTnT as a released FLASHApp backend.
