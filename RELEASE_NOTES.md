# v0.2.0.dev0

Initial standalone development release, extracted on 2026-09-30 from the Maize Anther Primordia Spatial Small-RNA Analysis working tree.

The runtime package is copied unchanged. It requires explicit library format and length bounds, class-labelled FASTAs, fixed Phasis product TSVs, reference identities, and molecular orientation. DA contrasts and retained-depth normalization are explicit. The release adds a standalone README, MIT license text matching the package's existing license declaration, citation metadata, public installation/data-download instructions, and automated software checks.

The software manifest records the paths, SHA-256 values and version of the standalone files. The test suite includes the generic CLI and count-resumption checks. Project-specific audits and cluster jobs remain in the original project; biological inputs and result files are supplied separately.

GitHub Actions checks installation, package metadata, distribution builds and the retained software tests on synthetic inputs. Consult the repository's Actions page for the result at the exact commit you use. No full 18-library analysis was run for this publication. The real-data tutorial remains incomplete until its processed libraries and separately frozen Phasis products are accessible.

Dependencies follow `pyproject.toml`, including PyDESeq2 `>=0.5,<0.6`; this release is not a fully locked environment. Record the release commit and installed package versions with each analysis.
