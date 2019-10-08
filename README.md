# pw85

Overlap test of two ellipsoids.

pw85 provides a C implementation of the “contact function” defined by
Perram and Wertheim (J. Comp. Phys. 58(3), 409–416,
[DOI:10.1016/0021-9991(85)90171-8][doi]) for two ellipsoids. Given two
ellipsoids, this function returns the *square* of the common factor by
which both ellipsoids must be scaled (their centers being fixed) in
order to be tangentially in contact.

pw85 is released under a BSD 3-Clause License.

- [documentation (HTML)][htmldoc]
- [documentation (PDF)][pdfdoc]

[doi]: https://doi.org/10.1016/0021-9991(85)90171-8 "Perram and Wertheim (1985)"
[htmldoc]: https://sbrisard.github.io/pw85/ "HTML documentation of the pw85 library"
[pdfdoc]: https://sbrisard.github.io/pw85/_downloads/7f04ad075e7d12a5ca017c9e5549f8f0/pw85.pdf "PDF documentation of the pw85 library"

## Release notes

- Version 1.0.1 (2019-10-08): new section “Overview” in the docs
- Version 1.0 (2019-09-25): first stable version
