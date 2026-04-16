# Overview of the functions in the resemble package

**maturing lifecycle**

Functions for spectral dissimilarity assessment, nearest-neighbour
search, memory-based learning, local expert libraries, and evolutionary
training subset selection in spectral chemometrics.

## Details

This is the version 3.0.0 – vortex of the package. It implements a
number of functions useful for modeling complex spectral spectra (e.g.
NIR, IR). The package includes functions for dimensionality reduction,
computing spectral dissimilarity matrices, nearest neighbor search, and
modeling spectral data using memory-based learning. This package builds
upon the methods presented in Ramirez-Lopez et al. (2013a)
[doi:10.1016/j.geoderma.2012.12.014](https://doi.org/10.1016/j.geoderma.2012.12.014)
, Ramirez-Lopez et al. (2026a) and Ramirez-Lopez et al. (2026b).

Development versions can be found in the github repository of the
package at <https://github.com/l-ramirez-lopez/resemble>.

The functions available for computing dissimilarity matrices are:

- [`dissimilarity`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md):
  Computes a dissimilarity matrix based on a specified method.

- [`diss_pca`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md):
  constructor for principal components-based dissimilarity method.

- [`diss_pls`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md):
  constructor for partial least squares-based dissimilarity method.

- [`diss_correlation`](https://l-ramirez-lopez.github.io/resemble/reference/diss_correlation.md):
  constructor for correlation-based dissimilarity method.

- [`diss_euclidean`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md):
  constructor for euclidean distance-based dissimilarity method.

- [`diss_mahalanobis`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md):
  constructor for Mahalanobis distance-based dissimilarity method.

- [`diss_cosine`](https://l-ramirez-lopez.github.io/resemble/reference/diss_cosine.md):
  constructor for cosine-based dissimilarity method.

The functions available for evaluating dissimilarity matrices are:

- [`diss_evaluate`](https://l-ramirez-lopez.github.io/resemble/reference/diss_evaluate.md):
  Evaluates the effectiveness of a dissimilarity matrix using side
  information.

The functions available for nearest neighbor search:

- [`search_neighbors`](https://l-ramirez-lopez.github.io/resemble/reference/search_neighbors.md):
  Search for nearest neighbors of a query spectrum in a reference
  dataset based on a specified dissimilarity method.

The functions available for modeling spectral data:

- [`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md):
  Memory-based learning for modeling spectral data.

- [`gesearch`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md):
  An evolutionary method to search optimal samples in large spectral
  datasets.

- [`liblex`](https://l-ramirez-lopez.github.io/resemble/reference/liblex.md):
  Builds a library of reusable localized models.

The functions available for dimensionality reduction are:

- [`ortho_projection`](https://l-ramirez-lopez.github.io/resemble/reference/ortho_projection.md):
  Computes an orthogonal projection of the data based on either
  principal components or partial least squares.

- [`predict.ortho_projection`](https://l-ramirez-lopez.github.io/resemble/reference/ortho_projection.md)

## References

Ramirez-Lopez, L., Viscarra Rossel, R., Behrens, T., Orellano, C.,
Perez-Fernandez, E., Kooijman, L., Wadoux, A. M. J.-C., Breure, T.,
Summerauer, L., Safanelli, J. L., & Plans, M. (2026a). When spectral
libraries are too complex to search: Evolutionary subset selection for
domain-adaptive calibration. *Analytica Chimica Acta*, under review.

Ramirez-Lopez, L., Metz, M., Lesnoff, M., Orellano, C., Perez-Fernandez,
E., Plans, M., Breure, T., Behrens, T., Viscarra Rossel, R., & Peng, Y.
(2026b). Rethinking local spectral modelling: From per-query refitting
to model libraries. *Analytica Chimica Acta*, under review.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte,
J.A.M., Scholten, T. (2013a). The spectrum-based learner: A new local
approach for modeling soil vis-NIR spectra of complex data sets.
*Geoderma* 195-196, 268-279.

## See also

Useful links:

- <https://github.com/l-ramirez-lopez/resemble>

- Report bugs at <https://github.com/l-ramirez-lopez/resemble/issues>

## Author

**Maintainer / Creator**: Leonardo Ramirez-Lopez
<ramirez.lopez.leo@gmail.com>

Authors:

- Leonardo Ramirez-Lopez
  ([ORCID](https://orcid.org/0000-0002-5369-5120))

- Antoine Stevens ([ORCID](https://orcid.org/0000-0002-1588-7519))

- Claudio Orellano
