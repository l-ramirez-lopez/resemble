# Classical memory-based learning (local modelling)

> *Think Globally, Fit Locally* – ([Saul and Roweis,
> 2003](#ref-saul2003think))

![](logo.png)

## 1 Memory-based learning

Memory-based learning (MBL) describes a family of local learning methods
that are well suited to complex and heterogeneous spectral datasets
([Ramirez-Lopez et al., 2013](#ref-ramirez2013spectrum)). In MBL,
instead of fitting a single global regression function, a local
regression model is fitted for each target observation using its nearest
neighbors identified in a calibration or reference set. Although the
global relationship between $X$ and $Y$ may be complex, MBL approximates
it through a collection of simpler local models, each assumed to be
valid within a restricted region of the predictor space ([Mitchell,
1997](#ref-mitchell1997machine)).

For a set of $m$ observations requiring prediction, this can be
expressed as

$${\widehat{y}}_{i} = {\widehat{f}}_{i}(x_{i};\theta_{i}),\qquad i = 1,\ldots,m$$

## References

Mitchell, T.M., 1997. Machine learning, volume 1 of 1.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Demattê, J.,
Scholten, T., 2013. The spectrum-based learner: A new local approach for
modeling soil vis–NIR spectra of complex datasets. Geoderma 195,
268–279.

Saul, L., Roweis, S., 2003. Think globally, fit locally: Unsupervised
learning of low dimensional manifolds. Journal of machine learning
research 4, 119–155.
