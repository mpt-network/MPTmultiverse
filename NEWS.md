
# MPTmultiverse 0.4-0

- Added a `NEWS.md` file to track changes to the package.
- Added within-subjects comparisons of parameter estimates that are stored in column
    `test_within` of the results object.
    - For the Bayesian methods, these comparisons rely on posterior distributions
        of parameter differences.
    - For the no-pooling maximum-likelihood methods, comparisons rely on paired *t* tests.
    - For the complete-pooling maximum-likelihood method, comparisons rely on the the point estimates
       of parameters and the estimates of the standard errors of differences calculated
       from the Hessian matrix.
