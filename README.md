# gaussian-composites
Functions for distributions which are Gaussian mixtures or quantile interpolations


There are two types of distributions described in this demo.  

The first type of distribution is a mixture.  To generate a sample of a mixture, choose one of the mixture components at random, then sample that distribution.  This library handles Gaussian mixtures.  

The second type of distribution is defined by interpolating quantiles.  Let $f$ be a generalized average, where for a vector input $x$, we have $\min(x)\leq f(x)$ and  $f(x) \leq \max(x)$.  Furthermore let $f$ be monotonic in every input.  If $Q_i$ are the quantile functions for some set of distributions, then we can define a new quantile function $Q: p \mapsto f( Q_1(p), \dots, Q_n(p))$.  This is monotonic since $f$ and $Q_i$ are monotonic, and its range is in $[0,1]$ becuase $f$ is bound by the $Q_i$ whose range is $[0,1]$.  

This library handles distributions whose quantiles are means of Gaussian quantiles. In other words, $f$ is a weighted mean and all of the $Q_i$ are Gaussian quantiles.

Graphically, the mixture corresponds to averaging distribution functions vertically  The $q$-interpolation corresponds to averaging distributions horizontally.  
