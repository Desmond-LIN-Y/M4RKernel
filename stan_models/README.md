# Stan models

This directory contains all of the Stan code we use to run simulation experiments, cross-sectional analyses, and longitudinal analyses.

| Naming            | Description                                                                       |
|------------------ |-----------------------------------------------------                              |
| gp-functions.stan | Contains general helper functions used in all HSGP models                         |
| 2d-functions.stan | 2D version HSGP general helper functions in all 2d HSGP models                    |
| 2d-functions.stan | 2D vectorized version HSGP general helper functions in all 2d HSGP models         |
| hsgp-             | 1D Hilbert space approximate Gaussian process models                              |
| 2dhsgp-           | 2D Hilbert space approximate Gaussian process models                              |
| HIV-              | RCCS variations                                                                   |
| eq-               | Models that use the exponential quadratic (squared exponential) covariance kernel |
| m12-              | Models that use the Matern 1/2 covariance kernel                                  |
| m32-              | Models that use the Matern 3/2 covariance kernel                                  |
| m52-              | Models that use the Matern 5/2 covariance kernel                                  |
| -cd               | Classical parameterisation (age-age) of the contact rate matrices                 |
| -rd               | Restructured parameterisation (difference-in-age) of the contact rate matrices    |
| -X- , -N-         | Product and sum kernels                                                           |
