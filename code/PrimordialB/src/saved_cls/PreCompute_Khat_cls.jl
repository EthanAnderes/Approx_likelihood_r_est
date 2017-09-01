#= ###################################################

* This pre-loads and saves the radial profiles of the bessel trasformations used in "likewhite.jl".
* These are used in `YXcov(x1, x2)`.
* This also saves the cls.
* The variables can be uploaded with
```
using JLD
@load "src/saved_cls/vKXnatxsq.jld"
@load "src/saved_cls/cls.jld"
```

!!!!!!!!!!
this needs updating
!!!!!!!!!


=# ################################################

using PrimordialB
using JLD

cls  = PrimordialB.class(
    r           = 0.2,
    r0          = 1.0,
    lmax        = 20_000,
    omega_b     = 0.0224567,
    omega_cdm   = 0.118489,
    tau_reio    = 0.128312,
    theta_s     = 0.0104098,
    logA_s_1010 = 3.29056,
    n_s         = 0.968602
)



# --- compute the radial profile
vKE2atxsq, vKB2atxsq, ugrid = PrimordialB.vecKn_atxsq(2, cls)
vKE3atxsq, vKB3atxsq, ugrid = PrimordialB.vecKn_atxsq(3, cls)
vKE4atxsq, vKB4atxsq, ugrid = PrimordialB.vecKn_atxsq(4, cls)
vKE5atxsq, vKB5atxsq, ugrid = PrimordialB.vecKn_atxsq(5, cls)
vKE6atxsq, vKB6atxsq, ugrid = PrimordialB.vecKn_atxsq(6, cls)

# --- save
@save "src/saved_cls/vKXnatxsq.jld" ugrid vKE2atxsq vKB2atxsq vKE3atxsq vKB3atxsq vKE4atxsq vKB4atxsq vKE5atxsq vKB5atxsq vKE6atxsq vKB6atxsq

# --- also save the cls
@save "src/saved_cls/cls.jld" cls
