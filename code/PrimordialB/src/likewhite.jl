

###############################

# patch extraction

#############################################


function extract_patchs(dqx, dux, x, y, dx, dy, σmod, sizbl)
    lengbl = sizbl[1]*sizbl[2]
    gnside = size(dqx,1)
    x1rng  = repmat(collect(1:gnside).',gnside,1)
    x2rng  = x1rng.'

    lnx1 = x + dx
    lnx2 = y + dy

    dqxblks  = Vector{Float64}[]
    duxblks  = Vector{Float64}[]
    lnx1blks = Vector{Float64}[]
    lnx2blks = Vector{Float64}[]
    x1blks = Vector{Float64}[]
    x2blks = Vector{Float64}[]
    σmodblks = Vector{Float64}[]
    iblocks = Vector{Int64}[] # the linear indices of where the local blocks are located withing (x,y)
    for ib in 1:sizbl[1]:gnside, jb in 1:sizbl[2]:gnside
        iblock = find((ib .<= x1rng .< (ib+sizbl[1])) & (jb .<= x2rng .< (jb+sizbl[2])))
        if !isempty(iblock)
            push!(dqxblks,   dqx[iblock])
            push!(duxblks,   dux[iblock])
            push!(lnx1blks, lnx1[iblock])
            push!(lnx2blks, lnx2[iblock])
            push!(x1blks, x[iblock])
            push!(x2blks, y[iblock])
            push!(σmodblks, σmod[iblock])
            push!(iblocks, iblock)
        end
    end
    return lnx1blks, lnx2blks, dqxblks, duxblks, x1blks, x2blks, σmodblks, iblocks
end




function extract_patch_indices(size_vec, sizbl)
    lengbl = sizbl[1]*sizbl[2]
    gnside = size_vec
    x2blks = Vector{Float64}[]
    order = 1:size_vec
    if rem(size_vec+1,lengbl) == 0
        blk_boundaries = collect(1:lengbl:size_vec+1)
    else
        blk_boundaries = vcat(1:lengbl:size_vec+1, size_vec+1)
    end
    block_indices = Vector{Int64}[]
    for k = 2:length(blk_boundaries)
        iblock = order[blk_boundaries[k-1]:blk_boundaries[k]-1]
        if !isempty(iblock)
            push!(block_indices,  iblock)
        end
    end
    return block_indices
end




#######################################################

# XYcov constructs all the covariance matrices

#######################################################



function scalebesselj(n::Int, u::Array{Float64,1})
    rtn = zeros(u)
    @inbounds @simd for i in eachindex(rtn)
        rtn[i] = besselj(n, u[i]) / (u[i] ^ n)
    end
    rtn[u .== 0.0] = (0.5 ^ n) * (1 / gamma(n+1))
    return rtn
end

# Note: this computes Khat(u^2)
function vecϕKn{dm, T}(n, cXXk, ugrid, g::FFTgrid{dm, T})
    #ugrid   = linspace(0.0, √(g.period), 5000) .^ 3
    #ugrid   = vcat(0.0, ugrid)
    # pre-compute the radial power
    ddkbar  = g.deltk ^ dm / ( (2π)^dm )
    magk = unique(g.r)
    Θint = zeros(magk)
    @inbounds for inx in eachindex(magk, Θint)
        Θint[inx] = ddkbar * (-abs2(magk[inx])/2)^n * sum( cXXk[g.r .== magk[inx]] )
    end
    # now loop over ugrid
    vecϕKn  = zeros(ugrid)
    @inbounds @simd for inx in eachindex(ugrid)
        vecϕKn[inx] = dot(Θint, scalebesselj(n, magk.*ugrid[inx]))
    end
    return vecϕKn
end




"""
`YXcov(x1, x2, radKE, radKB)`. Compute the CE, SE, CB, SE cov matrices. Arguments `radKE, radKB`
should be computed by something like the following...
```
cψEk    = mCls.cEEk ./ (g.r.^4)
cψEk[1] = 0.0
vKE2, ugrid = PrimordialB.vecϕKn(2, cψEk, g)
vKE3, ugrid = PrimordialB.vecϕKn(3, cψEk, g)
vKE4, ugrid = PrimordialB.vecϕKn(4, cψEk, g)
vKE5, ugrid = PrimordialB.vecϕKn(5, cψEk, g)
vKE6, ugrid = PrimordialB.vecϕKn(6, cψEk, g)
radKE = hcat(ugrid, vKE2, vKE3, vKE4, vKE5, vKE6)

cψBk    = mCls.cBB0k ./ (g.r.^4)
cψBk[1] = 0.0
vKB2, ugrid = PrimordialB.vecϕKn(2, cψBk, g)
vKB3, ugrid = PrimordialB.vecϕKn(3, cψBk, g)
vKB4, ugrid = PrimordialB.vecϕKn(4, cψBk, g)
vKB5, ugrid = PrimordialB.vecϕKn(5, cψBk, g)
vKB6, ugrid = PrimordialB.vecϕKn(6, cψBk, g)
radKB = hcat(ugrid, vKB2, vKB3, vKB4, vKB5, vKB6)
```
"""
function YXcov(x1, x2, radKE, radKB)
    sqx1   = x1.^2
    sqx2   = x2.^2
    fx1    = x1.^4
    fx2    = x2.^4
    normx  = √(sqx1 + sqx2)
    normx2  = normx.^2
    dnormx2 = sqx1 - sqx2
    dnormx4  = dnormx2.^2

    sKE2xsq = Dierckx.Spline1D(radKE[:,1], radKE[:,2]; k=4, bc="extrapolate", s=0.0)
    sKE3xsq = Dierckx.Spline1D(radKE[:,1], radKE[:,3]; k=4, bc="extrapolate", s=0.0)
    sKE4xsq = Dierckx.Spline1D(radKE[:,1], radKE[:,4]; k=4, bc="extrapolate", s=0.0)
    sKE5xsq = Dierckx.Spline1D(radKE[:,1], radKE[:,5]; k=4, bc="extrapolate", s=0.0)
    sKE6xsq = Dierckx.Spline1D(radKE[:,1], radKE[:,6]; k=4, bc="extrapolate", s=0.0)
    KE2 = reshape(sKE2xsq(vec(normx)), size(normx))
    KE3 = reshape(sKE3xsq(vec(normx)), size(normx))
    KE4 = reshape(sKE4xsq(vec(normx)), size(normx))
    KE5 = reshape(sKE5xsq(vec(normx)), size(normx))
    KE6 = reshape(sKE6xsq(vec(normx)), size(normx))

    sKB2xsq = Dierckx.Spline1D(radKB[:,1], radKB[:,2]; k=4, bc="extrapolate", s=0.0)
    sKB3xsq = Dierckx.Spline1D(radKB[:,1], radKB[:,3]; k=4, bc="extrapolate", s=0.0)
    sKB4xsq = Dierckx.Spline1D(radKB[:,1], radKB[:,4]; k=4, bc="extrapolate", s=0.0)
    sKB5xsq = Dierckx.Spline1D(radKB[:,1], radKB[:,5]; k=4, bc="extrapolate", s=0.0)
    sKB6xsq = Dierckx.Spline1D(radKB[:,1], radKB[:,6]; k=4, bc="extrapolate", s=0.0)
    KB2 = reshape(sKB2xsq(vec(normx)), size(normx))
    KB3 = reshape(sKB3xsq(vec(normx)), size(normx))
    KB4 = reshape(sKB4xsq(vec(normx)), size(normx))
    KB5 = reshape(sKB5xsq(vec(normx)), size(normx))
    KB6 = reshape(sKB6xsq(vec(normx)), size(normx))

    rtnCE, rtnCB         = zeros(normx), zeros(normx)
    rtnSE, rtnSB         = zeros(normx), zeros(normx)
    rtnSCE, rtnSCB       = zeros(normx), zeros(normx)
    rtn∂11CE, rtn∂11CB   = zeros(normx), zeros(normx)
    rtn∂11SE, rtn∂11SB   = zeros(normx), zeros(normx)
    rtn∂11SCE, rtn∂11SCB = zeros(normx), zeros(normx)
    rtn∂22CE, rtn∂22CB   = zeros(normx), zeros(normx)
    rtn∂22SE, rtn∂22SB   = zeros(normx), zeros(normx)
    rtn∂22SCE, rtn∂22SCB = zeros(normx), zeros(normx)
    rtn∂12CE, rtn∂12CB   = zeros(normx), zeros(normx)
    rtn∂12SE, rtn∂12SB   = zeros(normx), zeros(normx)
    rtn∂12SCE, rtn∂12SCB = zeros(normx), zeros(normx)
    @inbounds @simd for inx in eachindex(rtnCE)
        # CE, CB
        t2, t3, t4 = 16, 32*normx2[inx], 16*dnormx4[inx]
        rtnCE[inx] = KE2[inx]*t2 + KE3[inx]*t3 + KE4[inx]*t4
        rtnCB[inx] = KB2[inx]*t2 + KB3[inx]*t3 + KB4[inx]*t4
        # SE, SB
        t2, t3, t4 = 16, 32*normx2[inx], 64*sqx1[inx]*sqx2[inx]
        rtnSE[inx] = KE2[inx]*t2 + KE3[inx]*t3 + KE4[inx]*t4
        rtnSB[inx] = KB2[inx]*t2 + KB3[inx]*t3 + KB4[inx]*t4
        # SCE, SCB
        t4 = 32 * dnormx2[inx] * x1[inx] * x2[inx]
        rtnSCE[inx] =  KE4[inx] * t4
        rtnSCB[inx] = - KB4[inx] * t4

        # ∂11CE, ∂11CB
        t3, t4, t5, t6 = -32*3, -32*18*sqx1[inx], -32*(13fx1[inx]-6sqx1[inx]*sqx2[inx]+fx2[inx]), - 32*(2sqx1[inx]^3-4fx1[inx]*sqx2[inx]+2sqx1[inx]*fx2[inx])
        rtn∂11CE[inx]  = KE3[inx]*t3 + KE4[inx]*t4 + KE5[inx]*t5 + KE6[inx]*t6
        rtn∂11CB[inx]  = KB3[inx]*t3 + KB4[inx]*t4 + KB5[inx]*t5 + KB6[inx]*t6
        # ∂11SE, ∂11SB
        t3, t4, t5, t6 = -32*3, -32*6*(2sqx1[inx]+sqx2[inx]), -32*4*sqx1[inx]*(sqx1[inx]+6sqx2[inx]), -32*8*fx1[inx]*sqx2[inx]
        rtn∂11SE[inx]  = KE3[inx]*t3 + KE4[inx]*t4 + KE5[inx]*t5 + KE6[inx]*t6
        rtn∂11SB[inx]  = KB3[inx]*t3 + KB4[inx]*t4 + KB5[inx]*t5 + KB6[inx]*t6
        # ∂11SCE, ∂11SCB
        t4, t5, t6      = -64*x1[inx]*x2[inx]*3, -64*x1[inx]*x2[inx]*(7sqx1[inx]-3sqx2[inx]), -64*x1[inx]*x2[inx]*2*sqx1[inx]*(sqx1[inx]-sqx2[inx])
        rtn∂11SCE[inx]  = KE4[inx]*t4 + KE5[inx]*t5 + KE6[inx]*t6
        rtn∂11SCB[inx]  = -KB4[inx]*t4 - KB5[inx]*t5 - KB6[inx]*t6


        # ∂22CE, ∂22CB
        t3, t4, t5, t6 = -32*3, -32*18*sqx2[inx], -32*(13fx2[inx]-6sqx2[inx]*sqx1[inx]+fx1[inx]), - 32*(2sqx2[inx]^3-4fx2[inx]*sqx1[inx]+2sqx2[inx]*fx1[inx])
        rtn∂22CE[inx]  = KE3[inx]*t3 + KE4[inx]*t4 + KE5[inx]*t5 + KE6[inx]*t6
        rtn∂22CB[inx]  = KB3[inx]*t3 + KB4[inx]*t4 + KB5[inx]*t5 + KB6[inx]*t6
        # ∂22SE, ∂22SB
        t3, t4, t5, t6 = -32*3, -32*6*(2sqx2[inx]+sqx1[inx]), -32*4*sqx2[inx]*(sqx2[inx]+6sqx1[inx]), -32*8*fx2[inx]*sqx1[inx]
        rtn∂22SE[inx]  = KE3[inx]*t3 + KE4[inx]*t4 + KE5[inx]*t5 + KE6[inx]*t6
        rtn∂22SB[inx]  = KB3[inx]*t3 + KB4[inx]*t4 + KB5[inx]*t5 + KB6[inx]*t6
        # ∂22SCE, ∂22SCB
        t4, t5, t6     = -64*x1[inx]*x2[inx]*(-3), -64*x1[inx]*x2[inx]*(-7sqx2[inx]+3sqx1[inx]), -64*x1[inx]*x2[inx]*2*sqx2[inx]*(sqx1[inx]-sqx2[inx])
        rtn∂22SCE[inx]  = KE4[inx]*t4 + KE5[inx]*t5 + KE6[inx]*t6
        rtn∂22SCB[inx]  = -KB4[inx]*t4 - KB5[inx]*t5 - KB6[inx]*t6


        # ∂12CE, ∂12CB
        t4, t5, t6 = -64*x1[inx]*x2[inx]*3, -64*x1[inx]*x2[inx]*2*normx2[inx], -64*x1[inx]*x2[inx]*dnormx4[inx]
        rtn∂12CE[inx]  = KE4[inx]*t4 + KE5[inx]*t5 + KE6[inx]*t6
        rtn∂12CB[inx]  = KB4[inx]*t4 + KB5[inx]*t5 + KB6[inx]*t6
        # ∂12SE, ∂12SB
        t4, t5, t6 =  -64*x1[inx]*x2[inx]*9, -64*x1[inx]*x2[inx]*6*normx2[inx], -64*x1[inx]*x2[inx]*4*sqx1[inx]*sqx2[inx]
        rtn∂12SE[inx]  = KE4[inx]*t4 + KE5[inx]*t5 + KE6[inx]*t6
        rtn∂12SB[inx]  = KB4[inx]*t4 + KB5[inx]*t5 + KB6[inx]*t6
        # ∂12SCE, ∂12SCB
        t4, t5, t6 = -32*dnormx2[inx]*3, -32*dnormx2[inx]*2*normx2[inx] , -32*dnormx2[inx]*2*2*sqx1[inx]*sqx2[inx]
        rtn∂12SCE[inx]  = KE4[inx]*t4 + KE5[inx]*t5 + KE6[inx]*t6
        rtn∂12SCB[inx]  = -KB4[inx]*t4 - KB5[inx]*t5 - KB6[inx]*t6

    end
    return  rtnCE, rtnCB,
            rtnSE, rtnSB,
            rtnSCE, rtnSCB,
            rtn∂11CE, rtn∂11CB,
            rtn∂11SE, rtn∂11SB,
            rtn∂11SCE, rtn∂11SCB,
            rtn∂22CE, rtn∂22CB,
            rtn∂22SE, rtn∂22SB,
            rtn∂22SCE, rtn∂22SCB,
            rtn∂12CE, rtn∂12CB,
            rtn∂12SE, rtn∂12SB,
            rtn∂12SCE, rtn∂12SCB
end



# ---- variogram computation
"""
Variogram computation for uncertainty in phi. The  `radKϕ` argument should be constructed by
```
vϕK1, ugrid = PrimordialB.vecϕKn(1, varϕk, g)
vϕK2, ugrid = PrimordialB.vecϕKn(2, varϕk, g)
radKϕ = hcat(ugrid, vϕK1, vϕK2)
```
"""
function varioϕx(x1, x2, radKϕ)
    sx1   = abs2(x1)
    sx2   = abs2(x2)
    normx = √(sx1 + sx2)

    sϕK1 = Dierckx.Spline1D(radKϕ[:,1], radKϕ[:,2]; k=4, bc="extrapolate", s=0.0)
    sϕK2 = Dierckx.Spline1D(radKϕ[:,1], radKϕ[:,3]; k=4, bc="extrapolate", s=0.0)
    Kϕ1 = reshape(sϕK1(vec(normx)), size(normx))
    Kϕ2 = reshape(sϕK2(vec(normx)), size(normx))

    rtn∂11, rtn∂12, rtn∂22 = zeros(normx), zeros(normx), zeros(normx)
    @inbounds for inx in eachindex(rtn∂11)
        rtn∂11[inx]  = Kϕ1[inx]*(-2) + Kϕ2[inx]*(-4*sx1[inx])
        rtn∂22[inx]  = Kϕ1[inx]*(-2) + Kϕ2[inx]*(-4*sx2[inx])
        rtn∂12[inx]  = Kϕ2[inx]*(-4*x1[inx]*x2[inx])
    end
    rtn∂11 -= sϕK1(0.0)*(-2)
    rtn∂22 -= sϕK1(0.0)*(-2)
    return -rtn∂11, -rtn∂22, -rtn∂12
end





#######################################################

#  Generates the covariance matrices, with marginalization over nϕ.

#######################################################


function makeΣen_Σb_withoutnoise(xdiff, ydiff, lnxdiff, lnydiff, radKE, radKB, radKϕ)
    ΣCE, ΣCB, ΣSE, ΣSB, ΣSCE, ΣSCB,
    Σ∂11CE, Σ∂11CB, Σ∂11SE, Σ∂11SB, Σ∂11SCE, Σ∂11SCB,
    Σ∂22CE, Σ∂22CB, Σ∂22SE, Σ∂22SB, Σ∂22SCE, Σ∂22SCB,
    Σ∂12CE, Σ∂12CB, Σ∂12SE, Σ∂12SB, Σ∂12SCE, Σ∂12SCB = PrimordialB.YXcov(lnxdiff, lnydiff, radKE, radKB)

    v11Σϕ, v22Σϕ, v12Σϕ = PrimordialB.varioϕx(xdiff, ydiff, radKϕ)

    Σce, Σcb   = Array{Float64}(size(ΣCE)), Array{Float64}(size(ΣCE))
    Σse, Σsb   = Array{Float64}(size(ΣCE)), Array{Float64}(size(ΣCE))
    Σsce, Σscb = Array{Float64}(size(ΣCE)), Array{Float64}(size(ΣCE))
    @inbounds @simd for inx in eachindex(Σce)
        Σce[inx]   = ΣCE[inx]  - v11Σϕ[inx]*Σ∂11CE[inx] - v22Σϕ[inx]*Σ∂22CE[inx] - 2v12Σϕ[inx]*Σ∂12CE[inx]
        Σcb[inx]   = ΣCB[inx]  - v11Σϕ[inx]*Σ∂11CB[inx] - v22Σϕ[inx]*Σ∂22CB[inx] - 2v12Σϕ[inx]*Σ∂12CB[inx]
        Σse[inx]   = ΣSE[inx]  - v11Σϕ[inx]*Σ∂11SE[inx] - v22Σϕ[inx]*Σ∂22SE[inx] - 2v12Σϕ[inx]*Σ∂12SE[inx]
        Σsb[inx]   = ΣSB[inx]  - v11Σϕ[inx]*Σ∂11SB[inx] - v22Σϕ[inx]*Σ∂22SB[inx] - 2v12Σϕ[inx]*Σ∂12SB[inx]
        Σsce[inx]  = ΣSCE[inx] - v11Σϕ[inx]*Σ∂11SCE[inx] - v22Σϕ[inx]*Σ∂22SCE[inx] - 2v12Σϕ[inx]*Σ∂12SCE[inx]
        Σscb[inx]  = ΣSCB[inx] - v11Σϕ[inx]*Σ∂11SCB[inx] - v22Σϕ[inx]*Σ∂22SCB[inx] - 2v12Σϕ[inx]*Σ∂12SCB[inx]
    end
    Σb    = [Σsb Σscb; Σscb Σcb]
    Σen   = [Σce Σsce; Σsce Σse]

    return Σen, Σb
end


function addΣnoise!(Σ, σEErad, deltx, σmod)
    @assert length(σmod) == size(Σ, 1)
    diag_noise = abs2(σmod .* σEErad / deltx)
    @inbounds for k = 1:size(Σ, 1)
        Σ[k,k] += diag_noise[k]
    end
    return Σ
end


function make_v1Σenv2_v1Σbv2(blk1, blk2, vqb1, vqb2, vub1, vub2, x1b1, x1b2, x2b1, x2b2, lnx1b1, lnx1b2, lnx2b1, lnx2b2, radKE, radKB, radKϕ,  σEErad, deltx, σmod)
    x1diff    = x1b1 .- (x1b2).'
    x2diff    = x2b1 .- (x2b2).'
    lnx1diff  = lnx1b1 .- (lnx1b2).'
    lnx2diff  = lnx2b1 .- (lnx2b2).'
    v1  = vcat(vqb1, vub1)
    v2  = vcat(vqb2, vub2)
    Σen, Σb  = makeΣen_Σb_withoutnoise(x1diff, x2diff, lnx1diff, lnx2diff, radKE, radKB, radKϕ)
    if blk1 == blk2
        addΣnoise!(Σen, σEErad, deltx, σmod)
    end
    vΣen = transpose(v1) * Σen * v2
    vΣb  = transpose(v1) * Σb  * v2
    if blk1 != blk2
        vΣen += transpose(v2) * transpose(Σen) * v1
        vΣb  += transpose(v2) * transpose(Σb)  * v1
    end
    return vΣen, vΣb
end




#######################################################

#  likelihood r profile with marginalizing

#######################################################
function loglike_r(rng, qudata, D, V, Db, Vb)
    VLd    = At_mul_B(V, At_mul_B(Vb, qudata) ./ √(Db))
    rout   = Array{Float64}(length(rng))
    for i=1:length(rng)
        @inbounds rout[i] = - 0.5 * dot(VLd, VLd./(D+rng[i])) - 0.5 * sum(log(D+rng[i])) - 0.5 * length(VLd) * log(2π)
    end
    return rout
end



function loglike_r(rng, qudata, Σen, Σb)
    return loglike_r(rng, qudata, D, V, Db, Vb)
end

function get_D_V_Db_Vb(Σen, Σb)
    Db, Vb = eig(Symmetric(Σb))
    Ipos1  = Db .> 0.0
    Db     = Db[Ipos1]
    Vb     = Vb[:,Ipos1]
    B      = transpose(Vb) * Σen * Vb
    B      = transpose(1./√Db) .* B
    B      = B .* (1./√Db)
    D,V    = eig(Symmetric(B))
    # ------
    #Ipos2  = 1e10 .> D .> 0.0
    #Ipos2  = sort(D, rev=true)[10] .> D .> 0.0
    Ipos2  =  D .> 0.0
    # -----
    D      = D[Ipos2]
    V      = V[:,Ipos2]
    return D, V, Db, Vb
end



#######################################################

#  likelihood r profile with marginalizing and cholesky

#######################################################
# function loglike_r_chol(rng, qudata, Σen, Σb)
#     Ub = chol(Σb)
#     UΣenU = (Ub') \ (Σen / Ub)
#     D, V  = eig(Symmetric(UΣenU))
#     VUd    = At_mul_B(V, At_mul_B(Ub, qudata))
#     rout   = Array{Float64}(length(rng))
#     for i=1:length(rng)
#         @inbounds rout[i] = - 0.5 * dot(VUd, VUd./(D+rng[i])) - 0.5 * sum(log(D+rng[i])) - 0.5 * length(VUd) * log(2π)
#     end
#     return rout
# end
function loglike_r_chol(rng, qudata, Σen, Σb)
    Db, Vb = eig((Σen + Σen')./2)
    Ipos1  = Db .> 0.0
    Db     = Db[Ipos1]
    Vb     = Vb[:,Ipos1]
    B      = transpose(Vb) * Σb * Vb
    B      = transpose(1./√Db) .* B
    B      = B .* (1./√Db)
    D,V    = eig((B + B')./2)
    # ------
    #Ipos2  = 1e10 .> D .> 0.0
    #Ipos2  = sort(D, rev=true)[10] .> D .> 0.0
    Ipos2  =  D .> 0.0
    # -----
    D      = D[Ipos2]
    V      = V[:,Ipos2]
    VLd    = At_mul_B(V, At_mul_B(Vb, qudata) ./ √(Db))
    rout   = Array{Float64}(length(rng))
    for i=1:length(rng)
        @inbounds rout[i] = - 0.5 * dot(VLd, VLd./(1.0+rng[i]*D)) - 0.5 * sum(log(1.0+rng[i]*D)) - 0.5 * length(VLd) * log(2π)
    end
    return rout
end



"""
Exact low ell log likelihood for `r` with marginalized delensing residuals.
```
rout = lowℓ_loglike_marg_profile(rng, dbk, lmax, g, mCls, mat_clnbk)
```
"""
function lminlmax_loglike_marg_profile(rng, dbk, lmin, lmax, g, mCls, mat_clnbk)
    hermitianhalf = (g.k[1] .> 0) | ((g.k[1] .== 0.0) & (g.k[2] .> 0.0))
    #lowellbox     = (-lmax .<= g.k[2] .<= lmax) & (-lmax .<= g.k[1] .<= lmax)
    #highellbox    = !( (-lmin .<= g.k[2] .<= lmin) & (-lmin .< g.k[1] .<= lmin) )
    #lrestriction  = lowellbox & highellbox & hermitianhalf
    lrestriction  = (lmin .<= g.r .<= lmax) & hermitianhalf
    cBB0kmax   = mCls.cBB0k[lrestriction]
    cNNkmax    = mCls.cBBnoisek[lrestriction] + mat_clnbk[lrestriction]
    dbkmax     = dbk[lrestriction]
    δ0         = 1/g.deltk^2
    rout       = Array{Float64}(length(rng))
    for i=1:length(rng)
        cddkmax  = (rng[i]*cBB0kmax + cNNkmax)
        rout[i]  = - sum(abs2(dbkmax) ./ cddkmax) / δ0
        rout[i] += - sum(log(cddkmax))
    end
    return rout
end




######################################

#  periodize coordinates (for computing cov matrices) ... currently not used

######################################


function periodize(x, period)
    cst = 2π / (period)
    return angle(exp(im * cst .* x)) ./ cst
end
