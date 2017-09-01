# Non-stationary noise with beam


tic()
addprocs(30) #
#addprocs(7)
println("nproces = $(nprocs())")
@everywhere push!(LOAD_PATH, homedir())
using JLD, ProgressMeter
import PrimordialB
const seedstart = 0xfad38326fbba9e09  # 0x116f205f2f6defb2  # 0x3ef4a4ce1f1432a3  #rand(UInt64)
srand(seedstart)


#####################################################################################
# ------- figure 4 --------
savepath = joinpath("scripts/", "out_figure4.jld") # on Poisson 7/10/2017
μKarcminT          = 1.0
lense_unitr_sc(ek, bk, len, g, order) = PrimordialB.lense_unitr_sc_allorder(ek, bk, len, g, order)


#####################################################################################
###### fixed parameters of the simulation run  ############

nside              = 200
n_sims             = 500
nblks              = 4
sblks_approx       = 1000

pixel_side_length  = 10
phiest_pixel_side_length = 2
period                = pixel_side_length*nside*π/(180*60) # 5, nside*π/(180*60) corresponds to 1 arcmin pixels
phiest_nside          = ceil(Int64, pixel_side_length*nside/phiest_pixel_side_length) #nside * 5
phiest_μKarcminT      = μKarcminT
beamFWHM              = deg2rad(pixel_side_length/60)
phiest_beamFWHM       = deg2rad(phiest_pixel_side_length/60)
g                     = PrimordialB.FFTgrid(PrimordialB.Pix{period, nside}())
fsky_percent          = 100 * (period)^2 / (4π)
fsky_percent_forphi   = 100 * (phiest_pixel_side_length*phiest_nside*π/(180*60))^2 / (4π)
@everywhere rng_true  = [logspace(-3,-2, 20);linspace(0.0125,0.1, 20)] #[1:2:end] # 0.0005:.0025:.1
@everywhere rng       = 0:.0001:.2
order                 = 7 # lensing order
rand_us               = rand(30) # these will determine the location of the noise modulations
pix                   = PrimordialB.Pix{period, nside}()
σmod = let
    ncuts = 30
    σcut  = g.deltx * 2
    magcut = 5.0
    cutindx = rand(1:length(g.r), ncuts)
    kernel(x, y, μx, μy, σcut) = exp(-0.5 * ((x.-μx).^2 + (y.-μy).^2) ./ σcut^2 )
    σmod_tmp = ones(size(g.r))
    for i ∈ 1:ncuts
        σmod_tmp += magcut .* kernel(g.x[1], g.x[2], g.x[1][cutindx[i]], g.x[2][cutindx[i]], σcut)
    end
    σmod_tmp
end


lmin = 2.0
lmax = Inf #0.9*g.nyq
hermitianhalf = (g.k[1] .> 0) | ((g.k[1] .== 0.0) & (g.k[2] .> 0.0))
lrestriction  = (lmin .<= g.r .<= lmax) & hermitianhalf



####################################################################
############# get the cls and simulate phi #######################
cls, mCls, cls_no_B    = PrimordialB.get_spectra(μKarcminT, beamFWHM, period, nside) # for B survey
N0ℓ, Aℓ, mCls_cnϕk, cnϕk_mask = PrimordialB.get_N0ℓAℓ(cls_no_B, phiest_pixel_side_length, phiest_nside, phiest_μKarcminT, phiest_beamFWHM, period, nside)


ϕx, ϕk   = PrimordialB.sim_xk(mCls.cϕϕk, g)
len      = PrimordialB.LenseDecomp(ϕk, zeros(ϕk), g)
nϕx, nϕk = PrimordialB.sim_xk(mCls_cnϕk, g)
ϕestk    = ϕk + nϕk; ϕestk[cnϕk_mask] = 0.0
ρl       = PrimordialB.squash(mCls.cϕϕk ./ √(mCls.cϕϕk) ./ √(mCls_cnϕk + mCls.cϕϕk)); ρl[cnϕk_mask] = 0.0
μϕk      = PrimordialB.squash(mCls.cϕϕk./(mCls_cnϕk + mCls.cϕϕk)) .* ϕestk; μϕk[cnϕk_mask] = 0.0
varϕk    = PrimordialB.squash(mCls_cnϕk .* mCls.cϕϕk ./ (mCls_cnϕk + mCls.cϕϕk)); varϕk[cnϕk_mask] = mCls.cϕϕk[cnϕk_mask]
μlen     = PrimordialB.LenseDecomp(copy(μϕk), zeros(ϕk), g)



####################################################################

let σmod=σmod, varϕk=varϕk, μlen=μlen, mCls=mCls, order=order, pix=pix, ϕestk=ϕestk, ρl=ρl, lrestriction=lrestriction, μKarcminT=μKarcminT
    @everywhere const gl = PrimordialB.FFTgrid($pix)
    @everywhere const φ2_l     = 2angle.(gl.k[1] + im * gl.k[2])
    @everywhere const cos²φ    = abs2.(cos.(φ2_l))
    @everywhere const sin²φ    = abs2.(sin.(φ2_l))
    @everywhere const sincosφ  = cos.(φ2_l).*sin.(φ2_l)
    @everywhere const beamCovk =  $(mCls.cBBnoisek) ./ (2π)
    @everywhere const covϕk     = $varϕk ./ (2π)
    @everywhere const covCEk    = cos²φ .* $(mCls.cEEk) ./ (2π)
    @everywhere const covSEk    = sin²φ .* $(mCls.cEEk) ./ (2π)
    @everywhere const covSCEk   = sincosφ .* $(mCls.cEEk) ./ (2π)
    @everywhere const covCB0k   = cos²φ .* $(mCls.cBB0k) ./ (2π)
    @everywhere const covSB0k   = sin²φ .* $(mCls.cBB0k) ./ (2π)
    @everywhere const covSCB0k  = - sincosφ .* $(mCls.cBB0k) ./ (2π)
    @everywhere const local_order = $order
    @everywhere const local_μlen = $μlen
    @everywhere const local_lrestriction = $lrestriction
    @everywhere const local_ϕestk = $ϕestk
    @everywhere const local_ρl = $ρl
    @everywhere const local_mCls = $mCls
    @everywhere const local_pix = $pix
    @everywhere const local_μKarcminT = $μKarcminT
    @everywhere const local_σmod = $σmod

    @everywhere function cov_col(pix_index)
        exsy = exp.(- im .* gl.k[1] .* (gl.x[1][pix_index] + local_μlen.displx[pix_index]) - im .* gl.k[2] .* (gl.x[2][pix_index] + local_μlen.disply[pix_index]))::Array{Complex{Float64},2}
        exsy_nolen = exp.(- im .* gl.k[1] .* gl.x[1][pix_index] - im .* gl.k[2] .* gl.x[2][pix_index])::Array{Complex{Float64},2}

        ΣCE, ΣCB          = PrimordialB.lense(real.(gl.FFT \ (exsy.*covCEk)),  real.(gl.FFT \ (exsy.*covCB0k)),  local_μlen, gl, local_order)
        ΣSE, ΣSB          = PrimordialB.lense(real.(gl.FFT \ (exsy.*covSEk)),  real.(gl.FFT \ (exsy.*covSB0k)),  local_μlen, gl, local_order)
        ΣSCE, ΣSCB        = PrimordialB.lense(real.(gl.FFT \ (exsy.*covSCEk)), real.(gl.FFT \ (exsy.*covSCB0k)), local_μlen, gl, local_order)
        Σ∂11CE,  Σ∂11CB   = PrimordialB.lense(real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[1].*covCEk)),  real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[1].*covCB0k)),  local_μlen, gl, local_order)
        Σ∂11SE,  Σ∂11SB   = PrimordialB.lense(real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[1].*covSEk)),  real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[1].*covSB0k)),  local_μlen, gl, local_order)
        Σ∂11SCE, Σ∂11SCB  = PrimordialB.lense(real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[1].*covSCEk)), real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[1].*covSCB0k)), local_μlen, gl, local_order)
        Σ∂12CE,  Σ∂12CB   = PrimordialB.lense(real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[2].*covCEk)),  real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[2].*covCB0k)),  local_μlen, gl, local_order)
        Σ∂12SE,  Σ∂12SB   = PrimordialB.lense(real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[2].*covSEk)),  real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[2].*covSB0k)),  local_μlen, gl, local_order)
        Σ∂12SCE, Σ∂12SCB  = PrimordialB.lense(real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[2].*covSCEk)), real.(gl.FFT \ (exsy.*gl.k[1].*gl.k[2].*covSCB0k)), local_μlen, gl, local_order)
        Σ∂22CE,  Σ∂22CB   = PrimordialB.lense(real.(gl.FFT \ (exsy.*gl.k[2].*gl.k[2].*covCEk)),  real.(gl.FFT \ (exsy.*gl.k[2].*gl.k[2].*covCB0k)),  local_μlen, gl, local_order)
        Σ∂22SE,  Σ∂22SB   = PrimordialB.lense(real.(gl.FFT \ (exsy.*gl.k[2].*gl.k[2].*covSEk)),  real.(gl.FFT \ (exsy.*gl.k[2].*gl.k[2].*covSB0k)),  local_μlen, gl, local_order)
        Σ∂22SCE, Σ∂22SCB  = PrimordialB.lense(real.(gl.FFT \ (exsy.*gl.k[2].*gl.k[2].*covSCEk)), real.(gl.FFT \ (exsy.*gl.k[2].*gl.k[2].*covSCB0k)), local_μlen, gl, local_order)

        ∂12Σϕ   = real.(gl.FFT \ (exsy_nolen.*gl.k[1].*gl.k[2].*covϕk))
        ∂11Σϕ   = real.(gl.FFT \ (exsy_nolen.*gl.k[1].*gl.k[1].*covϕk))
        ∂22Σϕ   = real.(gl.FFT \ (exsy_nolen.*gl.k[2].*gl.k[2].*covϕk))

        v11Σϕ = ∂11Σϕ[pix_index] .- ∂11Σϕ
        v22Σϕ = ∂22Σϕ[pix_index] .- ∂22Σϕ
        v12Σϕ = ∂12Σϕ[pix_index] .- ∂12Σϕ

        Σce   = ΣCE  - v11Σϕ .* Σ∂11CE - v22Σϕ .* Σ∂22CE - 2v12Σϕ .* Σ∂12CE
        Σcb   = ΣCB  - v11Σϕ .* Σ∂11CB - v22Σϕ .* Σ∂22CB - 2v12Σϕ .* Σ∂12CB
        Σse   = ΣSE  - v11Σϕ .* Σ∂11SE - v22Σϕ .* Σ∂22SE - 2v12Σϕ .* Σ∂12SE
        Σsb   = ΣSB  - v11Σϕ .* Σ∂11SB - v22Σϕ .* Σ∂22SB - 2v12Σϕ .* Σ∂12SB
        Σsce  = ΣSCE - v11Σϕ .* Σ∂11SCE - v22Σϕ .* Σ∂22SCE - 2v12Σϕ .* Σ∂12SCE
        Σscb  = ΣSCB - v11Σϕ .* Σ∂11SCB - v22Σϕ .* Σ∂22SCB - 2v12Σϕ .* Σ∂12SCB

        ## stationary noise case....
        #Σbeam =  real.(gl.FFT \ (exsy_nolen.*beamCovk))
        #Σce += Σbeam
        #Σse += Σbeam

        # non-stationary noise case
        σ2xφxmy   = abs2.(local_σmod) .* real.(gl.FFT \ (exsy_nolen .* sqrt.(beamCovk) ))
        Σbeam = real.(gl.FFT \ (sqrt.(beamCovk) .* (gl.FFT * σ2xφxmy)))
        Σce += Σbeam
        Σse += Σbeam


        return Σce, Σcb, Σse, Σsb, Σsce, Σscb # vectorize these to get the column
    end
    @everywhere function VQcolQB(pix_index)
        qx = zeros(gl.nside, gl.nside)
        ux = zeros(gl.nside, gl.nside)
        qx[pix_index] = 1.0
        qk = gl.FFT * qx
        uk = gl.FFT * ux
        ek, bk,   = PrimordialB.qu2eb(qk, uk, gl)
        δ1residk  = bk - PrimordialB.estδ1BfromEk(ek,local_ϕestk, local_ρl, local_mCls, gl)
        return δ1residk[local_lrestriction]
    end
    @everywhere function VUcolQB(pix_index)
        qx = zeros(gl.nside, gl.nside)
        ux = zeros(gl.nside, gl.nside)
        ux[pix_index] = 1.0
        qk = gl.FFT * qx
        uk = gl.FFT * ux
        ek, bk,   = PrimordialB.qu2eb(qk, uk, gl)
        δ1residk  = bk - PrimordialB.estδ1BfromEk(ek, local_ϕestk, local_ρl, local_mCls, gl)
        return δ1residk[local_lrestriction]
    end

    @everywhere function VQcolE(pix_index)
        qx = zeros(gl.nside, gl.nside)
        ux = zeros(gl.nside, gl.nside)
        qx[pix_index] = 1.0
        qk = gl.FFT * qx
        uk = gl.FFT * ux
        ek, bk,   = PrimordialB.qu2eb(qk, uk, gl)
        return ek[local_lrestriction]
    end
    @everywhere function VUcolE(pix_index)
        qx = zeros(gl.nside, gl.nside)
        ux = zeros(gl.nside, gl.nside)
        ux[pix_index] = 1.0
        qk = gl.FFT * qx
        uk = gl.FFT * ux
        ek, bk,   = PrimordialB.qu2eb(qk, uk, gl)
        return ek[local_lrestriction]
    end
    @everywhere function VQcolB(pix_index)
        qx = zeros(gl.nside, gl.nside)
        ux = zeros(gl.nside, gl.nside)
        qx[pix_index] = 1.0
        qk = gl.FFT * qx
        uk = gl.FFT * ux
        ek, bk,   = PrimordialB.qu2eb(qk, uk, gl)
        return bk[local_lrestriction]
    end
    @everywhere function VUcolB(pix_index)
        qx = zeros(gl.nside, gl.nside)
        ux = zeros(gl.nside, gl.nside)
        ux[pix_index] = 1.0
        qk = gl.FFT * qx
        uk = gl.FFT * ux
        ek, bk,   = PrimordialB.qu2eb(qk, uk, gl)
        return bk[local_lrestriction]
    end
    @everywhere function block_VqVu(vqcol::Function, vucol::Function, indrng, nside)
        nb  = length(indrng)
        npx = nside^2
        vq = zeros(Complex{Float64}, nb, npx)
        vu = zeros(Complex{Float64}, nb, npx)
        ProgressMeter.@showprogress 4 "constructing Vq, Vu..." for i in 1:npx
            vq[:,i] = vqcol(i)[indrng]
            vu[:,i] = vucol(i)[indrng]
        end
        Vq = vcat(real.(vq), imag.(vq))
        Vu = vcat(real.(vu), imag.(vu))
        return Vq, Vu
    end
end


########## get projection matrices for local blocks ############
sblks = floor(Int64, min(sum(lrestriction)/nblks, sblks_approx))
# ------------ order the modes into signal to noise ---------
mat_clnbk  = PrimordialB.ResB_matpwr(mCls, ρl, g)
signal_to_noiseQB = mCls.cBB0k ./ (mCls.cBBnoisek + mat_clnbk)
signal_to_noiseE  = mCls.cEEk  ./ mCls.cEEnoisek
VqVu_SN_weightsQB = signal_to_noiseQB[lrestriction]
VqVu_SN_weightsE  = signal_to_noiseE[lrestriction]
permQB = sortperm(VqVu_SN_weightsQB, rev=true) # these are the modes that the quadratic estimate produces
permE  = sortperm(VqVu_SN_weightsE, rev=true) # these are the modes that the quadratic estimate produces
indxQB = (1:length(VqVu_SN_weightsQB))[permQB]
indxE  = (1:length(VqVu_SN_weightsE))[permE]
# ------------ find the projection matrices ------------
VqVu = @sync @parallel (vcat) for i in 1:nblks
    VqQB, VuQB = block_VqVu(VQcolQB,  VUcolQB, indxQB[(i-1)*sblks+(1:sblks)], nside)
    VqE, VuE   = block_VqVu(VQcolE,    VUcolE, indxE[(i-1)*sblks+(1:sblks)], nside)
    Vq = vcat(VqQB, VqE)
    Vu = vcat(VuQB, VuE)
    (i, Vq, Vu)
end
Vq = Dict(VqVu[i][1] => VqVu[i][2] for i in 1:nblks)
Vu = Dict(VqVu[i][1] => VqVu[i][3] for i in 1:nblks)
#clear memory
VqVu = 0.0
@everywhere gc()

# #TODO: try a Vecchia approach and have these blocks overlap.
for i = 1:nblks
    v = hcat(Vq[i], Vu[i]) # ------- now reduce vq, vu by finding an ortho basis of the range
    newv = transpose(PrimordialB.orth(transpose(v); cut = 1e-15))
    Vq[i] = newv[:,1:size(Vq[i],2)]
    Vu[i] = newv[:,size(Vu[i],2)+1:end]
    @everywhere gc()
end


# ------------ compute the covariance of the projections ------------
npx = nside^2
max_blk_size = 2*nprocs()
partition(x, max_blk_size) = [x[i:min(i+max_blk_size-1,length(x))] for i in 1:max_blk_size:length(x)]
pixl_partition = partition(1:npx, max_blk_size)
vΣev = Dict(i=>zeros(size(Vq[i],1), size(Vq[i],1)) for i=1:nblks)
vΣbv = Dict(i=>zeros(size(Vq[i],1), size(Vq[i],1)) for i=1:nblks)
@showprogress 4 "constructing vΣev, vΣbv..." for pixel_rnge in pixl_partition
    vecΣeb = @sync @parallel (hcat) for i in pixel_rnge
        Σce, Σcb, Σse, Σsb, Σsce, Σscb = cov_col(i)
        vcat(vec(Σce), vec(Σcb), vec(Σse), vec(Σsb), vec(Σsce), vec(Σscb))
    end
    tΣce, tΣcb, tΣse, tΣsb, tΣsce, tΣscb = vecΣeb[1:npx,:], vecΣeb[npx+1:2npx,:], vecΣeb[2npx+1:3npx,:], vecΣeb[3npx+1:4npx,:], vecΣeb[4npx+1:5npx,:], vecΣeb[5npx+1:6npx,:]
    for i in 1:nblks
        vΣev[i]  += (Vq[i] * tΣce  + Vu[i] * tΣsce) * transpose(Vq[i][:,pixel_rnge])
        vΣev[i]  += (Vq[i] * tΣsce + Vu[i] * tΣse)  * transpose(Vu[i][:,pixel_rnge])
        vΣbv[i]  += (Vq[i] * tΣsb  + Vu[i] * tΣscb) * transpose(Vq[i][:,pixel_rnge])
        vΣbv[i]  += (Vq[i] * tΣscb + Vu[i] * tΣcb)  * transpose(Vu[i][:,pixel_rnge])
    end
    @everywhere gc()
end

# ------------ get the eigen decomp of covariance of the projections ------------
@time vec_DVDbVb = @sync @parallel (vcat) for i = 1:nblks
    local_vΣev = remotecall_fetch(()->vΣev[i], 1)
    local_vΣbv = remotecall_fetch(()->vΣbv[i], 1)
    vDv, vVv, vDbv, vVbv = PrimordialB.get_D_V_Db_Vb(local_vΣev, local_vΣbv)
    i => (vDv, vVv, vDbv, vVbv)
end
dic_DVDbVb =  Dict(vec_DVDbVb)
# clear memory
vec_DVDbVb = 0
@everywhere gc()



# since the simulations below are done on the master node, we can release the workers
rmprocs(workers())


# ----------- simulate
rest_B      = zeros(length(rng_true), n_sims)
rest_modes  = zeros(length(rng_true), n_sims)
@showprogress 4 "simulations..." for ctr = 1:n_sims
    # simulate data of the form the form
    # dqx   = dncex + √(r) .* dsbx # <--- the data
    # dux   = dnsex + √(r) .* dcbx # <--- the data

    # # for stationary noise
    # nex, nek                   = PrimordialB.sim_xk(mCls.cEEnoisek, g)
    # nbx, nbk                   = PrimordialB.sim_xk(mCls.cBBnoisek, g)
    # nqk, nuk, nqx, nux         = PrimordialB.eb2qu(nek, nbk, g)

    # for non-stationary noise
    nex, nek                   = PrimordialB.sim_xk(1. + 0*mCls.cEEnoisek, g)
    nbx, nbk                   = PrimordialB.sim_xk(1. + 0*mCls.cBBnoisek, g)
    nqk, nuk, nqx, nux         = PrimordialB.eb2qu(nek, nbk, g)
    nqx .*= σmod
    nux .*= σmod
    nqx = real.(g.FFT \ (sqrt.(mCls.cEEnoisek) .* (g.FFT * nqx)))
    nux = real.(g.FFT \ (sqrt.(mCls.cEEnoisek) .* (g.FFT * nux)))

    ex, ek                     = PrimordialB.sim_xk(mCls.cEEk, g)
    #bx, bk                     = sim_xk(mCls.cBBk, g)
    bx, bk                     = PrimordialB.sim_xk(mCls.cBB0k, g)
    ln_cex, ln_sex, ln_cbx, ln_sbx  = lense_unitr_sc(ek, bk, len, g, order)
    dncex = ln_cex + nqx
    dnsex = ln_sex + nux
    dcbx  = ln_cbx
    dsbx  = ln_sbx

    # ---- blakes method all ell
    rll_B  = PrimordialB.full_quadratic_B_delenser_qxux(lmin, lmax, dncex, dnsex, dsbx, dcbx, mCls, rng, rng_true, ρl, ϕestk, period, nside)
    rll_B_maxind = mapslices(x->findmax(x)[2], rll_B, 2)[:]
    rest_B[:,ctr] = rng[rll_B_maxind]

    rll_modes = sum(1:nblks) do i
        vDv, vVv, vDbv, vVbv = dic_DVDbVb[i]
        qudata_en = Vq[i] * dncex[:] + Vu[i] * dnsex[:]
        qudata_b  = Vq[i] * dsbx[:]  + Vu[i] * dcbx[:]
        PrimordialB.loglike_r_rtrue(rng, rng_true, qudata_en, qudata_b, vDv, vVv, vDbv, vVbv)
    end
    rll_modes_maxind = mapslices(x->findmax(x)[2], rll_modes, 2)[:]
    rest_modes[:,ctr] = rng[rll_modes_maxind]
end


runtime = toc()

@save "$savepath" rest_B rest_modes rng rng_true n_sims len μKarcminT beamFWHM nside period seedstart N0ℓ Aℓ varϕk σmod pixel_side_length

message = """
script1.jl job completed in $(round(Int64, runtime/60)) minutes.
fsky is at $(round(fsky_percent,1)) %
μKarcminT = $(μKarcminT)
"""
run(pipeline(`echo $message`, `mail 7733153613@mms.att.net`))
println(message)
