module PrimordialB

using  Dierckx, JLD, PyCall #, SpecialFunctions
export
	# Custom types
	FFTgrid,
	MatrixCls,
	LenseDecomp,
    Pix,

	# methods
	lense,
	sim_xk,
	radial_power,
	qu2eb, eb2qu,
	squash, squash!


# FFTW.set_num_threads(CPU_CORES)
# FFTW.set_num_threads(3)

#=  To lint this file run:
using Lint
lintfile("src/PrimordialB.jl")
=#





#=##########################################################

Custom type Definitions

=##############################################################

# ---- Holds grid, model and planned FFT parameters for the quadratic estimate.
immutable FFTgrid{dm, T}
	period::Float64
	nside::Int64
	deltx::Float64
	deltk::Float64
	nyq::Float64
	x::Array{Array{Float64,dm},1}
	k::Array{Array{Float64,dm},1}
	r::Array{Float64,dm}
	FFT::T  # saved plan for fast fft
end


#---- Holds the cls expanded out to the 2 d spectral matrices.
immutable MatrixCls{dm}
	cϕϕk::Array{Float64,dm}
	cϕψk::Array{Float64,dm}
	cψψk::Array{Float64,dm}
	cTTk::Array{Float64,dm}
	cTEk::Array{Float64,dm}
	cEEk::Array{Float64,dm}
	cBBk::Array{Float64,dm}
	cBB0k::Array{Float64,dm}
	cTTnoisek::Array{Float64,dm}
	cEEnoisek::Array{Float64,dm}
	cBBnoisek::Array{Float64,dm}
end


# ----- Holds the decomposition of the lensing displacements
immutable LenseDecomp
	indcol::Array{Int64,2}
	indrow::Array{Int64,2}
	rdisplx::Array{Float64,2}
	rdisply::Array{Float64,2}
	displx::Array{Float64,2}
	disply::Array{Float64,2}
	ϕk::Array{Complex{Float64},2}
	ψk::Array{Complex{Float64},2}
end



##########################################################

# Type constructors

##############################################################

# TODO: fix this up for general dimension
immutable Pix{period, nside} end

@generated function FFTgrid{period, nside}(p::Pix{period, nside})
    deltx     = period / nside
	deltk     = 2π / period
    nyq       = 2π / (2deltx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * deltk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * deltx
    zrs    = zeros(nside, nside)
    k      = [k_side.' .+ zrs, k_side .+ zrs]#[reshape(k_side, 1, nside), reshape(k_side, nside, 1)]
    x      = [x_side.' .+ zrs, x_side .+ zrs]#[reshape(x_side, 1, nside), reshape(x_side, nside, 1)]
    r      = sqrt.(k[1].^2 + k[2].^2)
    FFT    =  ((deltx^2) / (2π)) * plan_fft(rand(nside,nside); flags=FFTW.ESTIMATE, timelimit=5)
    return :(FFTgrid{2, typeof($(FFT))}($period, $nside, $deltx, $deltk, $nyq, $x, $k, $r, $(FFT)))
end




function MatrixCls{dm,T}(g::FFTgrid{dm,T}, cls; σTTrad=0.0, σEErad=0.0,  σBBrad=0.0, beamFWHM=0.0)
	cϕϕk = cls_to_cXXk(cls[:ell], cls[:ϕϕ], g.r)
	cϕψk = cls_to_cXXk(cls[:ell], cls[:ϕψ], g.r)
	cψψk = cls_to_cXXk(cls[:ell], cls[:ψψ], g.r)
	cTTk = cls_to_cXXk(cls[:ell], cls[:tt], g.r)
	cTEk = cls_to_cXXk(cls[:ell], cls[:te], g.r)
	cEEk = cls_to_cXXk(cls[:ell], cls[:ee], g.r)
	cBBk = cls_to_cXXk(cls[:ell], cls[:bb], g.r)
	cBB0k= cls_to_cXXk(cls[:ell], cls[:bb0], g.r)
	cTTnoisek = cNNkgen(g.r; σrad=σTTrad, beamFWHM=beamFWHM)
	cEEnoisek = cNNkgen(g.r; σrad=σEErad, beamFWHM=beamFWHM)
	cBBnoisek = cNNkgen(g.r; σrad=σBBrad, beamFWHM=beamFWHM)
	MatrixCls{dm}(cϕϕk, cϕψk, cψψk, cTTk, cTEk, cEEk, cBBk, cBB0k, cTTnoisek, cEEnoisek, cBBnoisek)
end


function LenseDecomp(ϕk, ψk, g)
	#displx = real(g.FFT \ (im .* g.k[1] .* ϕk) +  g.FFT \ (im .* g.k[2] .* ψk))
	#disply = real(g.FFT \ (im .* g.k[2] .* ϕk) -  g.FFT \ (im .* g.k[1] .* ψk))
	displx, disply = LenseDecomp_helper1(ϕk, ψk, g)

	row, col  = size(g.x[1])
	indcol    = Array(Int64, row, col)
	indrow	  = Array(Int64, row, col)
	rdisplx   = Array(Float64, row, col)
	rdisply   = Array(Float64, row, col)
	@inbounds for j = 1:col, i = 1:row
		round_displx_deltx = round(Int64, displx[i,j]/g.deltx)
		round_disply_deltx = round(Int64, disply[i,j]/g.deltx)
	    indcol[i,j]  = indexwrap(j + round_displx_deltx, col)
	    indrow[i,j]  = indexwrap(i + round_disply_deltx, row)
		rdisplx[i,j] = displx[i,j] - g.deltx * round_displx_deltx
	    rdisply[i,j] = disply[i,j] - g.deltx * round_disply_deltx
	end
	return LenseDecomp(indcol, indrow, rdisplx, rdisply, displx, disply, copy(ϕk), copy(ψk))
end
function LenseDecomp_helper1(ϕk, ψk, g)
	tmpdxk = Array(Complex{Float64}, size(g.r))
	tmpdyk = Array(Complex{Float64}, size(g.r))
	@inbounds @simd for i in eachindex(tmpdxk, tmpdyk)
		tmpdxk[i] = complex(im * g.k[1][i] * ϕk[i] + im * g.k[2][i] * ψk[i])
		tmpdyk[i] = complex(im * g.k[2][i] * ϕk[i] - im * g.k[1][i] * ψk[i])
	end
	displx = real(g.FFT \ tmpdxk)
	disply = real(g.FFT \ tmpdyk)
	return displx, disply
end

LenseDecomp(len::LenseDecomp, g) = LenseDecomp(copy(len.ϕk), copy(len.ψk), g)





##########################################################

# Helper functions for the type constructors

##############################################################

indexwrap(ind::Int64, uplim)  = mod(ind - 1, uplim) + 1


function cls_to_cXXk{dm}(ell, cxxls, r::Array{Float64, dm})
		spl = Dierckx.Spline1D(ell, cxxls; k=1, bc="zero", s=0.0)
		rtn = squash(map(spl, r))::Array{Float64, dm}
		rtn[r.==0.0] = 0.0
		return squash(map(spl, r))::Array{Float64, dm}
end


function cNNkgen{dm}(r::Array{Float64,dm}; σrad=0.0, beamFWHM=0.0)
	beamSQ = exp(- (beamFWHM ^ 2) * abs2(r) ./ (8 * log(2)) )    # my misunderstanding ? anyway we alway use zero beam here
    rtn = ones(size(r)) .* abs2(σrad) ./ beamSQ
    rtn[r.==0.0] = 0.0
    return rtn
end


function getgrid{T}(g::FFTgrid{2,T})
	xco_side, kco_side = getxkside(g)
	kco1, kco2 = meshgrid(kco_side, kco_side)
	xco1, xco2 = meshgrid(xco_side, xco_side)
	kco    = Array{Float64,2}[kco1, kco2]
	xco    = Array{Float64,2}[xco1, xco2]
	return xco, kco
end


function getxkside{dm,T}(g::FFTgrid{dm,T})
	deltx    = g.period / g.nside
	deltk    = 2π / g.period
	xco_side = zeros(g.nside)
	kco_side = zeros(g.nside)
	for j in 0:(g.nside-1)
		xco_side[j+1] = (j < g.nside/2) ? (j*deltx) : (j*deltx - g.period)
		kco_side[j+1] = (j < g.nside/2) ? (j*deltk) : (j*deltk - 2*π*g.nside/g.period)
	end
	xco_side, kco_side
end


function meshgrid(side_x,side_y)
    	nx = length(side_x)
    	ny = length(side_y)
    	xt = repmat(vec(side_x).', ny, 1)
    	yt = repmat(vec(side_y)  , 1 , nx)
    	return xt, yt
end




###############################

# local likelihood functions

###############################
include("likewhite.jl")




###############################

# functions for simulation code

###############################
include("simulation_src.jl")



###############################

# the quadratic estimate

###############################


""" EB quadratic estimate """
function eb_quad_est(dek, dbk, cEEfid, cBBfid, cEEobsfid, cBBobsfid, g, Aℓ)
    CX = 1./cEEobsfid
    CY = 1./cBBobsfid
    estk  = -im .* g.k[1] .* sin_conv(im .* g.k[1] .* cEEfid .* dek, dbk, CX, CY, g)
    estk += -im .* g.k[2] .* sin_conv(im .* g.k[2] .* cEEfid .* dek, dbk, CX, CY, g)
    estk += -im .* g.k[1] .* sin_conv(dek,  im .* g.k[1] .* cBBfid .* dbk, CX, CY, g)
    estk += -im .* g.k[2] .* sin_conv(dek,  im .* g.k[2] .* cBBfid .* dbk, CX, CY, g)
    rtnk = Aℓ .* estk
    rtnx = real(g.FFT \ rtnk)
    return rtnk, rtnx
end

""" Aℓ for EB quadratic estimate"""
function getAℓ(cEEtru, cBBtru, cEEfid, cBBfid, cEEobsfid, cBBobsfid, g)
    tmp = zeros(Complex{Float64}, size(cEEtru))
    ons = ones(size(tmp))
    CX = 1./cEEobsfid
    CY = 1./cBBobsfid
    for p = 1:2, q = 1:2
        tmp += -g.k[p] .* g.k[q] .* sin2_conv( g.k[p] .* cEEfid, g.k[q] .* cBBtru, CX, CY, g)
        tmp += -g.k[p] .* g.k[q] .* sin2_conv( g.k[p] .* cEEtru, g.k[q] .* cBBfid, CX, CY, g)
        tmp +=  g.k[p] .* g.k[q] .* sin2_conv( g.k[p] .* g.k[q] .* cEEfid .* cEEtru, ons, CX, CY, g)
        tmp +=  g.k[p] .* g.k[q] .* sin2_conv( ons, g.k[p] .* g.k[q] .* cBBfid .* cBBtru, CX, CY, g)
    end
    tmp .*= 1 / 2 / π
    rtnk  = 1 ./ real(tmp)
    return squash(rtnk)
end

""" Aℓ for EB quadratic estimate"""
function getN0ℓ(cEEfid, cBBfid, cEEobsfid, cBBobsfid, cEEobstru, cBBobstru, cEBobstru, g, Aℓ)
    tmp = zeros(Complex{Float64}, size(cEEfid))
    ons = ones(size(tmp))

    CX = cEEobstru./abs2(cEEobsfid)
    CY = cBBobstru./abs2(cBBobsfid)
    for p = 1:2, q = 1:2
        tmp += -g.k[p] .* g.k[q] .* sin2_conv( g.k[p] .* cEEfid, g.k[q] .* cBBfid, CX, CY, g)
        tmp += -g.k[p] .* g.k[q] .* sin2_conv( g.k[p] .* cEEfid, g.k[q] .* cBBfid, CX, CY, g)
        tmp +=  g.k[p] .* g.k[q] .* sin2_conv( g.k[p] .* g.k[q] .* cEEfid .* cEEfid, ons, CX, CY, g)
        tmp +=  g.k[p] .* g.k[q] .* sin2_conv( ons, g.k[p] .* g.k[q] .* cBBfid .* cBBfid, CX, CY, g)
    end
    CX = cEBobstru ./ cEEobsfid./ cBBobsfid
    CY = cEBobstru ./ cEEobsfid./ cBBobsfid
    for p = 1:2, q = 1:2
        tmp +=  g.k[p] .* g.k[q] .* sin2_conv( g.k[p] .* cEEfid, g.k[q] .* cEEfid, CX, CY, g)
        tmp +=  g.k[p] .* g.k[q] .* sin2_conv( g.k[p] .* cBBfid, g.k[q] .* cBBfid, CX, CY, g)
        tmp += -g.k[p] .* g.k[q] .* sin2_conv( ons, g.k[p] .* g.k[q] .* cEEfid .* cBBfid, CX, CY, g)
        tmp += -g.k[p] .* g.k[q] .* sin2_conv( g.k[p] .* g.k[q] .* cEEfid .* cBBfid,  ons, CX, CY, g)
    end
    tmp .*= 1 / 2 / π
    rtnk  = abs2(Aℓ) .* real(tmp)
    return squash(rtnk)
end

function clenEBk(mCls, g)
    φ_k  = angle(g.k[1] + im * g.k[2])
    cQQk = cos(2φ_k).^2 .* mCls.cEEk + sin(2φ_k).^2 .* mCls.cBBk
    cUUk = sin(2φ_k).^2 .* mCls.cEEk + cos(2φ_k).^2 .* mCls.cBBk
    cQUk = cos(2φ_k) .* sin(2φ_k) .* (mCls.cEEk - mCls.cBBk)
    δ1Qδ1Qk = clenEBk_helper(cQQk, mCls.cϕϕk, g)
    δ1Uδ1Uk = clenEBk_helper(cUUk, mCls.cϕϕk, g)
    δ1Qδ1Uk = clenEBk_helper(cQUk, mCls.cϕϕk, g)
    rtnk   =  cos(2φ_k) .* sin(2φ_k) .* (δ1Uδ1Uk - δ1Qδ1Qk)
    rtnk  +=  (cos(2φ_k).^2 - sin(2φ_k).^2) .* δ1Qδ1Uk
    return real(rtnk)
end
# (cXYk, cϕϕk, g) → ∫dk̲̲ kp * kq * (k+l)p * (k+l)q * cXYk * cϕϕk
function clenEBk_helper(cXYk, cϕϕk, g)
    rtnk = zeros(Complex{Float64}, size(cϕϕk))
    for p = 1:2, q = 1:2
        rtnk += conv( g.k[p].*g.k[q].*cXYk,  g.k[p].*g.k[q].*cϕϕk, g)
    end
    rtnk .*= 1 / 2 / π
    return rtnk
end


function δ01E_δ01B(cEEk, cBBk, ϕk, g)
    φ_k     = angle(g.k[1] + im * g.k[2])
    ex, ek  = sim_xk(cEEk, g)
    bx, bk  = sim_xk(cBBk, g)
    qk, uk, qx, ux  = PrimordialB.eb2qu(ek, bk, g)
    δ0Ek = copy(ek)
    δ0Bk = copy(bk)
    δ1Ek = zeros(Complex{Float64}, size(g.x[1]))
    δ1Bk = zeros(Complex{Float64}, size(g.x[1]))
    for p = 1:2
        δ1Ek +=  - cos(2φ_k).*conv(im.*g.k[p].*qk, im.*g.k[p].*ϕk, g)  - sin(2φ_k).*conv(im.*g.k[p].*uk, im.*g.k[p].*ϕk, g)
        δ1Bk +=    sin(2φ_k).*conv(im.*g.k[p].*qk, im.*g.k[p].*ϕk, g)  - cos(2φ_k).*conv(im.*g.k[p].*uk, im.*g.k[p].*ϕk, g)
    end
    return δ0Ek, δ1Ek, δ0Bk, δ1Bk
end


function simulateEB_N1(itrs, cEEtru, cBBtru, cEEfid, cBBfid, cEEobsfid, cBBobsfid, mCls, g, Aℓ)
    N1k = zeros(Float64, size(g.x[1]))
    for cntr = 1:itrs
        ϕx, ϕk  = sim_xk(mCls.cϕϕk, g)
        δ0Ek_3, δ1Ek_3, δ0Bk_3, δ1Bk_3 = δ01E_δ01B(cEEtru, cBBtru, ϕk, g)
        δ0Ek_4, δ1Ek_4, δ0Bk_4, δ1Bk_4 = δ01E_δ01B(cEEtru, cBBtru, ϕk, g)
        ϕestk_3, = eb_quad_est(δ0Ek_4, δ1Bk_3, cEEfid, cBBfid, cEEobsfid, cBBobsfid, g, Aℓ)
        ϕestk_4, = eb_quad_est(δ0Ek_3, δ1Bk_4, cEEfid, cBBfid, cEEobsfid, cBBobsfid, g, Aℓ)
        N1k += ϕestk_3 .* conj(ϕestk_4) .* (g.deltk^2)
    end
    return real(N1k./ itrs)
end



""" (Xk, Yk, g) → ∫dk̲ X_{k+l} * Y_{-k}"""
function conv(Xk, Yk, g)
    Xx   = (g.FFT \ Xk)
    Yx   = (g.FFT \ Yk)
    rtnk = g.FFT * (Xx .* Yx)
    return  rtnk
end

""" (Xk, Yk, filtk, g) → ∫dk̲ im(k+l)X_{k+l} * im(-k)Y_{-k}  """
function ∂conv(Xk, Yk, g)
    rtnk = zeros(Complex{Float64},size(Xk))
    for p in 1:2
        Xx     = g.FFT \ (im .* g.k[p] .* Xk)
        Yx     = g.FFT \ (im .* g.k[p] .* Yk)
        rtnk .+= g.FFT * (Xx .* Yx)
    end
    return  rtnk
end




""" (Xk, Yk, g) → ∫dk̲ sin(2φ_{k+l} - 2φ_{k}) * (X_{k+l} * CX_{k+l}) * conj(Y_k * CY_{k})
where CX and CY are real """
function sin_conv(Xk, Yk, CXk, CYk, g)
    φ_k     = angle(g.k[1] + im * g.k[2])
    expi2φk = exp(im * 2 * φ_k)
    hatXx   = g.FFT \ (squash(Xk .* CXk) .* expi2φk)
    hatYx   = g.FFT \ (squash(Yk .* CYk) .* expi2φk)
    return  g.FFT * imag(hatXx .* conj(hatYx))
end




""" (Xk, Yk, g) → ∫dk̲ sin²(2φ_{k+l} - 2φ_{k}) * (X_{k+l} * CX_{k+l}) * conj(Y_k * CY_{k})
where CX and CY are real """
function sin2_conv(Xk, Yk, CXk, CYk,  g)
    φ_k     = angle(g.k[1] + im * g.k[2])
    expi4φk = exp(im * 4 * φ_k)
    Xx      = g.FFT \ squash(Xk .* CXk)
    Yx      = g.FFT \ squash(Yk .* CYk)
    checkXx = g.FFT \ (squash(Xk .* CXk) .* expi4φk)
    checkYx = g.FFT \ (squash(Yk .* CYk) .* expi4φk)
    rtn     = g.FFT * (Xx .* Yx - real(checkXx .* conj(checkYx)))
    rtn   .*= 0.5
    return  rtn
end





###############################

# De-lensing by treating observed E as unlensed template (i.e. Blakes method)

###############################

function full_blakes_method(lmin, lmax, dek, dbk, mCls, rng, ρl, g, ϕestk)
    mat_clnbk  = ResB_matpwr(mCls, ρl, g)
    δ1residk   = dbk - estδ1BfromEk(dek, ϕestk, ρl, mCls, g)
    rll        = lminlmax_loglike_marg_profile(rng, δ1residk, lmin, lmax, g, mCls, mat_clnbk)
    return rll
end

# new
function δ1BfromEk(_ek, _ϕk, g ) # was LenB1
	zr = zeros(size(g.r))
    qk, uk,    = PrimordialB.eb2qu(_ek, zr, g)
    δ1qk       = PrimordialB.∂conv(qk, _ϕk, g)
    δ1uk       = PrimordialB.∂conv(uk, _ϕk, g)
    δ1E_fromEk, δ1B_fromEk,   = PrimordialB.qu2eb(δ1qk, δ1uk, g)
    return δ1B_fromEk
end

# new
function estδ1BfromEk(dek, ϕestk, ρl, mCls, g)  # was LenB_cib
    μek  = dek .* squash(mCls.cEEk./(mCls.cEEk + mCls.cEEnoisek))
    μϕk  = ϕestk .* squash(ρl.^2)
    estδ1B_fromEk =  δ1BfromEk(μek, μϕk, g)
    return estδ1B_fromEk
end


function ResB_matpwr(mCls, ρl, g)
    mCls_ee_eff   = mCls.cEEk .* squash(mCls.cEEk ./ (mCls.cEEk + mCls.cEEnoisek))
    mCls_ϕϕ_eff   = mCls.cϕϕk .* squash(ρl.^2)
    mCl_lnbb      = lnB_matpwr(mCls.cEEk, mCls.cϕϕk, g)
    mCl_lnbb_cib  = lnB_matpwr(mCls_ee_eff, mCls_ϕϕ_eff, g)
    mCl_bb_res    = mCl_lnbb - mCl_lnbb_cib
    return mCl_bb_res
end
function lnB_matpwr(mCls_ee, mCls_ϕϕ, g)
    φ_l    = angle(g.k[1] + im * g.k[2])
    cosφ2  = cos(2φ_l)
    sinφ2  = sin(2φ_l)
    cosφ4  = cos(4φ_l)
    sinφ4  = sin(4φ_l)

    unt_t1 = lnB_pwr_helper(g.r.^2 .* mCls_ee,          g.r.^2 .*          mCls_ϕϕ, g, opt=:unit)
    unt_t2 = lnB_pwr_helper(g.r.^2 .* cosφ2 .* mCls_ee, g.r.^2 .* cosφ2 .* mCls_ϕϕ, g, opt=:unit)
    unt_t3 = lnB_pwr_helper(g.r.^2 .* sinφ2 .* mCls_ee, g.r.^2 .* sinφ2 .* mCls_ϕϕ, g, opt=:unit)

    cos_t1 = lnB_pwr_helper(g.r.^2 .* mCls_ee,          g.r.^2 .*          mCls_ϕϕ.* cosφ4, g, opt=:cos).* (-cosφ4)
    cos_t2 = lnB_pwr_helper(g.r.^2 .* cosφ2 .* mCls_ee, g.r.^2 .* cosφ2 .* mCls_ϕϕ.* cosφ4, g, opt=:cos).* (-cosφ4)
    cos_t3 = lnB_pwr_helper(g.r.^2 .* sinφ2 .* mCls_ee, g.r.^2 .* sinφ2 .* mCls_ϕϕ.* cosφ4, g, opt=:cos).* (-cosφ4)
    cos_t4 = lnB_pwr_helper(g.r.^2 .* mCls_ee,          g.r.^2 .*          mCls_ϕϕ.* sinφ4, g, opt=:cos).* (-sinφ4)
    cos_t5 = lnB_pwr_helper(g.r.^2 .* cosφ2 .* mCls_ee, g.r.^2 .* cosφ2 .* mCls_ϕϕ.* sinφ4, g, opt=:cos).* (-sinφ4)
    cos_t6 = lnB_pwr_helper(g.r.^2 .* sinφ2 .* mCls_ee, g.r.^2 .* sinφ2 .* mCls_ϕϕ.* sinφ4, g, opt=:cos).* (-sinφ4)

    sin_t1 = lnB_pwr_helper(g.r.^2 .* mCls_ee,          g.r.^2 .*          mCls_ϕϕ.* sinφ4, g, opt=:sin).* cosφ4
    sin_t2 = lnB_pwr_helper(g.r.^2 .* cosφ2 .* mCls_ee, g.r.^2 .* cosφ2 .* mCls_ϕϕ.* sinφ4, g, opt=:sin).* cosφ4
    sin_t3 = lnB_pwr_helper(g.r.^2 .* sinφ2 .* mCls_ee, g.r.^2 .* sinφ2 .* mCls_ϕϕ.* sinφ4, g, opt=:sin).* cosφ4
    sin_t4 = lnB_pwr_helper(g.r.^2 .* mCls_ee,          g.r.^2 .*          mCls_ϕϕ.* cosφ4, g, opt=:sin).* (-sinφ4)
    sin_t5 = lnB_pwr_helper(g.r.^2 .* cosφ2 .* mCls_ee, g.r.^2 .* cosφ2 .* mCls_ϕϕ.* cosφ4, g, opt=:sin).* (-sinφ4)
    sin_t6 = lnB_pwr_helper(g.r.^2 .* sinφ2 .* mCls_ee, g.r.^2 .* sinφ2 .* mCls_ϕϕ.* cosφ4, g, opt=:sin).* (-sinφ4)

    mCl = 0.25*(unt_t1 + unt_t2 + unt_t3
                +sin_t1 + sin_t2 + sin_t3 + sin_t4 + sin_t5 + sin_t6
                +cos_t1 + cos_t2 + cos_t3 + cos_t4 + cos_t5 + cos_t6)/(2pi)
    return mCl
end
function lnB_pwr_helper(Ak::Array{Float64,2}, Bk::Array{Float64,2}, g; opt = :sin)
    φ4_l    = 4angle(g.k[1] + im * g.k[2])
    Ak_hat  = Ak .* exp(im * φ4_l)
    Bk_hat  = Bk .* exp(im * φ4_l)

    if opt == :unit
        Ax = g.FFT \ Ak
        Bx = g.FFT \ Bk
        return real( g.FFT * (Ax .* conj(Bx)) )
    elseif opt == :sin
        Ax_hat  = g.FFT \ Ak_hat
        Bx_hat  = g.FFT \ Bk_hat
        return real( g.FFT * imag(Ax_hat .* conj(Bx_hat)) )
    elseif opt == :cos
        Ax_hat  = g.FFT \ Ak_hat
        Bx_hat  = g.FFT \ Bk_hat
        return real( g.FFT * real(Ax_hat .* conj(Bx_hat)) )
    else
        println("Wrong Option")
    end
end






##########################################################

# Lensing functions

##############################################################

""" Lense qx, ux:  `rqx, rux = lense(qx, ux, len, g, order=2, qk=g.FFT*qx, uk=g.FFT*ux)` """
function lense{T}(
			qx::Matrix{Float64},
			ux::Matrix{Float64},
			len,
			g::FFTgrid{2,T},
			order::Int64 = 2,
			qk::Matrix{Complex{Float64}} = g.FFT * qx,
			uk::Matrix{Complex{Float64}} = g.FFT * ux,
	)
	rqx, rux  = intlense(qx, ux, len)  # <--- return values
	@inbounds for n in 1:order, α₁ in 0:n
		# kα   = im ^ n .* g.k[1] .^ α₁ .* g.k[2] .^ (n - α₁)
		kα  = intlense_helper1(n, α₁, g.k[1], g.k[2])

		∂α_qx = real(g.FFT \ (kα .* qk))
		∂α_ux = real(g.FFT \ (kα .* uk))
		∂α_qx, ∂α_ux  = intlense(∂α_qx, ∂α_ux, len)

		# xα   = len.rdisplx .^ α₁ .* len.rdisply .^ (n - α₁) ./ factorial(α₁) ./ factorial(n - α₁)
		# rqx += xα .* ∂α_qx
		# rux += xα .* ∂α_ux
		intlense_helper2!(rqx, rux, n, α₁, len.rdisplx, len.rdisply, ∂α_qx, ∂α_ux)
    end
    return rqx, rux
end
function intlense_helper1(n, α₁, k1, k2)
	rtn  = Array(Complex{Float64}, size(k1))
	imn, nmα₁  = im ^ n, n - α₁
	@inbounds @simd for i in eachindex(rtn)
		rtn[i] = complex(imn * k1[i] ^ α₁ * k2[i] ^ nmα₁)
	end
	return rtn
end
function intlense_helper2!(rqx, rux, n, α₁, rx, ry, ∂qx, ∂ux)
	fα₁, fnmα₁, nmα₁ = factorial(α₁), factorial(n - α₁), n - α₁
	@inbounds @simd for i in eachindex(rqx, rux)
		xα      = rx[i] ^ α₁ * ry[i] ^ nmα₁ / fα₁ / fnmα₁
		rqx[i] += xα * ∂qx[i]
		rux[i] += xα * ∂ux[i]
	end
	return nothing
end

function intlense(qx, ux, len)
	rqx  = similar(qx)
	rux  = similar(ux)
    @inbounds for i in eachindex(rqx, rux)
            rqx[i] = qx[len.indrow[i], len.indcol[i]]
            rux[i] = ux[len.indrow[i], len.indcol[i]]
    end
    return rqx, rux
end




##########################################################

# Miscellaneous functions

##############################################################
classy = try
	pyimport(:classy)
end
function class(;ϕscale  = 1.0, ψscale = 0.0, lmax = 6_000, r = 0.005, r0 = 1.0, omega_b = 0.0224567, omega_cdm =0.118489, tau_reio = 0.128312, theta_s = 0.0104098, logA_s_1010 = 3.29056, n_s = 0.968602)
	cosmo = classy[:Class]()
	cosmo[:struct_cleanup]()
	cosmo[:empty]()
	params              = Dict(
	   		"output"        => "tCl, pCl, lCl",
	   		"modes"         => "s,t",
	   		"lensing"       => "yes",
			"l_max_scalars" => lmax + 500,
			"l_max_tensors" => lmax + 500, #lmax + 500,
	      	"omega_b"       => omega_b,
	    	"omega_cdm"     => omega_cdm,
	      	"tau_reio"      => tau_reio,
	      	"100*theta_s"   => 100*theta_s,
	      	"ln10^{10}A_s"  => logA_s_1010,
	      	"n_s"           => n_s,
			"r"             => r,
	        #"k_pivot"      => 0.05,
			#"k_step_trans" => 0.1, # 0.01 for super high resolution
	   		#"l_linstep"    => 10, # 1 for super high resolution
	   		)
	cosmo[:set](params)
	cosmo[:compute]()
	cls_ln              = cosmo[:lensed_cl](lmax)
	cls                 = cosmo[:raw_cl](lmax)
	rtn                 = Dict{Symbol, Array{Float64,1}}(
			:ell        => cls["ell"],
			:ln_tt      => cls_ln["tt"] * (10^6 * cosmo[:T_cmb]()) ^ 2,
			:ln_ee      => cls_ln["ee"] * (10^6 * cosmo[:T_cmb]()) ^ 2,
			:ln_bb      => cls_ln["bb"] * (10^6 * cosmo[:T_cmb]()) ^ 2,
			:ln_te      => cls_ln["te"] * (10^6 * cosmo[:T_cmb]()) ^ 2,
			:ln_tϕ      => cls_ln["tp"] * (10^6 * cosmo[:T_cmb]()),
			:tt 	    => cls["tt"] * (10^6 * cosmo[:T_cmb]()) ^ 2,
			:ee 	    => cls["ee"] * (10^6 * cosmo[:T_cmb]()) ^ 2,
			:bb 	    => cls["bb"] * (10^6 * cosmo[:T_cmb]()) ^ 2,
			:te 	    => cls["te"] * (10^6 * cosmo[:T_cmb]()) ^ 2,
			:tϕ 	    => cls["tp"] * (10^6 * cosmo[:T_cmb]()),
			:ϕϕ         => ϕscale.*cls["pp"],
			:ϕψ         => 0.0.*cls["pp"],
			:ψψ         => ψscale.*cls["pp"],
			:bb0        => cls["bb"] * (10^6 * cosmo[:T_cmb]()) ^ 2 * (r0 / r),
		)
	return rtn
end




""" Convert qu to eb:  `ek, bk, ex, bx = qu2eb(qk, uk, g)` """
function qu2eb(qk, uk, g)
	φ2_l = 2angle(g.k[1] + im * g.k[2])
	ek   = - qk .* cos(φ2_l) - uk .* sin(φ2_l)
	bk   =   qk .* sin(φ2_l) - uk .* cos(φ2_l)
	ex   = real(g.FFT \ ek)
	bx   = real(g.FFT \ bk)
	return ek, bk, ex, bx
end


""" Convert eb to qu: `qk, uk, qx, ux = eb2qu(ek, bk, g)` """
function eb2qu(ek, bk, g)
	φ2_l = 2angle(g.k[1] + im * g.k[2])
	qk   = - ek .* cos(φ2_l) + bk .* sin(φ2_l)
	uk   = - ek .* sin(φ2_l) - bk .* cos(φ2_l)
    qx   = real(g.FFT \ qk)
	ux   = real(g.FFT \ uk)
	return qk, uk, qx, ux
end



function ekbk2_lense_cex_sex_cbx_sbx(ek, bk, len, g, order::Int64 = 7)  # scaled b mode to 1
	cex, sex, cbx, sbx = ekbk2_cex_sex_cbx_sbx(ek, bk, g)
	ln_cex, ln_sex         = PrimordialB.lense(cex, sex, len, g, order)
	ln_sbx, ln_cbx         = PrimordialB.lense(sbx, cbx, len, g, order)
	return ln_cex, ln_sex, ln_cbx, ln_sbx
end

function ekbk2_cex_sex_cbx_sbx(ek, bk, g)
	zr = zeros(size(g.r))
	cek, sek, cex, sex = PrimordialB.eb2qu(ek, zr, g)
	sbk, cbk, sbx, cbx = PrimordialB.eb2qu(zr, bk, g)
	return cex, sex, cbx, sbx
end



function pixelsize_to_muK(pixel_side_length, nside)
    fsky_percent = 100 .* (pixel_side_length .* nside .* π ./ (180*60)) .^2 ./ (4π)
    μKrange      = √(fsky_percent ./ 4.0)
    return μKrange
end



""" kbins, rtnk  = radial_power(fk, smooth, g) """
function radial_power(fk, smooth, g)
	rtnk = Float64[]
	dk = g.deltk
	dm =  length(g.x)
	kbins = collect((smooth*dk):(smooth*dk):(g.nyq))
	for wavenumber in kbins
		indx = (wavenumber-smooth*dk) .< g.r .<= (wavenumber+smooth*dk)
		push!(rtnk, sum(abs2(fk[indx]).* (dk.^dm)) / sum(indx))
	end
	return kbins, rtnk
end


""" kbins, rtnk  = radial_sum(fk, smooth, g) """
function radial_ave(fk, smooth, g)
	rtnk  = Float64[]
	dk    = g.deltk
	kbins = collect((smooth*dk):(smooth*dk):(g.nyq))
	for wavenumber in kbins
		indx = (wavenumber-smooth*dk) .< g.r .<= (wavenumber+smooth*dk)
		push!(rtnk, sum(fk[indx]) / sum(indx))
	end
	return kbins, rtnk
end




# -------- Simulate a mean zero Gaussian random field in the pixel domain given a spectral density.
function sim_xk{dm, T}(cXXk::Array{Float64,dm}, g::FFTgrid{dm, T})
	wx, wk = white_wx_wk(g)
	zk = √(cXXk) .* wk
	zx = real(g.FFT \ zk)
	return zx, zk
end

# ----- white noise
function white_wx_wk{dm, T}(g::FFTgrid{dm, T})
	dx  = g.deltx ^ dm
	wx = randn(size(g.r)) ./ √(dx)
	wk = g.FFT * wx
	return wx, wk
end


squash{T<:Number}(x::T)  = isnan(abs(x)) ? zero(T) : isfinite(abs(x)) ? x : zero(T)
function squash!{dm,T}(x::Array{T,dm}, mask::BitArray{dm}=trues(size(x)))
	@inbounds @simd for i in eachindex(x)
		if !mask[i] || !isfinite(x[i]) || isnan(x[i])
			x[i] = zero(T)
		end
	end
	return x
end
function squash{dm,T}(x::Array{T,dm}, mask::BitArray{dm}=trues(size(x)))
	y = copy(x)
	return squash!(y, mask)
end






end # module
