# used by ... scripts/Sigma_detection.jl






function get_spectra(μKarcminT, beamFWHM, period, nside)
    g = FFTgrid(Pix{period, nside}())
    if !isfile("src/saved_cls/cls.jld") | !isfile("src/saved_cls/cls_no_B.jld")
        cls  = PrimordialB.class(
            r           = 0.005,
            r0          = 1.0,
            lmax        = 15_000,
            omega_b     = 0.0224567,
            omega_cdm   = 0.118489,
            tau_reio    = 0.128312,
            theta_s     = 0.0104098,
            logA_s_1010 = 3.29056,
            n_s         = 0.968602
            )
        cls_no_B  = PrimordialB.class(
            r           = 1e-15,
            r0          = 1.0,
            lmax        = 15_000,
            omega_b     = 0.0224567,
            omega_cdm   = 0.118489,
            tau_reio    = 0.128312,
            theta_s     = 0.0104098,
            logA_s_1010 = 3.29056,
            n_s         = 0.968602
            )
        JLD.@save "src/saved_cls/cls.jld" cls
        JLD.@save "src/saved_cls/cls_no_B.jld" cls_no_B
    end
    cls = JLD.load("src/saved_cls/cls.jld", "cls")
    cls_no_B = JLD.load("src/saved_cls/cls_no_B.jld", "cls_no_B")
    σEEarcmin = √2 * μKarcminT
    σBBarcmin = √2 * μKarcminT
    σEErad    = σEEarcmin * (π / 180 / 60)
    σBBrad    = σBBarcmin * (π / 180 / 60)
    mCls      = MatrixCls(g, cls; σEErad = σEErad, σBBrad = σBBrad, beamFWHM=beamFWHM)
    return cls, mCls, cls_no_B
end




function get_radKEB(mCls, period, nside)
    g = FFTgrid(Pix{period, nside}())
    if !isfile("src/saved_cls/radKX_$(nside)_$(period).jld")
        ugridE   = vcat(0.0, linspace((π/180/60/4).^(1/2), (g.period).^(1/2), 4000) .^ 3)
        cψEk    = mCls.cEEk ./ (g.r.^4)
        cψEk[1] = 0.0
        vKE2 = vecϕKn(2, cψEk, ugridE, g)
        vKE3 = vecϕKn(3, cψEk, ugridE, g)
        vKE4 = vecϕKn(4, cψEk, ugridE, g)
        vKE5 = vecϕKn(5, cψEk, ugridE, g)
        vKE6 = vecϕKn(6, cψEk, ugridE, g)

        ugridB  = vcat(0.0, linspace((π/180/60/4).^(1/2), (g.period).^(1/2), 4000) .^ 3)
        cψBk    = mCls.cBB0k ./ (g.r.^4)
        cψBk[1] = 0.0
        vKB2 = vecϕKn(2, cψBk, ugridB, g)
        vKB3 = vecϕKn(3, cψBk, ugridB, g)
        vKB4 = vecϕKn(4, cψBk, ugridB, g)
        vKB5 = vecϕKn(5, cψBk, ugridB, g)
        vKB6 = vecϕKn(6, cψBk, ugridB, g)

        JLD.@save "src/saved_cls/radKX_$(nside)_$(period).jld" ugridE ugridB vKE2 vKE3 vKE4 vKE5 vKE6 vKB2 vKB3 vKB4 vKB5 vKB6
    end
    d =  JLD.load("src/saved_cls/radKX_$(nside)_$(period).jld")
    radKE = hcat(d["ugridE"], d["vKE2"], d["vKE3"], d["vKE4"], d["vKE5"], d["vKE6"])
    radKB = hcat(d["ugridB"], d["vKB2"], d["vKB3"], d["vKB4"], d["vKB5"], d["vKB6"])

    return radKE, radKB
end





#############################################################
function get_N0ℓAℓ(cls_no_B, phiest_pixel_side_length, phiest_nside, phiest_μKarcminT, phiest_beamFWHM, period, nside)
    phiest_period = phiest_pixel_side_length*phiest_nside*π/(180*60)
    gnew = FFTgrid(Pix{phiest_period, phiest_nside}())
    mask = (gnew.r .> 0.95*gnew.nyq) # | (g.r .< lmin # quadratic estimate mask
    phiest_σEErad    = √2 * phiest_μKarcminT * (π / 180 / 60)
    phiest_σBBrad    = √2 * phiest_μKarcminT * (π / 180 / 60)
    mCls_no_B  = MatrixCls(gnew, cls_no_B; σEErad = phiest_σEErad, σBBrad = phiest_σBBrad, beamFWHM=phiest_beamFWHM)
    cEE     = mCls_no_B.cEEk
    cBB     = mCls_no_B.cBBk
    cEElen  = cls_to_cXXk(cls_no_B[:ell], cls_no_B[:ln_ee], gnew.r)
    cBBlen  = cls_to_cXXk(cls_no_B[:ell], cls_no_B[:ln_bb], gnew.r)
    cEEobs  = cEElen + mCls_no_B.cEEnoisek
    cBBobs  = cBBlen + mCls_no_B.cBBnoisek
    cEEobs_mask = copy(cEEobs)
    cBBobs_mask = copy(cBBobs)
    cEEobs_mask[mask] = Inf
    cBBobs_mask[mask] = Inf
    cEBlen  = zeros(cBBobs)
    Aℓ  = getAℓ(cEE, cBB, cEE, cBB, cEEobs_mask, cBBobs_mask, gnew)
    N0ℓ = getN0ℓ(cEE, cBB, cEEobs_mask, cBBobs_mask, cEEobs, cBBobs_mask, cEBlen, gnew, Aℓ)
    g = FFTgrid(Pix{period, nside}())
    N0ℓ_save, Aℓ_save = zeros(size(g.r)), zeros(size(g.r))
    for i in eachindex(g.k[1])
        _, indx = findmin((g.k[1][i] .- gnew.k[1]).^2 + (g.k[2][i] .- gnew.k[2]).^2)
        N0ℓ_save[i] = N0ℓ[indx]
        Aℓ_save[i] = Aℓ[indx]
    end

    mCls_cnϕk            = copy(N0ℓ_save)
    cnϕk_mask            = (g.r .< 1) | (N0ℓ_save .<= 0.0)
    mCls_cnϕk[cnϕk_mask] = Inf

    N0ℓ_save, Aℓ_save, mCls_cnϕk, cnϕk_mask
end








#  .oooooo..o  o8o                                oooo                .    o8o
# d8P'    `Y8  `"'                                `888              .o8    `"'
# Y88bo.      oooo  ooo. .oo.  .oo.   oooo  oooo   888   .oooo.   .o888oo oooo   .ooooo.  ooo. .oo.    .oooo.o
#  `"Y8888o.  `888  `888P"Y88bP"Y88b  `888  `888   888  `P  )88b    888   `888  d88' `88b `888P"Y88b  d88(  "8
#      `"Y88b  888   888   888   888   888   888   888   .oP"888    888    888  888   888  888   888  `"Y88b.
# oo     .d8P  888   888   888   888   888   888   888  d8(  888    888 .  888  888   888  888   888  o.  )88b
# 8""88888P'  o888o o888o o888o o888o  `V88V"V8P' o888o `Y888""8o   "888" o888o `Y8bod8P' o888o o888o 8""888P'
#
#
#

function full_quadratic_B_delenser_qxux(lmin, lmax, dncex, dnsex, dsbx, dcbx,  mCls, rng, rng_true, ρl, ϕestk, period, nside)
    g = FFTgrid(Pix{period, nside}())
    hermitianhalf = (g.k[1] .> 0) | ((g.k[1] .== 0.0) & (g.k[2] .> 0.0))
    lrestriction  = (lmin .<= g.r .<= lmax) & hermitianhalf
    mat_clnbk       = ResB_matpwr(mCls, ρl, g)

    cBB0kmax    = mCls.cBB0k[lrestriction]
    cNNkmax     = mCls.cBBnoisek[lrestriction] + mat_clnbk[lrestriction]
    δ0          = 1/g.deltk^2

    # for the dncex, dnsex part
    dek_e, dbk_e, = qu2eb(g.FFT * dncex, g.FFT * dnsex, g)
    dbkmax_e      = (dbk_e - estδ1BfromEk(dek_e, ϕestk, ρl, mCls, g))[lrestriction]

    # for the dsbx, dcbx part
    dek_b, dbk_b, = qu2eb(g.FFT * dsbx, g.FFT * dcbx, g)
    dbkmax_b      = (dbk_b - estδ1BfromEk(dek_b, ϕestk, ρl, mCls, g))[lrestriction]

    abs2_dbkmax_e = abs2(dbkmax_e)
    abs2_dbkmax_b = abs2(dbkmax_b)
    cross_dbkmax_eb = real(dbkmax_e .* conj(dbkmax_b) + dbkmax_b .* conj(dbkmax_e))

    rng_mat = collect(rng)' .+ zeros(rng_true)
    rng_true_mat = collect(rng_true) .+ zeros(rng_mat)
    rout  = zeros(size(rng_mat))
    for j in eachindex(rng_mat)
        cddkmax  = (rng_mat[j]*cBB0kmax + cNNkmax)
        rout[j]  = - sum((abs2_dbkmax_e + rng_true_mat[j] .* abs2_dbkmax_b + √(rng_true_mat[j]) .* cross_dbkmax_eb) ./ cddkmax) / δ0
        rout[j] += - sum(log(cddkmax))
    end
    return rout
end



################ determines how the simulation are generated ##################


# len_E = E + (δ1E_fromE + δ2E_fromE + ⋯) + (δ1E_fromB + δ2E_fromB + ⋯)
# len_B = B + (δ1B_fromE + δ2B_fromE + ⋯) + (δ1B_fromB + δ2B_fromB + ⋯)
function lense_unitr_sc_allorder(ek, bk, len, g, order::Int64 = 7)
    b1k = copy(bk) # Note: input bk must have r = 1.0 for this to work right !!!
    ln_cex_rtn, ln_sex_rtn, ln_cb1x_rtn, ln_sb1x_rtn  = PrimordialB.ekbk2_lense_cex_sex_cbx_sbx(ek, b1k, len, g, order)
	return ln_cex_rtn, ln_sex_rtn, ln_cb1x_rtn, ln_sb1x_rtn
end


# len_E = E + (0) + (0)
# len_B = B + (δ1B_fromE + 0) + (0)
function lense_unitr_sc_null(ek, bk, len, g, order::Int64 = 7)
    b1k = copy(bk) # Note: input bk must have r = 1.0 for this to work right !!!
    zr = zeros(size(g.r))

    ln_cex_rtn = zeros(size(g.r))
    ln_sex_rtn = zeros(size(g.r))
    ln_cb1x_rtn = zeros(size(g.r))
    ln_sb1x_rtn = zeros(size(g.r))

    ############# contributions to len_B ###############
    #--------add ------------------------------------
    # len_E = 0 + (0) + (0)
    # len_B = B + (0) + (0)
    #----------------------------------------------
    _, _, ln_cb1x_tmp, ln_sb1x_tmp    = PrimordialB.ekbk2_cex_sex_cbx_sbx(zr, b1k, g)
    ln_cb1x_rtn += ln_cb1x_tmp
    ln_sb1x_rtn += ln_sb1x_tmp


    # #--------add------------------------------------
    # # len_E = 0 + (0) + (0)
    # # len_B = 0 + (δ1B_fromE + δ2B_fromE + ⋯) + (0)
    # #----------------------------------------------
    # ln_cex_tmp, ln_sex_tmp,  = PrimordialB.ekbk2_lense_cex_sex_cbx_sbx(ek, zr, len, g, order)
    # δ012⋯E_fromEk, δ12⋯B_fromEk,   = PrimordialB.qu2eb(g.FFT * ln_cex_tmp, g.FFT * ln_sex_tmp, g)
    # _, _, ln_cex_tmp, ln_sex_tmp = PrimordialB.eb2qu(zr,  δ12⋯B_fromEk, g)
    # ln_cex_rtn += ln_cex_tmp
    # ln_sex_rtn += ln_sex_tmp

    #----- add ---------------------------------------
    # len_E = 0 + (0) + (0)
    # len_B = 0 + (δ1B_fromE + 0) + (0)
    # --------------------------------------------
    δ1B_fromEk = PrimordialB.δ1BfromEk(ek, len.ϕk, g )
    _, _, ln_cex_tmp, ln_sex_tmp = PrimordialB.eb2qu(zr,  δ1B_fromEk, g)
    ln_cex_rtn += ln_cex_tmp
    ln_sex_rtn += ln_sex_tmp


    ############# contributions to len_E ###############
    #--------add------------------------------------
    # len_E = E + (0) + (0)
    # len_B = 0 + (0) + (0)
    #----------------------------------------------
    ln_cex_tmp, ln_sex_tmp,     = PrimordialB.ekbk2_cex_sex_cbx_sbx(ek, zr, g)
    ln_cex_rtn += ln_cex_tmp
    ln_sex_rtn += ln_sex_tmp

    #
    # #--------add ------------------------------------
    # # len_E = E + (δ1E_fromE + δ2E_fromE + ⋯) + (0)
    # # len_B = 0 + (0) + (0)
    # #----------------------------------------------
    # ln_cex_tmp, ln_sex_tmp,   = PrimordialB.ekbk2_lense_cex_sex_cbx_sbx(ek, zr, len, g, order)
    # δ012⋯E_fromEk, δ12⋯B_fromEk,   = PrimordialB.qu2eb(g.FFT * ln_cex_tmp, g.FFT * ln_sex_tmp, g)
    # _, _, ln_cex_tmp, ln_sex_tmp = PrimordialB.eb2qu(δ012⋯E_fromEk,  zr, g)
    # ln_cex_rtn += ln_cex_tmp
    # ln_sex_rtn += ln_sex_tmp



    # #----- add ---------------------------------------
    # # len_E = 0 + (δ1E_fromE) + (0)
    # # len_B = 0 + (0) + (0)
    # # --------------------------------------------
    # qk, uk,    = PrimordialB.eb2qu(ek, zr, g)
    # δ1qk       = PrimordialB.∂conv(qk, len.ϕk, g)
    # δ1uk       = PrimordialB.∂conv(uk, len.ϕk, g)
    # δ1E_fromEk, δ1B_fromEk,   = PrimordialB.qu2eb(δ1qk, δ1uk, g)
    # ln_cex_tmp, ln_sex_tmp,   = PrimordialB.ekbk2_cex_sex_cbx_sbx(δ1E_fromEk, zr, g)
    # ln_cex_rtn += ln_cex_tmp
    # ln_sex_rtn += ln_sex_tmp


    return ln_cex_rtn, ln_sex_rtn, ln_cb1x_rtn, ln_sb1x_rtn
end





# All orders lensing
# len_E = E + (δ1E_fromE + δ2E_fromE + ⋯) + (0)
# len_B = B + (δ1B_fromE + 0) + (0)
function lense_unitr_sc_allEfromE(ek, bk, len, g, order::Int64 = 7)
    b1k = copy(bk) # Note: input bk must have r = 1.0 for this to work right !!!
    zr = zeros(size(g.r))

    ln_cex_rtn = zeros(size(g.r))
    ln_sex_rtn = zeros(size(g.r))
    ln_cb1x_rtn = zeros(size(g.r))
    ln_sb1x_rtn = zeros(size(g.r))

    ############# contributions to len_B ###############
    #--------add ------------------------------------
    # len_E = 0 + (0) + (0)
    # len_B = B + (0) + (0)
    #----------------------------------------------
    _, _, ln_cb1x_tmp, ln_sb1x_tmp    = PrimordialB.ekbk2_cex_sex_cbx_sbx(zr, b1k, g)
    ln_cb1x_rtn += ln_cb1x_tmp
    ln_sb1x_rtn += ln_sb1x_tmp


    # #--------add------------------------------------
    # # len_E = 0 + (0) + (0)
    # # len_B = 0 + (δ1B_fromE + δ2B_fromE + ⋯) + (0)
    # #----------------------------------------------
    # ln_cex_tmp, ln_sex_tmp,  = PrimordialB.ekbk2_lense_cex_sex_cbx_sbx(ek, zr, len, g, order)
    # δ012⋯E_fromEk, δ12⋯B_fromEk,   = PrimordialB.qu2eb(g.FFT * ln_cex_tmp, g.FFT * ln_sex_tmp, g)
    # _, _, ln_cex_tmp, ln_sex_tmp = PrimordialB.eb2qu(zr,  δ12⋯B_fromEk, g)
    # ln_cex_rtn += ln_cex_tmp
    # ln_sex_rtn += ln_sex_tmp

    #----- add ---------------------------------------
    # len_E = 0 + (0) + (0)
    # len_B = 0 + (δ1B_fromE + 0) + (0)
    # --------------------------------------------
    δ1B_fromEk = PrimordialB.δ1BfromEk(ek, len.ϕk, g )
    _, _, ln_cex_tmp, ln_sex_tmp = PrimordialB.eb2qu(zr,  δ1B_fromEk, g)
    ln_cex_rtn += ln_cex_tmp
    ln_sex_rtn += ln_sex_tmp


    ############# contributions to len_E ###############
    # #--------add------------------------------------
    # # len_E = E + (0) + (0)
    # # len_B = 0 + (0) + (0)
    # #----------------------------------------------
    # ln_cex_tmp, ln_sex_tmp,     = PrimordialB.ekbk2_cex_sex_cbx_sbx(ek, zr, g)
    # ln_cex_rtn += ln_cex_tmp
    # ln_sex_rtn += ln_sex_tmp


    #--------add ------------------------------------
    # len_E = E + (δ1E_fromE + δ2E_fromE + ⋯) + (0)
    # len_B = 0 + (0) + (0)
    #----------------------------------------------
    ln_cex_tmp, ln_sex_tmp,   = PrimordialB.ekbk2_lense_cex_sex_cbx_sbx(ek, zr, len, g, order)
    δ012⋯E_fromEk, δ12⋯B_fromEk,   = PrimordialB.qu2eb(g.FFT * ln_cex_tmp, g.FFT * ln_sex_tmp, g)
    _, _, ln_cex_tmp, ln_sex_tmp = PrimordialB.eb2qu(δ012⋯E_fromEk,  zr, g)
    ln_cex_rtn += ln_cex_tmp
    ln_sex_rtn += ln_sex_tmp



    # #----- add ---------------------------------------
    # # len_E = 0 + (δ1E_fromE) + (0)
    # # len_B = 0 + (0) + (0)
    # # --------------------------------------------
    # qk, uk,    = PrimordialB.eb2qu(ek, zr, g)
    # δ1qk       = PrimordialB.∂conv(qk, len.ϕk, g)
    # δ1uk       = PrimordialB.∂conv(uk, len.ϕk, g)
    # δ1E_fromEk, δ1B_fromEk,   = PrimordialB.qu2eb(δ1qk, δ1uk, g)
    # ln_cex_tmp, ln_sex_tmp,   = PrimordialB.ekbk2_cex_sex_cbx_sbx(δ1E_fromEk, zr, g)
    # ln_cex_rtn += ln_cex_tmp
    # ln_sex_rtn += ln_sex_tmp


    return ln_cex_rtn, ln_sex_rtn, ln_cb1x_rtn, ln_sb1x_rtn
end







# len_E = E + (0) + (0)
# len_B = B + (δ1B_fromE + δ2B_fromE + ⋯) + (0)
function lense_unitr_sc_allBfromE(ek, bk, len, g, order::Int64 = 7)
    b1k = copy(bk) # Note: input bk must have r = 1.0 for this to work right !!!
    zr = zeros(size(g.r))

    ln_cex_rtn = zeros(size(g.r))
    ln_sex_rtn = zeros(size(g.r))
    ln_cb1x_rtn = zeros(size(g.r))
    ln_sb1x_rtn = zeros(size(g.r))

    ############# contributions to len_B ###############
    #--------add ------------------------------------
    # len_E = 0 + (0) + (0)
    # len_B = B + (0) + (0)
    #----------------------------------------------
    _, _, ln_cb1x_tmp, ln_sb1x_tmp    = PrimordialB.ekbk2_cex_sex_cbx_sbx(zr, b1k, g)
    ln_cb1x_rtn += ln_cb1x_tmp
    ln_sb1x_rtn += ln_sb1x_tmp


    #--------add------------------------------------
    # len_E = 0 + (0) + (0)
    # len_B = 0 + (δ1B_fromE + δ2B_fromE + ⋯) + (0)
    #----------------------------------------------
    ln_cex_tmp, ln_sex_tmp,  = PrimordialB.ekbk2_lense_cex_sex_cbx_sbx(ek, zr, len, g, order)
    δ012⋯E_fromEk, δ12⋯B_fromEk,   = PrimordialB.qu2eb(g.FFT * ln_cex_tmp, g.FFT * ln_sex_tmp, g)
    _, _, ln_cex_tmp, ln_sex_tmp = PrimordialB.eb2qu(zr,  δ12⋯B_fromEk, g)
    ln_cex_rtn += ln_cex_tmp
    ln_sex_rtn += ln_sex_tmp

    # #----- add ---------------------------------------
    # # len_E = 0 + (0) + (0)
    # # len_B = 0 + (δ1B_fromE + 0) + (0)
    # # --------------------------------------------
    # δ1B_fromEk = PrimordialB.δ1BfromEk(ek, len.ϕk, g )
    # _, _, ln_cex_tmp, ln_sex_tmp = PrimordialB.eb2qu(zr,  δ1B_fromEk, g)
    # ln_cex_rtn += ln_cex_tmp
    # ln_sex_rtn += ln_sex_tmp


    ############# contributions to len_E ###############
    #--------add------------------------------------
    # len_E = E + (0) + (0)
    # len_B = 0 + (0) + (0)
    #----------------------------------------------
    ln_cex_tmp, ln_sex_tmp,     = PrimordialB.ekbk2_cex_sex_cbx_sbx(ek, zr, g)
    ln_cex_rtn += ln_cex_tmp
    ln_sex_rtn += ln_sex_tmp

    #
    # #--------add ------------------------------------
    # # len_E = E + (δ1E_fromE + δ2E_fromE + ⋯) + (0)
    # # len_B = 0 + (0) + (0)
    # #----------------------------------------------
    # ln_cex_tmp, ln_sex_tmp,   = PrimordialB.ekbk2_lense_cex_sex_cbx_sbx(ek, zr, len, g, order)
    # δ012⋯E_fromEk, δ12⋯B_fromEk,   = PrimordialB.qu2eb(g.FFT * ln_cex_tmp, g.FFT * ln_sex_tmp, g)
    # _, _, ln_cex_tmp, ln_sex_tmp = PrimordialB.eb2qu(δ012⋯E_fromEk,  zr, g)
    # ln_cex_rtn += ln_cex_tmp
    # ln_sex_rtn += ln_sex_tmp



    # #----- add ---------------------------------------
    # # len_E = 0 + (δ1E_fromE) + (0)
    # # len_B = 0 + (0) + (0)
    # # --------------------------------------------
    # qk, uk,    = PrimordialB.eb2qu(ek, zr, g)
    # δ1qk       = PrimordialB.∂conv(qk, len.ϕk, g)
    # δ1uk       = PrimordialB.∂conv(uk, len.ϕk, g)
    # δ1E_fromEk, δ1B_fromEk,   = PrimordialB.qu2eb(δ1qk, δ1uk, g)
    # ln_cex_tmp, ln_sex_tmp,   = PrimordialB.ekbk2_cex_sex_cbx_sbx(δ1E_fromEk, zr, g)
    # ln_cex_rtn += ln_cex_tmp
    # ln_sex_rtn += ln_sex_tmp


    return ln_cex_rtn, ln_sex_rtn, ln_cb1x_rtn, ln_sb1x_rtn
end







# ooooooooo.                          o8o                         .    o8o                                                               .o8
# `888   `Y88.                        `"'                       .o8    `"'                                                              "888
#  888   .d88' oooo d8b  .ooooo.     oooo  .ooooo.   .ooooo.  .o888oo oooo   .ooooo.  ooo. .oo.        ooo. .oo.  .oo.    .ooooo.   .oooo888   .ooooo.   .oooo.o           .ooooo.   .ooooo.  oooo    ooo
#  888ooo88P'  `888""8P d88' `88b    `888 d88' `88b d88' `"Y8   888   `888  d88' `88b `888P"Y88b       `888P"Y88bP"Y88b  d88' `88b d88' `888  d88' `88b d88(  "8          d88' `"Y8 d88' `88b  `88.  .8'
#  888          888     888   888     888 888ooo888 888         888    888  888   888  888   888        888   888   888  888   888 888   888  888ooo888 `"Y88b.           888       888   888   `88..8'
#  888          888     888   888     888 888    .o 888   .o8   888 .  888  888   888  888   888        888   888   888  888   888 888   888  888    .o o.  )88b .o.      888   .o8 888   888    `888'
# o888o        d888b    `Y8bod8P'     888 `Y8bod8P' `Y8bod8P'   "888" o888o `Y8bod8P' o888o o888o      o888o o888o o888o `Y8bod8P' `Y8bod88P" `Y8bod8P' 8""888P' Y8P      `Y8bod8P' `Y8bod8P'     `8'
#                                     888                                                                                                                         '
#                                 .o. 88P
#                                 `Y888P
#
#
#






function orth(v; cut = 0.0)
    U, S, V = svd(v)
    return U[:, S .> cut]
end







function quadratic_B_delenser(dqx, dux, ϕestk, ρl, mCls, period, nside)
    g = FFTgrid(Pix{period, nside}())
    dqk = g.FFT * dqx
    duk = g.FFT * dux
    dek, dbk,  = qu2eb(dqk, duk, g)
    δ1residk   = dbk - estδ1BfromEk(dek, ϕestk, ρl, mCls, g)
    return δ1residk
end



# ooooo                                      oooo       oooo   o8o  oooo                  oooo   o8o  oooo                                  .o8
# `888'                                      `888       `888   `"'  `888                  `888   `"'  `888                                 "888
#  888          .ooooo.   .ooooo.   .oooo.    888        888  oooo   888  oooo   .ooooo.   888  oooo   888 .oo.    .ooooo.   .ooooo.   .oooo888
#  888         d88' `88b d88' `"Y8 `P  )88b   888        888  `888   888 .8P'   d88' `88b  888  `888   888P"Y88b  d88' `88b d88' `88b d88' `888
#  888         888   888 888        .oP"888   888        888   888   888888.    888ooo888  888   888   888   888  888   888 888   888 888   888
#  888       o 888   888 888   .o8 d8(  888   888        888   888   888 `88b.  888    .o  888   888   888   888  888   888 888   888 888   888
# o888ooooood8 `Y8bod8P' `Y8bod8P' `Y888""8o o888o      o888o o888o o888o o888o `Y8bod8P' o888o o888o o888o o888o `Y8bod8P' `Y8bod8P' `Y8bod88P"
#
#
#
#






function loglike_r_rtrue(rng, rng_true, qudata_en, qudata_b, D, V, Db, Vb)
    #qudata_en = vcat(dncex, dnsex)
    #qudata_b = vcat(dcbx, dsbx)
    VLden   = At_mul_B(V, At_mul_B(Vb, qudata_en) ./ √(Db))
    VLdb    = At_mul_B(V, At_mul_B(Vb, qudata_b) ./ √(Db))
    Vlden_sq = VLden .^ 2
    Vldb_sq  = VLdb .^ 2
    Vlden_b  = VLdb .* VLden
    rng_mat = collect(rng)' .+ zeros(rng_true)
    rng_true_mat = collect(rng_true) .+ zeros(rng_mat)
    rout  = zeros(size(rng_mat))
    for j in eachindex(rng_mat)
        for k in eachindex(Vlden_sq)
            rout[j] += - 0.5 * (Vlden_sq[k] + rng_true_mat[j]*Vldb_sq[k] + 2*sqrt(rng_true_mat[j])*Vlden_b[k])/(D[k] + rng_mat[j]) - (0.5*log(D[k] + rng_mat[j]))
        end
    end
    return rout
end
