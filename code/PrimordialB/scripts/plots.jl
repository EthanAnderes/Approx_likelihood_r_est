#' plots.jl: make plots for the paper
#' ================================================

#' To generate a document from the code output, run the following at the Julia terminal:
#' `using Weave` then
#' `weave("scripts/plots.jl", doctype = "github", out_path = "scripts")`

#= ------------ for loading a simulation run from a remote server
local_filename   = "out_figure4.jld"
remote_filename  = "out_figure4.jld"
remote_machine   =  "poisson" #
@eval run(`
    scp
    anderes@$remote_machine.ucdavis.edu:~/PrimordialB.jl/scripts/$remote_filename
    /Users/ethananderes/Dropbox/PrimordialB/scripts/$local_filename
`)
=# #---------------------------------------------



#   .oooooo.                   o8o            oooo                        oooo                .
#  d8P'  `Y8b                  `"'            `888                        `888              .o8
# 888      888    oooo  oooo  oooo   .ooooo.   888  oooo       oo.ooooo.   888   .ooooo.  .o888oo
# 888      888    `888  `888  `888  d88' `"Y8  888 .8P'         888' `88b  888  d88' `88b   888
# 888      888     888   888   888  888        888888.          888   888  888  888   888   888
# `88b    d88b     888   888   888  888   .o8  888 `88b.        888   888  888  888   888   888 .
#  `Y8bood8P'Ybd'  `V88V"V8P' o888o `Y8bod8P' o888o o888o       888bod8P' o888o `Y8bod8P'   "888"
#                                                               888
#                                                              o888o
#
#


#= Quick plot for viewing results

using JLD, PyPlot
loadpath = "scripts/saved_run_0x3f77869c073f2346/out_figure2c.jld"
rng_true = load(loadpath, "rng_true")
rest_B_1 = load(loadpath, "rest_B")
rest_1   = load(loadpath, "rest_modes")
bias_B_1  = mean((rest_B_1 .- rng_true),2)
bias_L_1  = mean((rest_1 .- rng_true),2)
standard_error_B_1 = sqrt.(mean((rest_B_1 .- rng_true).^2, 2))
standard_error_L_1 = sqrt.(mean((rest_1 .- rng_true).^2, 2))


fig, ax = plt[:subplots](1,1, figsize=(5,3.7), sharex=true,  sharey=true)
fontsize = 10

ax[:plot](rng_true, rng_true./standard_error_B_1, label = "quadratic delenser")
ax[:plot](rng_true, rng_true./standard_error_L_1, label = "likelihood approx")
ax[:set_ylabel](L"$\frac{r}{\sigma(r)}$       ", fontsize=fontsize+4,  rotation="horizontal")
ax[:set_xlabel](L"$r$", fontsize=fontsize)
ax[:set_xticks]([.01, .03, .05, .07, .09])
ax[:tick_params](axis="y", labelsize = fontsize-2)
ax[:tick_params](axis="x", labelsize = fontsize-2)
ax[:legend](loc="best", fontsize=fontsize)
ax[:grid](true)
ax[:set_title](L"\mu K arcmin = 1.0", fontsize=fontsize)


fig, ax = plt[:subplots](1,1, figsize=(5,3.7), sharex=true,  sharey=true)
fontsize = 10

ax[:plot](rng_true, bias_B_1./rng_true, label = "quadratic delenser")
ax[:plot](rng_true, bias_L_1./rng_true, label = "likelihood approx")
ax[:set_ylabel](L"$\frac{Bias(r)}{r}$       ", fontsize=fontsize+4,  rotation="horizontal")
ax[:set_xlabel](L"$r$", fontsize=fontsize)
ax[:set_xticks]([.01, .03, .05, .07, .09])
ax[:tick_params](axis="y", labelsize = fontsize-2)
ax[:tick_params](axis="x", labelsize = fontsize-2)
ax[:legend](loc="best", fontsize=fontsize)
ax[:grid](true)
ax[:set_title](L"\mu K arcmin = 1.0", fontsize=fontsize)
ax[:axhline](y=.0, color = "k", linestyle="dashed")


=#








#  .o88o.  o8o                                                   .o
#  888 `"  `"'                                                 o888
# o888oo  oooo   .oooooooo oooo  oooo  oooo d8b  .ooooo.        888
#  888    `888  888' `88b  `888  `888  `888""8P d88' `88b       888
#  888     888  888   888   888   888   888     888ooo888       888
#  888     888  `88bod8P'   888   888   888     888    .o       888
# o888o   o888o `8oooooo.   `V88V"V8P' d888b    `Y8bod8P'      o888o
#               d"     YD
#               "Y88888P'
#
#
#
# Plot of phi estimation noise used in the simulations


using PyPlot, JLD
savefigs = true
savepath = "scripts/saved_run_0x3f77869c073f2346/"



#' f1.pdf : phi N0ℓ for CMB-S3 vrs CMB-S4
#' ================================================
include("src/PrimordialB.jl")
kbins, N0ℓ_a, cϕϕ = let
    nside     = load(savepath*"out_figure2a.jld", "nside")
    period    = load(savepath*"out_figure2a.jld", "period")
    μKarcminT = load(savepath*"out_figure2a.jld", "μKarcminT")
    pixel_side_length = load(savepath*"out_figure2a.jld", "pixel_side_length")

    phiest_pixel_side_length = 2 # fixed
    phiest_nside          = ceil(Int64, pixel_side_length*nside/phiest_pixel_side_length) #nside * 5
    phiest_μKarcminT      = μKarcminT
    beamFWHM              = deg2rad(pixel_side_length/60)
    phiest_beamFWHM       = deg2rad(phiest_pixel_side_length/60)

    N0ℓ       = load(savepath*"out_figure2a.jld", "N0ℓ")
    g = PrimordialB.FFTgrid(PrimordialB.Pix{period, nside}())
    cls, mCls, cls_no_B    = PrimordialB.get_spectra(μKarcminT, beamFWHM, period, nside) # for B survey

    kbins, N0ℓ_tmp  = PrimordialB.radial_ave((g.r.^4) .* N0ℓ, 1, g)
    _, cϕϕ    = PrimordialB.radial_ave((g.r.^4) .* mCls.cϕϕk, 1, g)
    kbins, N0ℓ_tmp, cϕϕ
end
kbins, N0ℓ_b, cϕϕ = let
    nside     = load(savepath*"out_figure2b.jld", "nside")
    period    = load(savepath*"out_figure2b.jld", "period")
    μKarcminT = load(savepath*"out_figure2b.jld", "μKarcminT")
    pixel_side_length = load(savepath*"out_figure2b.jld", "pixel_side_length")

    phiest_pixel_side_length = 2 # fixed
    phiest_nside          = ceil(Int64, pixel_side_length*nside/phiest_pixel_side_length) #nside * 5
    phiest_μKarcminT      = μKarcminT
    beamFWHM              = deg2rad(pixel_side_length/60)
    phiest_beamFWHM       = deg2rad(phiest_pixel_side_length/60)

    N0ℓ       = load(savepath*"out_figure2b.jld", "N0ℓ")
    g = PrimordialB.FFTgrid(PrimordialB.Pix{period, nside}())
    cls, mCls, cls_no_B    = PrimordialB.get_spectra(μKarcminT, beamFWHM, period, nside) # for B survey

    kbins, N0ℓ_tmp  = PrimordialB.radial_ave((g.r.^4) .* N0ℓ, 1, g)
    _, cϕϕ    = PrimordialB.radial_ave((g.r.^4) .* mCls.cϕϕk, 1, g)
    kbins, N0ℓ_tmp, cϕϕ
end
kbins, N0ℓ_c, cϕϕ = let
    nside     = load(savepath*"out_figure2c.jld", "nside")
    period    = load(savepath*"out_figure2c.jld", "period")
    μKarcminT = load(savepath*"out_figure2c.jld", "μKarcminT")
    pixel_side_length = load(savepath*"out_figure2c.jld", "pixel_side_length")

    phiest_pixel_side_length = 2 # fixed
    phiest_nside          = ceil(Int64, pixel_side_length*nside/phiest_pixel_side_length) #nside * 5
    phiest_μKarcminT      = μKarcminT
    beamFWHM              = deg2rad(pixel_side_length/60)
    phiest_beamFWHM       = deg2rad(phiest_pixel_side_length/60)

    N0ℓ       = load(savepath*"out_figure2c.jld", "N0ℓ")
    g = PrimordialB.FFTgrid(PrimordialB.Pix{period, nside}())
    cls, mCls, cls_no_B    = PrimordialB.get_spectra(μKarcminT, beamFWHM, period, nside) # for B survey

    kbins, N0ℓ_tmp  = PrimordialB.radial_ave((g.r.^4) .* N0ℓ, 1, g)
    _, cϕϕ    = PrimordialB.radial_ave((g.r.^4) .* mCls.cϕϕk, 1, g)
    kbins, N0ℓ_tmp, cϕϕ
end






#' Plot
fig, ax = plt[:subplots](1,1, figsize=(4.5,3.5))
fontsize = 10
ax[:plot](kbins, cϕϕ, label = L"$\ell^4 C_\ell^{\phi\phi}$")
ax[:plot](kbins, N0ℓ_a, label = L"\mu K arcmin = 10")
ax[:plot](kbins, N0ℓ_b, label = L"\mu K arcmin = 1")
ax[:plot](kbins, N0ℓ_c, label = L"\mu K arcmin = 0.5")
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:set_title](L"$\phi$ uncertainty", fontsize=fontsize)
ax[:set_xlabel](L"multiple $\ell$", fontsize=fontsize)
ax[:legend](loc="best", fontsize=fontsize-2)
#ax[:set_ylim]([1e-8, 1.5e-6])
ax[:grid](true)


#' Now save
savename  = "f1.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end



# using PyPlot
# loglog(cls[:ell], cls[:ell].^2 .* cls[:tt])
# loglog(cls[:ell], cls[:ell].^2 .* cls[:ee])
# loglog(cls[:ell], cls[:ell].^2 .* cls[:bb])
# #loglog(cls[:ell], cls[:ell].^2 .* cls_no_B[:bb])
# loglog(cls[:ell], cls[:ell].^2 .* cls_no_B[:ln_bb])
# ksl = g.k[2][1:25,1]
# loglog(ksl, ksl.^2 .* mCls.cBBnoisek[1:25,1])
# loglog(ksl, ksl.^2 .* mCls.cBBk[1:25,1])
# loglog(ksl, ksl.^2 .* mCls.cEEk[1:25,1])
# loglog(ksl, ksl.^2 .* mCls.cTTk[1:25,1])



#  .o88o.  o8o                                                   .oooo.
#  888 `"  `"'                                                 .dP""Y88b
# o888oo  oooo   .oooooooo oooo  oooo  oooo d8b  .ooooo.             ]8P'
#  888    `888  888' `88b  `888  `888  `888""8P d88' `88b          .d8P'
#  888     888  888   888   888   888   888     888ooo888        .dP'
#  888     888  `88bod8P'   888   888   888     888    .o      .oP     .o
# o888o   o888o `8oooooo.   `V88V"V8P' d888b    `Y8bod8P'      8888888888
#               d"     YD
#               "Y88888P'
#
#
#
#



using PyPlot, JLD
savefigs = true
savepath = "scripts/saved_run_0x3f77869c073f2346/"



#' f2a.pdf, f2b.pdf, f2c.pdf:A comparison of
#' ================================================

rng_true = load(savepath*"out_figure2a.jld", "rng_true")

# low resolution survey for ϕ est
rest_B_1 = load(savepath*"out_figure2a.jld", "rest_B")
rest_1   = load(savepath*"out_figure2a.jld", "rest_modes")
rest_B_2 = load(savepath*"out_figure2b.jld", "rest_B")
rest_2   = load(savepath*"out_figure2b.jld", "rest_modes")
rest_B_3 = load(savepath*"out_figure2c.jld", "rest_B")
rest_3   = load(savepath*"out_figure2c.jld", "rest_modes")


# compute bias and standard error
bias_B_1  = mean((rest_B_1 .- rng_true),2)
bias_B_2  = mean((rest_B_2 .- rng_true),2)
bias_B_3  = mean((rest_B_3 .- rng_true),2)
bias_L_1  = mean((rest_1 .- rng_true),2)
bias_L_2  = mean((rest_2 .- rng_true),2)
bias_L_3  = mean((rest_3 .- rng_true),2)

standard_error_B_1 = sqrt.(mean((rest_B_1 .- rng_true).^2, 2))
standard_error_B_2 = sqrt.(mean((rest_B_2 .- rng_true).^2, 2))
standard_error_B_3 = sqrt.(mean((rest_B_3 .- rng_true).^2, 2))
standard_error_L_1 = sqrt.(mean((rest_1 .- rng_true).^2, 2))
standard_error_L_2 = sqrt.(mean((rest_2 .- rng_true).^2, 2))
standard_error_L_3 = sqrt.(mean((rest_3 .- rng_true).^2, 2))

maxylab = maximum(vcat(
    rng_true./standard_error_B_1,
    rng_true./standard_error_B_2,
    rng_true./standard_error_B_3,
    rng_true./standard_error_L_1,
    rng_true./standard_error_L_2,
    rng_true./standard_error_L_3
    ))

minbias, maxbias = vcat(
         bias_B_1./standard_error_B_1,
         bias_B_2./standard_error_B_2,
         bias_B_3./standard_error_B_3,
         bias_L_1./standard_error_L_1,
         bias_L_2./standard_error_L_2,
         bias_L_3./standard_error_L_3,
         ) |> x-> (minimum(x), maximum(x))


#' f2a.pdf : plot r / σ(r)
#' ================================================


fig, ax = plt[:subplots](1,3, figsize=(10,3.7), sharex=true,  sharey=true)
fontsize = 10

#  CMB-S3 for ϕ est
ax[1][:set_title](L"\mu K arcmin = 10", fontsize=fontsize)
ax[1][:plot](rng_true, rng_true./standard_error_B_1, label = "quadratic delenser")
ax[1][:plot](rng_true, rng_true./standard_error_L_1, label = "likelihood approx")
ax[1][:set_ylabel](L"$\frac{r}{\sigma(r)}$       ", fontsize=fontsize+4,  rotation="horizontal")
ax[1][:set_xlabel](L"$r$", fontsize=fontsize)
ax[1][:set_xticks]([.01, .03, .05, .07, .09])
ax[1][:tick_params](axis="y", labelsize = fontsize-2)
ax[1][:tick_params](axis="x", labelsize = fontsize-2)
ax[1][:legend](loc="best", fontsize=fontsize)
ax[1][:grid](true)

#  CMB-S4 for ϕ est
ax[2][:set_title](L"\mu K arcmin = 1", fontsize=fontsize)
ax[2][:plot](rng_true, rng_true./standard_error_B_2, label = "quadratic delenser")
ax[2][:plot](rng_true, rng_true./standard_error_L_2, label = "likelihood approx")
ax[2][:set_xlabel](L"$r$", fontsize=fontsize)
ax[2][:tick_params](axis="y", left="off", right="off")#labelsize = fontsize-2, pad = 5, labelright="on")
ax[2][:tick_params](axis="x", labelsize = fontsize-2)
ax[2][:axis]("tight")
ax[2][:grid](true)


#  CMB-S4 for ϕ est
ax[3][:set_title](L"\mu K arcmin = 0.5", fontsize=fontsize)
ax[3][:plot](rng_true, rng_true./standard_error_B_3, label = "quadratic delenser")
ax[3][:plot](rng_true, rng_true./standard_error_L_3, label = "likelihood approx")
ax[3][:set_xlabel](L"$r$", fontsize=fontsize)
ax[3][:tick_params](axis="y", labelsize = fontsize-2, labelright="on")
ax[3][:tick_params](axis="x", labelsize = fontsize-2, pad = 5)
ax[3][:axis]("tight")
ax[3][:grid](true)


#fig[:suptitle]("    Detection", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f2a.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end



#' f2b.pdf : plot  Bias(r) / r
#' ================================================

fig, ax = plt[:subplots](1,3, figsize=(10,3.7), sharex=true,  sharey=true)
fontsize = 10

#  CMB-S3 for ϕ est
ax[1][:set_title](L"\mu K arcmin = 10", fontsize=fontsize)
ax[1][:plot](rng_true, bias_B_1./rng_true, label = "quadratic delenser")
ax[1][:plot](rng_true, bias_L_1./rng_true, label = "likelihood approx")
ax[1][:set_ylabel](L"$\frac{Bias(r)}{r}$        ", fontsize=fontsize+4,  rotation="horizontal")
ax[1][:set_xlabel](L"$r$", fontsize=fontsize)
ax[1][:set_xticks]([.01, .03, .05, .07, .09])
ax[1][:tick_params](axis="y", labelsize = fontsize-2)
ax[1][:tick_params](axis="x", labelsize = fontsize-2)
ax[1][:legend](loc="best", fontsize=fontsize)
ax[1][:grid](true)
ax[1][:axhline](y=.0, color = "k", linestyle="dashed")

#  CMB-S4 for ϕ est
ax[2][:set_title](L"\mu K arcmin = 1", fontsize=fontsize)
ax[2][:plot](rng_true, bias_B_2./rng_true, label = "quadratic delenser")
ax[2][:plot](rng_true, bias_L_2./rng_true, label = "likelihood approx")
ax[2][:set_xlabel](L"$r$", fontsize=fontsize)
ax[2][:tick_params](axis="y", left="off", right="off")
ax[2][:tick_params](axis="x", labelsize = fontsize-2)
ax[2][:axis]("tight")
ax[2][:grid](true)
ax[2][:axhline](y=.0, color = "k", linestyle="dashed")


#  CMB-S4 for ϕ est
ax[3][:set_title](L"\mu K arcmin = 0.5", fontsize=fontsize)
ax[3][:plot](rng_true, bias_B_3./rng_true, label = "quadratic delenser")
ax[3][:plot](rng_true, bias_L_3./rng_true, label = "likelihood approx")
ax[3][:set_xlabel](L"$r$", fontsize=fontsize)
ax[3][:tick_params](axis="y", labelsize = fontsize-2, labelright="on")
ax[3][:tick_params](axis="x", labelsize = fontsize-2)
ax[3][:axis]("tight")
ax[3][:grid](true)
ax[3][:axhline](y=.0, color = "k", linestyle="dashed")

#fig[:suptitle]("    Bias", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f2b.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end










#  .o88o.  o8o                                              .oooo.
#  888 `"  `"'                                            .dP""Y88b
# o888oo  oooo   .oooooooo oooo  oooo  oooo d8b  .ooooo.        ]8P'
#  888    `888  888' `88b  `888  `888  `888""8P d88' `88b     <88b.
#  888     888  888   888   888   888   888     888ooo888      `88b.
#  888     888  `88bod8P'   888   888   888     888    .o o.   .88P
# o888o   o888o `8oooooo.   `V88V"V8P' d888b    `Y8bod8P' `8bd88P'
#               d"     YD
#               "Y88888P'
#
#
#
# Explore higher order bias impact on the quadraitc delenser



using PyPlot, JLD
savefigs = true
savepath = "scripts/saved_run_0x3f77869c073f2346/"


rng_true = load(savepath*"out_figure3a.jld", "rng_true")

# low resolution survey for ϕ est
rest_B_1 = load(savepath*"out_figure3a.jld", "rest_B")
rest_1   = load(savepath*"out_figure3a.jld", "rest_modes")
rest_B_2 = load(savepath*"out_figure3b.jld", "rest_B")
rest_2   = load(savepath*"out_figure3b.jld", "rest_modes")

rest_B_3 = load(savepath*"out_figure3c.jld", "rest_B")
rest_3   = load(savepath*"out_figure3c.jld", "rest_modes")


# compute bias and standard error
bias_B_1  = mean((rest_B_1 .- rng_true),2)
bias_B_2  = mean((rest_B_2 .- rng_true),2)
bias_B_3  = mean((rest_B_3 .- rng_true),2)
bias_L_1  = mean((rest_1 .- rng_true),2)
bias_L_2  = mean((rest_2 .- rng_true),2)
bias_L_3  = mean((rest_3 .- rng_true),2)

standard_error_B_1 = sqrt.(mean((rest_B_1 .- rng_true).^2, 2))
standard_error_B_2 = sqrt.(mean((rest_B_2 .- rng_true).^2, 2))
standard_error_B_3 = sqrt.(mean((rest_B_3 .- rng_true).^2, 2))
standard_error_L_1 = sqrt.(mean((rest_1 .- rng_true).^2, 2))
standard_error_L_2 = sqrt.(mean((rest_2 .- rng_true).^2, 2))
standard_error_L_3 = sqrt.(mean((rest_3 .- rng_true).^2, 2))

maxylab = maximum(vcat(
    rng_true./standard_error_B_1,
    rng_true./standard_error_B_2,
    rng_true./standard_error_B_3,
    rng_true./standard_error_L_1,
    rng_true./standard_error_L_2,
    rng_true./standard_error_L_3
    ))

minbias, maxbias = vcat(
         bias_B_1./standard_error_B_1,
         bias_B_2./standard_error_B_2,
         bias_B_3./standard_error_B_3,
         bias_L_1./standard_error_L_1,
         bias_L_2./standard_error_L_2,
         bias_L_3./standard_error_L_3,
         ) |> x-> (minimum(x), maximum(x))


#' f2a.pdf : plot r / σ(r)
#' ================================================


fig, ax = plt[:subplots](1,1, figsize=(4.2,3.7), sharex=true,  sharey=true)
fontsize = 10

#  CMB-S3 for ϕ est
ax[:set_title](L"Quadratic delenser ($\mu K arcmin = 1$) ", fontsize=fontsize-1)
ax[:plot](rng_true, rng_true./standard_error_B_1, "C0", label = "Linear order (null test)")
ax[:plot](rng_true, rng_true./standard_error_B_2, "C3.", label = "All order E from E")
ax[:plot](rng_true, rng_true./standard_error_B_3, "C4", label = "All order E from B")
ax[:set_ylabel](L"$\frac{r}{\sigma(r)}$       ", fontsize=fontsize+4,  rotation="horizontal")
ax[:set_xlabel](L"$r$", fontsize=fontsize)
ax[:set_xticks]([.01, .03, .05, .07, .09])
ax[:tick_params](axis="y", labelsize = fontsize-2)
ax[:tick_params](axis="x", labelsize = fontsize-2)
ax[:legend](loc="best", fontsize=fontsize)
ax[:grid](true)


#fig[:suptitle]("    Detection", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f3a.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end





#' f2b.pdf : plot  Bias(r) / σ(r)
#' ================================================



fig, ax = plt[:subplots](1,1, figsize=(4.2,3.7), sharex=true,  sharey=true)
fontsize = 10

#  CMB-S3 for ϕ est
ax[:set_title](L"Quadratic delenser ($\mu K arcmin = 1$) ", fontsize=fontsize-1)
ax[:plot](rng_true, bias_B_1./rng_true, "C0", label = "Linear order (null test)")
ax[:plot](rng_true, bias_B_2./rng_true, "C3.", label = "All order E from E")
ax[:plot](rng_true, bias_B_3./rng_true, "C4", label = "All order E from B")
ax[:set_ylabel](L"$\frac{Bias(r)}{r}$        ", fontsize=fontsize+4,  rotation="horizontal")
ax[:set_xlabel](L"$r$", fontsize=fontsize)
ax[:set_xticks]([.01, .03, .05, .07, .09])
ax[:tick_params](axis="y", labelsize = fontsize-2)
ax[:tick_params](axis="x", labelsize = fontsize-2)
ax[:legend](loc="best", fontsize=fontsize)
ax[:grid](true)
ax[:axhline](y=.0, color = "k", linestyle="dashed")


#fig[:suptitle]("    Detection", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f3b.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end





#  .o88o.  o8o                                                       .o
#  888 `"  `"'                                                     .d88
# o888oo  oooo   .oooooooo oooo  oooo  oooo d8b  .ooooo.         .d'888
#  888    `888  888' `88b  `888  `888  `888""8P d88' `88b      .d'  888
#  888     888  888   888   888   888   888     888ooo888      88ooo888oo
#  888     888  `88bod8P'   888   888   888     888    .o           888
# o888o   o888o `8oooooo.   `V88V"V8P' d888b    `Y8bod8P'          o888o
#               d"     YD
#               "Y88888P'
#
#
#
#
#
#
#
#
#
#
#

using PyPlot, JLD
savefigs = true
savepath = "scripts/saved_run_0x3f77869c073f2346/"


rng_true = load(savepath*"out_figure4.jld", "rng_true")

# low resolution survey for ϕ est
rest_B_1 = load(savepath*"out_figure4.jld", "rest_B")
rest_1   = load(savepath*"out_figure4.jld", "rest_modes")


# compute bias and standard error
bias_B_1  = mean((rest_B_1 .- rng_true),2)
bias_L_1  = mean((rest_1 .- rng_true),2)

standard_error_B_1 = sqrt.(mean((rest_B_1 .- rng_true).^2, 2))
standard_error_L_1 = sqrt.(mean((rest_1 .- rng_true).^2, 2))

maxylab = maximum(vcat(
    rng_true./standard_error_B_1,
    rng_true./standard_error_L_1,
    ))

minbias, maxbias = vcat(
         bias_B_1./standard_error_B_1,
         bias_L_1./standard_error_L_1,
         ) |> x-> (minimum(x), maximum(x))


#' f4a.pdf : plot r / σ(r)
#' ================================================


fig, ax = plt[:subplots](1,1, figsize=(4.2,3.7), sharex=true,  sharey=true)
fontsize = 10

#  CMB-S3 for ϕ est
ax[:set_title](L"Nonstationary noise (baseline $\mu K arcmin = 1$)", fontsize=fontsize-1)
ax[:plot](rng_true, rng_true./standard_error_B_1, label = "quadratic delenser")
ax[:plot](rng_true, rng_true./standard_error_L_1, label = "likelihood approx")
ax[:set_ylabel](L"$\frac{r}{\sigma(r)}$       ", fontsize=fontsize+4,  rotation="horizontal")
ax[:set_xlabel](L"$r$", fontsize=fontsize)
ax[:set_xticks]([.01, .03, .05, .07, .09])
ax[:tick_params](axis="y", labelsize = fontsize-2)
ax[:tick_params](axis="x", labelsize = fontsize-2)
ax[:legend](loc="best", fontsize=fontsize)
ax[:grid](true)


#fig[:suptitle]("    Detection", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f4a.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end



#' f4b.pdf : plot  Bias(r) / σ(r)
#' ================================================

fig, ax = plt[:subplots](1,1, figsize=(4.2,3.7), sharex=true,  sharey=true)
fontsize = 10

#  CMB-S3 for ϕ est
ax[:set_title](L"Nonstationary noise (baseline $\mu K arcmin = 1$)", fontsize=fontsize-1)
ax[:plot](rng_true, bias_B_1./rng_true, label = "quadratic delenser")
ax[:plot](rng_true, bias_L_1./rng_true, label = "likelihood approx")
ax[:set_ylabel](L"$\frac{Bias(r)}{r}$        ", fontsize=fontsize+4,  rotation="horizontal")
ax[:set_xlabel](L"$r$", fontsize=fontsize)
ax[:set_xticks]([.01, .03, .05, .07, .09])
ax[:tick_params](axis="y", labelsize = fontsize-2)
ax[:tick_params](axis="x", labelsize = fontsize-2)
ax[:legend](loc="best", fontsize=fontsize)
ax[:grid](true)
ax[:axhline](y=.0, color = "k", linestyle="dashed")

#fig[:suptitle]("    Bias", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f4b.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end




#' f4c.pdf :
#' ================================================


fig = figure(figsize=(4.2,3.7))
fontsize = 10
σmod       = load(savepath*"out_figure4.jld", "σmod")
μKarcminT  = load(savepath*"out_figure4.jld", "μKarcminT")
#imshow(fftshift(σmod .* √2 .* μKarcminT))
imshow(fftshift(σmod))
colorbar()
title("Non-stationary noise modulation", fontsize=fontsize)
axis("tight")
axis("off")

#fig[:suptitle]("    Bias", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f4c.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end
