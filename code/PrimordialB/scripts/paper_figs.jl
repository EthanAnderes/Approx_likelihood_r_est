using PyPlot, JLD
savefigs = true
savepath = "scripts/saved_run_0x3f77869c073f2346/"


#figure 1
#' Plot
fig, ax = plt[:subplots](1,1, figsize=(6,3.5))
fontsize = 20
x_grid = [x for x in range(-10,21) for y in range(0,21)]
y_grid = [y for x in range(-10,21) for y in range(0,21)]
ax[:scatter](x_grid, y_grid, s=1.3, color="k")
Θ = linspace(0., π, 100)
r0= 4.3
r1= sqrt(2)*r0
r2= sqrt(3)*r0
r3= sqrt(4)*r0
ax[:plot](r0*cos(Θ), r0*sin(Θ))
ax[:plot](r1*cos(Θ), r1*sin(Θ))
ax[:plot](r2*cos(Θ), r2*sin(Θ))
ax[:plot](r3*cos(Θ), r3*sin(Θ))
ax[:text](0, -1, L"$\ell_1$",  horizontalalignment="center", fontsize=fontsize)
ax[:text](-12, 5, L"$\ell_2$",  horizontalalignment="center", rotation=90, fontsize=fontsize)
ax[:text](0, 1, L"$^0B^{\rm del}_{\bf l}$",  horizontalalignment="center", fontsize=fontsize)
ax[:text](0, 4.5, L"$^1B^{\rm del}_{\bf l}$",  horizontalalignment="center", fontsize=fontsize-1)
ax[:text](0, 6.2, L"$^2B^{\rm del}_{\bf l}$",  horizontalalignment="center", fontsize=fontsize-2)
ax[:text](0, 7.6, L"$^3B^{\rm del}_{\bf l}$",  horizontalalignment="center", fontsize=fontsize-3)
ax[:set_xlim]([-11, 11])
ax[:set_ylim]([0, 11])
ax[:set_aspect]("equal")
ax[:axis]("off")

fig[:tight_layout]()

#' Now save
savename  = "f1.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end



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


#figure 2

fig, ax = plt[:subplots](1,1, figsize=(4.5,3.5))
fontsize = 12
ax[:plot](kbins, cϕϕ, lw=2, label = L"$\ell^4 C_\ell^{\phi\phi}$")
ax[:plot](kbins, N0ℓ_a, lw=2, label = L"$\Delta_{\rm T} = 10 \ \mu{\rm K-arcmin}$")
ax[:plot](kbins, N0ℓ_b, lw=2, label = L"$\Delta_{\rm T} = 1 \ \mu{\rm K-arcmin}$")
ax[:plot](kbins, N0ℓ_c, lw=2, label = L"$\Delta_{\rm T} = 0.5 \ \mu{\rm K-arcmin}$")
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:set_xlim]([10.,1.2e3])
ax[:set_ylim]([1.e-8,4e-6])
ax[:set_title](L"$\phi$ uncertainty", fontsize=fontsize)
ax[:set_xlabel](L"multiple $\ell$", fontsize=fontsize)
ax[:legend](loc="best", fontsize=fontsize-2)
#ax[:set_ylim]([1e-8, 1.5e-6])
ax[:grid](true)


#' Now save
savename  = "f2.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end


#figure 3

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

#' f3.pdf : plot r / σ(r) and Bias(r)/r
#' ================================================


fig, ax = plt[:subplots](2,3, figsize=(10,6), sharex=true)
fontsize = 13

#  CMB-S3 for ϕ est
ax[1,1][:set_title](L"$\Delta_{\rm T} = 10 \ \mu{\rm K-arcmin}$", fontsize=fontsize)
ax[1,1][:plot](rng_true, rng_true./standard_error_L_1, "b", lw = 2,label = L"${\rm likelihood\ approx}$")
ax[1,1][:plot](rng_true, rng_true./standard_error_B_1, "r.", lw = 2,label = L"${\rm quadratic\ delenser}$")
ax[1,1][:set_ylabel](L"$\frac{r}{\sigma(r)}$       ", fontsize=fontsize+4, rotation="horizontal")
#ax[1,1][:set_xlabel](L"$r$", fontsize=fontsize)
ax[1,1][:set_ylim]([0,17])
ax[1,1][:set_xticks]([.01, .03, .05, .07, .09])
ax[1,1][:tick_params](axis="y", labelsize = fontsize-2)
ax[1,1][:tick_params](axis="x", labelsize = fontsize-2)
ax[1,1][:legend](loc="best", fontsize=fontsize)
ax[1,1][:grid](true)

#  CMB-S4 for ϕ est
ax[1,2][:set_title](L"$\Delta_{\rm T} = 1 \ \mu{\rm K-arcmin}$", fontsize=fontsize)
ax[1,2][:plot](rng_true, rng_true./standard_error_B_2, "r.",lw = 2,label = L"${\rm quadratic\ delenser}$")
ax[1,2][:plot](rng_true, rng_true./standard_error_L_2, "b-",lw = 2,label = "likelihood approx")
#ax[1,2][:set_xlabel](L"$r$", fontsize=fontsize)
ax[1,2][:set_ylim]([0,17])
ax[1,2][:tick_params](axis="y", left="off", right="off")#labelsize = fontsize-2, pad = 5, labelright="on")
ax[1,2][:tick_params](axis="x", labelsize = fontsize-2)
ax[1,2][:grid](true)


#  CMB-S4 for ϕ est
ax[1,3][:set_title](L"$\Delta_{\rm T} = 0.5 \ \mu{\rm K-arcmin}$", fontsize=fontsize)
ax[1,3][:plot](rng_true, rng_true./standard_error_B_3, "r.",lw = 2,label = "quadratic delenser")
ax[1,3][:plot](rng_true, rng_true./standard_error_L_3, "b-",lw = 2,label = "likelihood approx")
#ax[1,3][:set_xlabel](L"$r$", fontsize=fontsize)
ax[1,3][:set_ylim]([0,17])
ax[1,3][:tick_params](axis="y", labelsize = fontsize-2, labelright="on")
ax[1,3][:tick_params](axis="x", labelsize = fontsize-2, pad = 5)
ax[1,3][:axis]("tight")
ax[1,3][:grid](true)


#  CMB-S3 for ϕ est
#ax[2,1][:set_title](L"$\Delta_{\rm T} = 10 \ \mu{\rm K-arcmin}$", fontsize=fontsize)
ax[2,1][:plot](rng_true, bias_L_1./rng_true, "b-", lw =2, label = L"${\rm likelihood\ approx}$")
ax[2,1][:plot](rng_true, bias_B_1./rng_true, "r.", lw = 2, label = L"${\rm quadratic\ delenser}$")
ax[2,1][:set_ylabel](L"$\frac{{\rm Bias}(r)}{r}$        ", fontsize=fontsize+4, rotation="horizontal")
ax[2,1][:set_xlabel](L"$r$", fontsize=fontsize)
ax[2,1][:set_ylim]([-0.02,0.2])
ax[2,1][:set_xticks]([.01, .03, .05, .07, .09])
ax[2,1][:tick_params](axis="y", labelsize = fontsize-2)
ax[2,1][:tick_params](axis="x", labelsize = fontsize-2)
ax[2,1][:legend](loc="best", fontsize=fontsize)
ax[2,1][:grid](true)
ax[2,1][:axhline](y=.0, "k.")

#  CMB-S4 for ϕ est
#ax[2,2][:set_title](L"$\Delta_{\rm T} = 1 \ \mu{\rm K-arcmin}$", fontsize=fontsize)
ax[2,2][:plot](rng_true, bias_B_2./rng_true, "r.", lw = 2,label = "quadratic delenser")
ax[2,2][:plot](rng_true, bias_L_2./rng_true, "b-", lw = 2,label = "likelihood approx")
ax[2,2][:set_xlabel](L"$r$", fontsize=fontsize)
ax[2,2][:set_ylim]([-0.02,0.2])
ax[2,2][:tick_params](axis="y", left="off", right="off")
ax[2,2][:tick_params](axis="x", labelsize = fontsize-2)
ax[2,2][:grid](true)


#  CMB-S4 for ϕ est
#ax[2,3][:set_title](L"$\Delta_{\rm T} = 0.5 \ \mu{\rm K-arcmin}$", fontsize=fontsize)
ax[2,3][:plot](rng_true, bias_B_3./rng_true, "r.", lw = 2,label = "quadratic delenser")
ax[2,3][:plot](rng_true, bias_L_3./rng_true, "b-", lw = 2,label = "likelihood approx")
ax[2,3][:set_xlabel](L"$r$", fontsize=fontsize)
ax[2,3][:set_ylim]([-0.02,0.2])
ax[2,3][:tick_params](axis="y", labelsize = fontsize-2, labelright="on")
ax[2,3][:tick_params](axis="x", labelsize = fontsize-2)
#ax[2,3][:axis]("tight")
ax[2,3][:grid](true)



#fig[:suptitle]("    Detection", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f3.pdf"
if savefigs
      fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end



#for figure 4 use
dtec_level_full = rng_true./standard_error_B_2
bias_level_full = bias_B_2./rng_true



# figure 4: Explore higher order bias impact on the quadraitc delenser


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


#' f4a.pdf : plot r / σ(r)
#' ================================================


fig, ax = plt[:subplots](1,1, figsize=(4.2,3.7), sharex=true,  sharey=true)
fontsize = 13

#  CMB-S3 for ϕ est
ax[:set_title](L"$\Delta_{\rm T} = 1 \ \mu{\rm K-arcmin}$", fontsize=fontsize-1)
ax[:plot](rng_true, rng_true./standard_error_B_1, "k-", lw=2, label = L"${\rm Linear\ order\ (null\ test)}$")
ax[:plot](rng_true, rng_true./standard_error_B_2, "g:", lw=6, label = L"${\rm All\ order\ E\ from\ E}$")
ax[:plot](rng_true, rng_true./standard_error_B_3, "b--", lw=2, label = L"${\rm All\ order\ B\ from\ E}$")
ax[:plot](rng_true, dtec_level_full, "r.", lw=2, label = L"${\rm All\ order\ E/B\ from\ E}$")
ax[:set_ylabel](L"$\frac{r}{\sigma(r)}$       ", fontsize=fontsize+4,  rotation="horizontal")
ax[:set_xlabel](L"$r$", fontsize=fontsize)
ax[:set_xticks]([.01, .03, .05, .07, .09])
ax[:tick_params](axis="y", labelsize = fontsize-2)
ax[:tick_params](axis="x", labelsize = fontsize-2)
#ax[:legend](loc="best", fontsize=fontsize)
ax[:grid](true)
ax[:set_xlim]([0.,0.1])

#fig[:suptitle]("    Detection", fontsize=fontsize+3)
fig[:tight_layout]()


#' Now save
savename  = "f4a.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end


fig, ax = plt[:subplots](1,1, figsize=(4.5,3.7), sharex=true,  sharey=true)
fontsize = 13

#  CMB-S3 for ϕ est
ax[:set_title](L"$\Delta_{\rm T} = 1 \ \mu{\rm K-arcmin}$", fontsize=fontsize-1)
ax[:plot](rng_true, bias_B_1./rng_true, "k-", lw=2, label = L"${\rm Linear\ order\ (null\ test)}$")
ax[:plot](rng_true, bias_B_2./rng_true, "g:", lw=6, label = L"${\rm All\ order\ E\ from\ E}$")
ax[:plot](rng_true, bias_B_3./rng_true, "b--", lw=2, label = L"${\rm All\ order\ B\ from\ E}$")
ax[:plot](rng_true, bias_level_full, "r.", lw=2, label = L"${\rm All\ order\ E/B\ from\ E}$")
ax[:set_ylabel](L"$\frac{{\rm Bias}(r)}{r}$        ", fontsize=fontsize+4,  rotation="horizontal")
ax[:set_xlabel](L"$r$", fontsize=fontsize)
ax[:set_xticks]([.01, .03, .05, .07, .09])
ax[:tick_params](axis="y", labelsize = fontsize-2)
ax[:tick_params](axis="x", labelsize = fontsize-2)
ax[:legend](loc="best", fontsize=fontsize-2)
ax[:grid](true)
ax[:set_xlim]([0.,0.1])
#ax[:axhline](y=.0, color = "k", linestyle="dashed")


#fig[:suptitle]("    Detection", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f4b.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end




#figure 5

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




#' f5a.pdf :
#' ================================================


fig = figure(figsize=(4.5,3.7))
fontsize = 13
σmod       = load(savepath*"out_figure4.jld", "σmod")
μKarcminT  = load(savepath*"out_figure4.jld", "μKarcminT")
#imshow(fftshift(σmod .* √2 .* μKarcminT))
imshow(fftshift(σmod))
colorbar()
#title("Non-stationary noise modulation", fontsize=fontsize)
text(0., 220., "Non-stationary noise modulation", fontsize=fontsize)
axis("tight")
axis("off")
#fig[:suptitle]("    Bias", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f5a.pdf"

if savefigs
  fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end

#' f5b.pdf : plot r / σ(r)
#' ================================================


fig, ax = plt[:subplots](1,1, figsize=(4.2,3.7), sharex=true,  sharey=true)
fontsize = 13

#  CMB-S3 for ϕ est
#ax[:set_title](L"Nonstationary noise (baseline $\mu K arcmin = 1$)", fontsize=fontsize-1)
ax[:plot](rng_true, rng_true./standard_error_L_1, "b-", lw=2, label = L"${\rm likelihood\ approx}$")
ax[:plot](rng_true, rng_true./standard_error_B_1, "r.", lw=2, label = L"${\rm quadratic\ delenser}$")
ax[:set_ylabel](L"$\frac{r}{\sigma(r)}$       ", fontsize=fontsize+4,  rotation="horizontal")
ax[:set_xlabel](L"$r$", fontsize=fontsize)
ax[:set_xticks]([.01, .03, .05, .07, .09])
ax[:tick_params](axis="y", labelsize = fontsize-2)
ax[:tick_params](axis="x", labelsize = fontsize-2)
ax[:legend](loc="best", fontsize=fontsize)
ax[:grid](true)
ax[:set_xlim]([0.,0.1])


#fig[:suptitle]("    Detection", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f5b.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end



#' f5c.pdf : plot  Bias(r) / σ(r)
#' ================================================

fig, ax = plt[:subplots](1,1, figsize=(4.5,3.7), sharex=true,  sharey=true)
fontsize = 13

#  CMB-S3 for ϕ est
#ax[:set_title](L"Nonstationary noise (baseline $\mu K arcmin = 1$)", fontsize=fontsize-1)
ax[:plot](rng_true, bias_L_1./rng_true, "b-", lw=2, label = L"${\rm likelihood\ approx}$")
ax[:plot](rng_true, bias_B_1./rng_true, "r.", lw=2, label =L"${\rm quadratic\ delenser}$")
ax[:set_ylabel](L"$\frac{{\rm Bias}(r)}{r}$        ", fontsize=fontsize+4,  rotation="horizontal")
ax[:set_xlabel](L"$r$", fontsize=fontsize)
ax[:set_xticks]([.01, .03, .05, .07, .09])
ax[:tick_params](axis="y", labelsize = fontsize-2)
ax[:tick_params](axis="x", labelsize = fontsize-2)
ax[:legend](loc="best", fontsize=fontsize)
ax[:grid](true)
ax[:set_xlim]([0.,0.1])
#ax[:axhline](y=.0, color = "k", linestyle="dashed")

#fig[:suptitle]("    Bias", fontsize=fontsize+3)
fig[:tight_layout]()

#' Now save
savename  = "f5c.pdf"
if savefigs
     fig[:savefig](joinpath(savepath, savename),bbox_inches="tight")
end
