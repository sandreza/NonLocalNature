@info "diffusivity kernel"
using FourierCore, FourierCore.Grid, FourierCore.Domain
using FFTW, LinearAlgebra, BenchmarkTools, Random, JLD2
# using GLMakie, HDF5
using ProgressBars
rng = MersenneTwister(1234)
Random.seed!(123456789)

#=
scaleit = 2^3
tstart = 2^5 * scaleit
tend = 2^6 * scaleit


load_psi!(ψ)
P * ψ;
ζ .= Δ .* ψ
P⁻¹ * ζ
are defined outside
=#

maxind = minimum([40, floor(Int, N[1] / 2)])
index_choices = 2:maxind

start_index = floor(Int, tstart / Δt)

sθ .= 0.0
θ .= 0

for index_choice in ProgressBar(index_choices)
    kᶠ = kˣ[index_choice]
    @. θ += cos(kᶠ * x) / (kᶠ)^2 / κ # scaling so that source is order 1
end
P * θ # in place fft
@. 𝒟θ = 𝒟κ * θ
P⁻¹ * 𝒟θ # in place fft
P⁻¹ * θ # in place fft
sθ .+= -𝒟θ # add to source

t = [0.0]
iend = ceil(Int, tend / Δt)

# new realization of flow
rand!(rng, φ) # between 0, 1
φ .*= 2π # to make it a random phase

θ̄ = arraytype(zeros(ComplexF64, N, N, N_ens))

iter = ProgressBar(1:iend)
ke_list = Float64[]
t .= 0.0
for i = iter
    step!(S, S̃, φ, φ̇, k₁, k₂, k₃, k₄, Δt, rng, t, parameters)
    if i % mod_index == 0
        θ_min, θ_max = extrema(real.(θ))
        ζ_min, ζ_max = extrema(real.(ζ))
        set_multiline_postfix(iter, "θ_min: $θ_min \nθ_max: $θ_max \nζ_min: $ζ_min \nζ_max: $ζ_max")
    end
    if i > start_index
        θ̄ .+= Δt .* θ
    end
    if i % mod_index == 0
        push!(ke_list, real(mean(u .* u + v .* v)))
    end
end

θ̄ ./= (tend - tstart)
θ̄_A = Array(real.(θ̄))

tmpsave = Array(real.(mean(θ̄, dims=2)))[:, 1, :]
tmp = Array(real.(fft(mean(θ̄, dims=(2, 3))[:]))) # tmp = real.(fft(Array(mean(θ[:,:,1:10], dims = (2,3)))[:]))
kxa = Array(kˣ)[:]
effective_diffusivities = ((Ns[1] / 2) ./ tmp) ./ (kxa .^ 2) .- κ # (((N[1] / 2) ./ tmp) .- λ) ./ (kxa .^ 2) .- κ
effective_diffusivities = effective_diffusivities[index_choices]

# estimate kernel on grid
kernel = real.(fft([0.0, effective_diffusivities..., zeros(65)..., reverse(effective_diffusivities)...]))
kernel = kernel .- mean(kernel[63:65])
kernel = circshift(kernel, 64)

@info "done with diffusivity kernel"
#=
fig = Figure()
ax = Axis(fig[1, 1], xlabel="x", ylabel="y", aspect=1)
scatter!(ax, kernel)
display(fig)
=#
#=
list1 = Float64[]
for i in 1:128
    tmp = Array(real.(fft(mean(θ̄[:,:,i], dims=2)[:]))) # tmp = real.(fft(Array(mean(θ[:,:,1:10], dims = (2,3)))[:]))
    kxa = Array(kˣ)[:]
    eff = (((N[1] / 2) ./ tmp) .- λ) ./ (kxa .^ 2) .- κ
    println(eff[2])
    push!(list1, eff[2])
end
=#
