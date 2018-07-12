include("laser.jl")
using Laser
# using Plots
using JLD
using LaTeXStrings

function averagex(H::Array{Float64,2}, nx::Int)
    for j in 1:nx
        for i in 3:nx - 2
            H[i, j] = (H[i - 2, j] + H[i - 1, j] + H[i + 1, j] + H[i + 2, j]) / 4.0
        end
    end
    H
end

function averagey(H::Array{Float64,2}, nx::Int)
    for j in 3:nx - 2
        for i in 1:nx
            H[i, j] = (H[i, j - 2] + H[i, j - 1] + H[i, j + 1] + H[i, j + 2]) / 4.0
        end
    end
    H
end


function sigma(Ein::Float64, paras::Array{Float64,1})
    E = Ein * Esi
    E0 = paras[1]
    sigma0 = paras[2]
    ya = paras[3]
    P = paras[4]
    yw = paras[5]
    y0 = paras[6]
    y1 = paras[7]

    x = E / E0 - y0
    y = sqrt(x^2 + y1^2)
    F = ((x - 1)^2 + yw^2) * y^(0.5 * P - 5.5) * (1 + sqrt(y / ya))^(-P)

    return sigma0 * F
end

function getFxFy(FxFy::Array{Complex128,2}, ele1::Array{Float64,1}, ele2::Array{Float64,1}, ele3::Array{Float64,1}, t::Array{Float64,1}, coeff::Float64, coeff2::Float64, nt::Int)
    # const coeff2 = 1.0
    for j in 1:nt, i in 1:nt
        # FxFy[i, j] = ele1[i] * ele1[j] * (t[i] < t[j]) * cis(-coeff * (t[j] * (ele2[j])^2 - t[i] * ele2[i]^2))
        # FxFy[i, j] = ele1[i] * ele1[j] * (t[i] < t[j]) * cis(-coeff2 * (t[j] * (ele3[j])^2 - t[i] * ele3[i]^2))
        FxFy[i, j] = ele1[i] * ele1[j] * (t[i] < t[j]) * cis(-coeff * (t[j] * (ele2[j])^2 - t[i] * ele2[i]^2) - coeff2 * (t[j] * (ele3[j])^2 - t[i] * ele3[i]^2))
    end
end 

println("nt =", nt)

paras1 = [1.361e1, 9.492e2, 1.469e0, 3.188e0, 2.039e0, 4.434e-1, 2.136e0]
paras2 = [1.720e0, 1.369e4, 3.288e1, 2.963e0, 0e0, 0e0, 0e0]

dw = 2 * pi / nt / dt
t = collect((1:nt) * dt)
###  4τ --- 5τ  --- 0.18 to 0.23
###  4τ ---- 4.5 τ ------ 0.4
FxFy = zeros(Complex128, nt, nt)
## 0.04 可以
const coeff = 0.03
const coeff2 = coeff * 3#2.7

@time getFxFy(FxFy, ele1, ele2, ele3, t, coeff, coeff2, nt)

p_fft! = plan_fft!(similar(FxFy))#, flags=FFTW.MEASURE)

@time p_fft! \ FxFy

const np = 1000
dp = 2.0 / np
#  energy = [1:np]*dp
p = (-np / 2:np / 2 - 1) * dp
const Esi = 27.2114
Ei = -79.0 / Esi
Eg = -54.4 / Esi

spectrum = zeros(np, np)
σ = zeros(4)
Δ = coeff * sum(ele2.^2) * dt / (t[end] - t[1])

@time for n2 in 1:np
    for n1 in 1:np
        E1 = (p[n1])^2 / 2
        E2 = (p[n2])^2 / 2

        Ea = E1 + Eg
        Eb = E2 + Eg
        Ef = E1 + E2

        σ[1] = sigma(Ea - Ei, paras1)
        σ[2] = sigma(Ef - Ea, paras2)
        σ[3] = sigma(Eb - Ei, paras1)
        σ[4] = sigma(Ef - Eb, paras2)
        Ka = FxFy[round(Int, (Ea - Ei) / dw + 1), round(Int, (Ef - Ea) / dw + 1)]
        Kb = FxFy[round(Int, (Eb - Ei) / dw + 1), round(Int, (Ef - Eb) / dw + 1)]

        spectrum[n1, n2] = abs(p[n1] * p[n2]) * abs(sqrt(σ[1] * σ[2] / (Ea - Ei) / (Ef - Ea)) * Ka + sqrt(σ[3] * σ[4] / (Eb - Ei) / (Ef - Eb)) * Kb)^2

        # spectrum[n1, n2] = abs(p[n1]*p[n2])*abs(sqrt(σ[1]*σ[2]/(Ea - Ei + Δ)/(Ef - Ea))*Ka + sqrt(σ[3]*σ[4]/(Ea - Ei + Δ)/(Ef - Eb))*Kb)^2

        # spectrum[n2, n1] = spectrum[n1, n2]
    end
end

spectrum ./= maximum(spectrum)
for i in 1:5
    spectrum = averagex(spectrum, np)
    spectrum = averagey(spectrum, np)
end
spectrum ./= maximum(spectrum)

save("data/sp.jld", "sp", spectrum)
