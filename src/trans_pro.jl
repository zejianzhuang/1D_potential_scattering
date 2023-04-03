

using CSV
using DataFrames



const me = 9.1094e-31 # Unit: kg
const L = 8e-10 # Unit: m
const hbar = 1.0546e-34 # Unit: Js
#const V0 = 20 # Unit: eV
const q = 1.602e-19; # Unit: J


function trans_pro!(w::Vector; V0=20.0)
    if length(w) == 2
        E = w[1] + w[2] * 1im # 
    else
        E = w[1] # Bound energy
    end

    eps = 1e-34
    kp = sqrt(2 * me * q * (E + V0) + 0im) / hbar
    k = sqrt(2 * me * q * E + 0im) / hbar # Unit: m^-1
    denoT = cos(kp * L) - 1im / 2.0 * (k / (kp + eps) + kp / (k + eps)) * sin(kp * L)

    T = 1 / (denoT + eps)
    return abs2(T)
end

w = 0:0.001:60
tp = [trans_pro!([ww]) for ww in w]

df = DataFrame(E = w, prob=tp)
CSV.write("output/trans_pro.csv", df)

println("💖💖💖success💖💖💖")