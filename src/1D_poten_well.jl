

using CSV
using DataFrames
using NLsolve


const me = 9.1094e-31 # Unit: kg
const L = 8e-10 # Unit: m
const hbar = 1.0546e-34 # Unit: Js
#const V0 = 20 # Unit: eV
const q = 1.602e-19; # Unit: J

k(E) = sqrt(2*0.510e6*E + 0im)

function deno_of_T!(w::Vector; V0=20.0)
    if length(w) == 2
        E = w[1] + w[2] * 1im # 
    else
        E = w[1] # Bound energy
    end

    eps = 1e-34
    kp = sqrt(2 * me * q * (E + V0) + 0im) / hbar
    k = sqrt(2 * me * q * E + 0im) / hbar # Unit: m^-1
    denoT = cos(kp * L) - 1im / 2.0 * (k / (kp + eps) + kp / (k + eps)) * sin(kp * L)

    if abs(real(denoT)) <= 1e-8 && abs(imag(denoT)) <= 1e-8
        return [0.0, 0.0]
    elseif abs(imag(denoT)) <= 1e-8
        return [real(denoT), 0.0]
    else
        return [real(denoT), imag(denoT)]
    end
end


function T!(w::Vector; V0=20)
    eps = 1e-8
    denoT = deno_of_T!(w, V0=V0)

    return 1 / (denoT[1] + denoT[2] * 1im + eps)
end

# Generate the data to plot. The diagram will be done by matplotlib package.
# From the plot, we can figure out some positions of the poles to bound states.
E = -19.9:0.001:-0.1
df1 = DataFrame(
    reE = E,
    real_denoT = [deno_of_T!([w])[1] for w in E],
                imag_denoT = [deno_of_T!([w])[2] for w in E] )
# Two-dims data to plot. This shows the positions of the pole structures to resonances.
# The figure will be plotted by Mathematica.
magn_T = []
reE = []
imE = []
dim = 100
for rew in range(0.01, 70, dim)
    append!(reE, ones(dim)*rew)
    append!(imE, [imw for imw in range(-15, 0, dim)])
    append!(magn_T, [abs(T!([rew, imw])) for imw in range(-15, 0, dim)] )
end
df2 = DataFrame(real_of_E = reE, imag_of_E = imE, magnitude_of_T = magn_T)


# Store the data to a csv file
CSV.write("./output/denominator_of_T_with_E_less_than_0.csv", df1)
CSV.write("./output/magnitude_of_T.csv", df2)

# Find the poles to bound states
init_w = [[-3.0011, 0.], [-8.1100, 0.], [-12.00, 0.], [-15.1, 0.], [-17.88777, 0.], [-19.9, 0.]]

tex_pole_to_BS = Vector{Union{Missing, String}}(missing, 7)
pole_to_BS = tex_pole_to_BS
BS_in_k = Vector{Union{Missing, String}}(missing, 7)
converge = []

for (i, w) in enumerate(init_w)
    r = nlsolve(deno_of_T!, w)
    if converged(r) == true
        BSk = k(r.zero[1])
        append!(converge, ["true"])
        tex_pole_to_BS[i] = "\$$(round(r.zero[1], digits=3))\$"
        pole_to_BS[i] = "$(round(r.zero[1], digits=3) )"
        BS_in_k[i] = "$(imag(BSk))j"
        #append!(pole_to_BS, ["\$$(round(r.zero[1], digits=3))\$"] )
    end
end

# Find the poles to resonaces
init_cw = [[0.633333, -0.9],
    [8.0, -3.0],
    [16.8, -4.9],
    [26.0, -6.0],
    [37.0, -8.0],
    [49.0, -10.0],
    [63.0, -12.0]]
tex_pole_to_R = Vector{String}()
pole_to_resonance = Vector{String}()
R_in_k = Vector{String}()

for w in init_cw
    r = nlsolve(deno_of_T!, w)
    Rk = k(r.zero[1]+r.zero[2]*1im)
    append!(tex_pole_to_R, ["\$$(round(r.zero[1], digits=3))$(+round(r.zero[2], digits=3) )i\$"])
    append!(pole_to_resonance, ["$(r.zero[1])$(+r.zero[2])j"] )
    append!(R_in_k, ["$(real(Rk))$(+imag(Rk))j"] )
    #println("$(round(r.zero[1], digits=3)) $(round(r.zero[2], digits=3) )i")
end


# Create a csv file used for generating tex table code.
tex_pole = DataFrame()
tex_pole[!, "Pole to bound state"] = tex_pole_to_BS
tex_pole[!, "Pole to resonace"] = tex_pole_to_R

pole = DataFrame()
pole[!, :pole_to_BS] = pole_to_BS
pole[!, :pole_to_R] = pole_to_resonance
pole[!, :BS_in_k] = BS_in_k
pole[!, :R_in_k] = R_in_k

CSV.write("./_to_latex/pole_csv_to_tex.csv", tex_pole)
CSV.write("./output/pole_structure.csv", pole)

println("💖💖💖success💖💖💖")
