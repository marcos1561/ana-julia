using Plots

include("anajulia.jl")
import .anajulia

# Configurações da corda
dx = 0.01
dt = (300 * 1/dx)^-1
k = 1
l = Int(1/dx)

# Configurações inicias (gaussian)
center = 0.8
sigma = 0.01
height = 1

# Integração
num_steps = 80


x = Vector{Float32}(0:l-1)
y = anajulia.init_state.gaussian_package(
    x, center, sigma)
# y = anajulia.init_state.pluck_string(
#     x, center, 1)

wave = anajulia.Wave(y)

plot(wave.y, label="inicial")
for i in 1:num_steps
    anajulia.wave_step!(wave, k)
end
plot!(wave.y, label="final")