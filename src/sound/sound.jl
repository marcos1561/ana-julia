using Plots

include("../anajulia.jl")
import .anajulia

# Configurações da corda
dx = 0.01
dt = (300 * 1/dx)^-1
k = 1
l = Int(1/dx)

# Configurações inicias (gaussian)
center = 0.5
sigma = 0.01
height = 1

# Condifurações do ponto de coleta e série
point_pos = 0.05
num_freqs = 15

x = Vector{Float32}(0:l-1)
# y = anajulia.init_state.pluck_string(
#     x, center, height)
y = anajulia.init_state.gaussian_package(
    x, center, sigma)

wave = anajulia.Wave(y)
p_init_state = plot(wave.y, label="Estado inicial")
    
time, waveform = anajulia.watch_point(wave, point_pos, k, dt, 10)

using CSV, Tables
# CSV.write("data/waveform.csv", Tables.table(waveform))

as_t, bs_t, freqs_t = anajulia.fourier_transform(waveform, time, 3000, 1)
power_t = anajulia.get_power(as_t, bs_t)

speed = anajulia.get_speed(k, dt)
as, bs, freqs = anajulia.fourier_coefficients(waveform, time, speed, wave.l, 20)
power = anajulia.get_power(as, bs)

p_power = plot(freqs_t, power_t/maximum(power_t), legend=false)
xlabel!(p_power, "Frequência (Hz)")
ylabel!(p_power, "Potência")
scatter!(p_power, freqs, power/maximum(power))

plot(p_init_state, p_power, layout=(2,1))
