using Plots

include("../anajulia.jl")
import .anajulia

# Configurações da corda
int_pars = anajulia.IntPars(
    k  = 1,
    c  = 500,
    dx = 0.01,
    b  = 0.1,
    epsilon = 1e-5,
    n = Int(1/dx),
)

# dt = anajulia.get_dt(dx, speed)
# dt = (300 * 1/dx)^-1
# k = 1

# Configurações inicias (gaussian)
center = 0.2
sigma = 0.01
height = 1

# Condifurações do ponto de coleta e série
point_pos = 0.05
time = 1

x = Vector{Float32}(0:l-1)
y = anajulia.init_state.pluck_string(
    x, center, height)
# y = anajulia.init_state.gaussian_package(
#     x, center, sigma)
    
wave = anajulia.Wave(y)
p_init_state = plot(wave.y, label="Estado inicial")

period = (wave.l * dx) / speed
num_periods = trunc(Int, time/period)
time, waveform = anajulia.watch_point(wave, point_pos, k, dt, num_periods)

using WAV
fs = trunc(Int, 1/dt)

p = plot(time[1:500], waveform[1:500])
display(p)
wavplay(waveform, fs)

# using CSV, Tables
# CSV.write("data/waveform.csv", Tables.table(waveform))

# as_t, bs_t, freqs_t = anajulia.fourier_transform(waveform, time, 3000, 1)
# power_t = anajulia.get_power(as_t, bs_t)

# # speed = anajulia.get_speed(k, dt)
# # as, bs, freqs = anajulia.fourier_coefficients(waveform, time, speed, wave.l, 20)
# # power = anajulia.get_power(as, bs)

# p_power = plot(freqs_t, power_t/maximum(power_t), legend=false)
# xlabel!(p_power, "Frequência (Hz)")
# ylabel!(p_power, "Potência")
# # scatter!(p_power, freqs, power/maximum(power))

# plot(p_init_state, p_power, layout=(2,1))
