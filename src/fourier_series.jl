using Plots

include("anajulia.jl")
import .anajulia

# Configurações da corda
dx = 0.01
dt = (300 * 1/dx)^-1
k = 1
l = Int(1/dx)

# Configurações inicias (gaussian)
center = 0.1
sigma = 0.01


# Condifurações do ponto de coleta e série
point_pos = 0.8
num_freqs = 15

x = Vector{Float32}(0:l-1)
y = anajulia.init_state.gaussian_package(
    x, center, sigma)

wave = anajulia.Wave(y)

time, waveform = anajulia.watch_point(wave, point_pos, k, dt, 2)

speed = anajulia.get_speed(k, dt)

@gif for num_f in 1:num_freqs
    as, bs, freqs = anajulia.fourier_coefficients(waveform, time, speed, wave.l, num_f)
    y_forier = anajulia.forier_sum(time, as, bs, freqs)
    
    # plot(time, waveform, label="Sinal")
    # plot!(time, y_forier, label="Série de fourier ($(num_freqs) frequências)")
    plot(time, waveform)
    ylims!(-0.2, 0.5)
    plot!(time, y_forier)
end fps=5
