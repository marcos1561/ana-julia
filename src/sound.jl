using Plots

include("anajulia.jl")
import .anajulia

# Configurações da corda
int_pars = anajulia.IntPars(
    k  = 1/4,
    c  = 880,
    dx = 0.01,
    b  = 5,
    epsilon = 1e-8*1,
    l = 1,
)

# Configurações inicias (gaussian)
center = 0.2
sigma = 0.01
height = 1

# Condifurações do ponto de coleta e série
point_pos = 0.05
time = 1

# Estado inicial
x = Vector{Float32}(0:int_pars.n-1)
y = anajulia.init_state.pluck_string(
    x, center, height)
# y = anajulia.init_state.gaussian_package(
#     x, center, sigma)

wave = anajulia.Wave(y)
p_init_state = plot(wave.y, label="Estado inicial")

# Colentado o sinal
period = (int_pars.n * int_pars.dx) / int_pars.c
num_periods = trunc(Int, time/period)
time, waveform = anajulia.watch_point(wave, point_pos, int_pars, num_periods)

print(time[end])

# Áudio
using WAV
fs = trunc(Int, 1/int_pars.dt)
wavplay(waveform, fs)
# wavwrite(waveform, fs, "sound/la_note_0_2.wav")

# ==
# Grapfico do Sinal
# ==
# p = plot(time[1:500], waveform[1:500])
# display(p)

# ==
# Espectro
# ==
mask = time .< time[end]*0.2
as_t, bs_t, freqs_t = anajulia.fourier_transform(waveform[mask], time[mask], 3000, 5)
power_t = anajulia.get_power(as_t, bs_t)

# ==
# Gráfico do Espectro
# ==
p_power = plot(freqs_t, power_t/maximum(power_t), legend=false)
xlabel!(p_power, "Frequência (Hz)")
ylabel!(p_power, "Potência (Un. arbitrária)")

# Harmônicos
harmonic = int_pars.c / ( 2 * int_pars.n * int_pars.dx)
num_harms = 5
for n in 1:num_harms
    freq = harmonic * n 
    plot!(p_power, [freq, freq], [0, 0.5], linestyle=:dash, c="black")
end
savefig("animacoes/la_power_high_stiff_ana_julia.png")
display(p_power)

# ==
# Série de Fourier com os harmônicos
# ==
# speed = anajulia.get_speed(k, dt)
# as, bs, freqs = anajulia.fourier_coefficients(waveform, time, speed, wave.l, 20)
# power = anajulia.get_power(as, bs)
# # scatter!(p_power, freqs, power/maximum(power))

# plot(p_init_state, p_power, layout=(2,1))
