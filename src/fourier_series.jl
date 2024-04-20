using Plots

include("anajulia.jl")
import .anajulia

# Configurações da corda
int_pars = anajulia.IntPars(
    k = 1,
    c = 300,
    dx = 0.01,
    b = 0,
    epsilon = 0,
    l = 1,
)


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

time, waveform = anajulia.watch_point(wave, point_pos, int_pars, 2)

speed = anajulia.get_speed(k, dt)

wave_length = int_pars.n * int_pars.dx
# anim = @animate for num_f in 1:num_freqs
@gif for num_f in 1:num_freqs
    as, bs, freqs = anajulia.fourier_coefficients(waveform, time, int_pars.c, wave_length, num_f)
    y_forier = anajulia.forier_sum(time, as, bs, freqs)
    
    # plot(time, waveform, label="Sinal")
    # plot!(time, y_forier, label="Série de fourier ($(num_freqs) frequências)")
    plot(time, waveform, legend=false, title="N=$(num_f-1)")
    ylims!(-0.2, 0.5)
    plot!(time, y_forier)
end fps=5
# gif(anim, "animacoes/serie_fourier.gif", fps=5)
