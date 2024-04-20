using Plots

include("anajulia.jl")
import .anajulia

file_name = "anim.gif"

# Configurações da corda
int_pars = anajulia.IntPars(
    k  = 1/4,
    c  = 880,
    dx = 0.01,
    b  = 5,
    epsilon = 1e-9*1,
    l = 1,
)

# Configurações inicias (gaussian)
center = 0.5
sigma = 0.01
height = 1

# Configurações da animação
fps = 30
time = 17
sim_time = 0.02


x = Vector{Float32}(0:int_pars.n-1)
# y = anajulia.init_state.gaussian_package(
#     x, center, sigma)
y = anajulia.init_state.pluck_string(
    x, center, height)


num_steps = trunc(Int, sim_time / int_pars.dt)
num_frames = trunc(Int, fps*time)
num_steps_per_frame = trunc(Int, num_steps / num_frames)

wave = anajulia.Wave(y)

int_pars_ideal = anajulia.IntPars(
    k=int_pars.k, c=int_pars.c , dx=int_pars.dx, 
    b=0, epsilon=0, l=int_pars.n*int_pars.dx,
)
wave_ideal = anajulia.Wave(y)

y_init = copy(wave.y)
time = 0
anim = @animate for i in 1:num_frames
    for j in 1:num_steps_per_frame
        anajulia.wave_step!(wave, int_pars)
        anajulia.wave_step!(wave_ideal, int_pars_ideal)
        global time += int_pars.dt
    end
    t_str = string(round(time * 1000, digits=2))
    plot(y_init, label="Inicial", title="t=$(t_str) ms")
    ylims!(-height, height)
    xlims!(-3, 130)
    xlabel!("(cm)")
    ylabel!("(cm)")
    plot!(wave_ideal.y, label="Ideal")
    plot!(wave.y, label="Realista")
end 
# fps=fps

gif(anim, "animacoes/$(file_name)", fps=fps)