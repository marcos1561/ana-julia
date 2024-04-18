using Plots

include("anajulia.jl")
import .anajulia

# Configurações da corda
int_pars = anajulia.IntPars(
    k  = 1/4,
    c  = 500,
    dx = 0.01,
    b  = 10,
    epsilon = 1e-8,
    l = 1,
)

println(int_pars.dt)

# Configurações inicias (gaussian)
center = 0.3
sigma = 0.01
height = 1


x = Vector{Float32}(0:int_pars.n-1)
# y = anajulia.init_state.gaussian_package(
#     x, center, sigma)
y = anajulia.init_state.pluck_string(
    x, center, 1)


fps = 30
time = 10
sim_time = 0.02

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
@gif for i in 1:num_frames
    for j in 1:num_steps_per_frame
        anajulia.wave_step!(wave, int_pars)
        anajulia.wave_step!(wave_ideal, int_pars_ideal)
        global time += int_pars.dt
    end
    t_str = string(round(time * 1000, digits=2))
    plot(y_init, label="inicial", title="t=$(t_str) ms")
    ylims!(-height, height)
    plot!(wave_ideal.y, label="ideal")
    plot!(wave.y, label="realista")
end fps=fps