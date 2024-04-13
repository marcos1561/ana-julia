using Plots, Integrals

struct Wave
    y::Vector{Float32}
    y_last::Vector{Float32}
    l::Int

    Wave(y) = new(y, copy(y), length(y)) 
end

function wave_step!(wave, k)
    y = wave.y
    y_last = wave.y_last

    l = length(y)
    y_copy = copy(y)
    for i in 2:l-1
        y[i] = 2 * (1 - k^2) * y_copy[i] - y_last[i] + k^2 * (y_copy[i+1] + y_copy[i-1])
        y_last[i] = y_copy[i]
    end
end

function get_speed(k, dt)
    return k / dt    
end

function watch_point(wave, xo, k, dt)
    speed = get_speed(k, dt)
    tf =  wave.l / speed
    
    watch_idx = Int(wave.l * xo)
    time = Vector{Float32}(0:dt:tf)
    y_data = Vector{Float32}(undef, length(time))
    for i in 1:length(time)
        y_data[i] = wave.y[watch_idx]
        wave_step!(wave, k)
    end

    return time, y_data
end

function get_amplitudes(y_data, time, c, l, num_freq)
    method = TrapezoidalRule()

    as = Vector{Float32}(undef, num_freq)
    bs = Vector{Float32}(undef, num_freq)
    freqs = Vector{Float32}(undef, num_freq)
    for i in 1:num_freq
        omega = i * π * c / l
        freqs[i] = omega / (2 * π)

        int_a = y_data .* cos.(i * π * c / l * time)
        int_b = y_data .* sin.(i * π * c / l * time)

        prob_a = SampledIntegralProblem(int_a, time)
        prob_b = SampledIntegralProblem(int_b, time)

        as[i] = c/l * solve(prob_a, method).u
        bs[i] = c/l * solve(prob_b, method).u
    end

    return as, bs, freqs
end

function get_power(as, bs)
    power = zeros(length(as))
    @. power = as^2 * bs^2
    return power
end

function forier_sum(time, as, bs, c, l)
    y = zeros(length(time))
    for i in 1:length(as)
        omega = i * π * c / l
        # for j in 1:length(y)
        #     y[j] += as[i] * cos(omega * time[j]) + bs[i] * sin(omega * time[j])
        # end
        @. y += as[i] * cos(omega * time) + bs[i] * sin(omega * time)
    end

    return y
end


l = 500
k = 1
dt = 1
num_steps = 40

sigma = 0.05
x0 = l/4

y = Vector{Float32}(undef, l)
x = Vector{Float32}(0:l-1)

@. y = exp(-sigma*(x - x0)^2)

# max_id = Int(l/4)
# @. y[1:max_id] = x[1:max_id]/x[max_id]  
# @. y[max_id:end] = -1 * (x[max_id:end] - x[max_id]) * 1/(x[end] - x[max_id])  + 1


wave = Wave(y)

speed = get_speed(k, dt)

time, y_data = watch_point(wave, 0.3, k, dt)
plot(time, y_data)

as, bs, freqs = get_amplitudes(y_data, time, speed, wave.l, 100)
power = get_power(as, bs)

y_forier = forier_sum(time, as, bs, speed, wave.l)

# plot(freqs, power)
plot!(time, y_forier)

# plot(wave.y, label="inicial")
# for i in 1:num_steps
#     wave_step!(wave, k)
# end
# plot!(wave.y, label="final")



# y_last = copy(y)
# y_init = copy(y)

# plot(y)
# @gif for i in 1:num_steps
#     wave_step!(y, y_last, k)
#     plot(y_init, label="inicial")
#     plot!(y, label="final")
#     ylims!(-1, 1)
# end



