module anajulia

using Integrals

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

function watch_point(wave, xo, k, dt, num_periods=1)
    speed = get_speed(k, dt)
    tf =  wave.l / speed
    
    watch_idx = Int(wave.l * xo)
    time = Vector{Float32}(0:dt:tf*num_periods)
    y_data = Vector{Float32}(undef, length(time))
    for i in 1:length(time)
        y_data[i] = wave.y[watch_idx]
        wave_step!(wave, k)
    end

    return time, y_data
end

function fourier_coefficients(y_data_in, time_in, c, l, num_freq)
    mask = time_in .< 4*l/c
    y_data = y_data_in[mask]
    time = time_in[mask]

    method = TrapezoidalRule()

    as = Vector{Float32}(undef, num_freq)
    bs = Vector{Float32}(undef, num_freq)
    freqs = Vector{Float32}(undef, num_freq)
    for i in 1:num_freq
        omega = (i-1) * π * c / l
        freqs[i] = omega / (2 * π)

        int_a = y_data .* cos.(omega * time)
        int_b = y_data .* sin.(omega * time)

        prob_a = SampledIntegralProblem(int_a, time)
        prob_b = SampledIntegralProblem(int_b, time)

        if i==1
            as[i] = 1/(time[end] - time[1]) * solve(prob_a, method).u
        else
            as[i] = c/l * solve(prob_a, method).u
        end

        bs[i] = c/l * solve(prob_b, method).u
    end
    
    return as, bs, freqs
end

function fourier_transform(y_data, time, max_freq, df)
    method = TrapezoidalRule()

    freq_arr = 0:df:max_freq
    num_freqs = length(freq_arr)

    as = Vector{Float32}(undef, num_freqs)
    bs = Vector{Float32}(undef, num_freqs)
    
    for (i, freq) in enumerate(freq_arr)
        omega = freq * 2 * π
        
        int_a = y_data .* cos.(omega * time)
        int_b = y_data .* sin.(omega * time)
        
        prob_a = SampledIntegralProblem(int_a, time)
        prob_b = SampledIntegralProblem(int_b, time)
    
        as[i] = 1/(2*π)^.5 * solve(prob_a, method).u
        bs[i] = 1/(2*π)^.5 * solve(prob_b, method).u
    end

    return as, bs, freq_arr
end

function get_power(as, bs)
    power = zeros(length(as))
    @. power = as^2 + bs^2
    return power
end

function forier_sum(time, as, bs, freqs)
    y = zeros(length(time))
    for i in 1:length(as)
        omega = freqs[i] * 2 * π
        # for j in 1:length(y)
        #     y[j] += as[i] * cos(omega * time[j]) + bs[i] * sin(omega * time[j])
        # end
        @. y += as[i] * cos(omega * time) + bs[i] * sin(omega * time)
    end

    return y
end

module init_state
    function gaussian_package(x, center_rel, sigma)
        y = Vector{Float32}(undef, length(x))
        
        center = x[1] + center_rel * (x[end] - x[1])
        @. y = exp(-sigma*(x - center)^2)
        return y
    end
    
    function pluck_string(x, center_rel, height)
        y = Vector{Float32}(undef, length(x))    
        
        center = x[1] + center_rel * (x[end] - x[1])
        mask = x .< center
        @. y[mask] = height/(center - x[1]) * (x[mask] -x[1])  
        @. y[.!mask] = -height/(x[end] -center) * (x[.!mask]-x[end]) 
        
        return y
    end
end

end