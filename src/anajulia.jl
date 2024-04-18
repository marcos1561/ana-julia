module anajulia

using Integrals

struct Wave
    y::Vector{Float32}
    y_last::Vector{Float32}
    l::Int

    Wave(y) = new(copy(y), copy(y), length(y)) 
end

mutable struct IntPars
    k::Float32
    c::Float32
    dx::Float32
    dt::Float32

    b::Float32
    epsilon::Float32
    n::Int

    a::Vector{Float32}

    IntPars(;k, c, dx, b, epsilon, l) = (
        # dt = get_dt(dx, c);
        dt = k * dx / c;
        n = trunc(Int, l/dx);

        d = 1 + b*dt;
        
        a = Vector{Float32}(undef, 4);
        a[1] = (2 - 2*k^2 - 6*epsilon*(n*k)^2)/d;
        a[2] = (-1 + b*dt)/d;
        a[3] = k^2*(1+4*epsilon*n^2)/d;
        a[4] = (epsilon*(n*k)^2)/d;

        new(k, c, dx, dt, b, epsilon, n, a)
    )
end

Base.copy(p::IntPars) = IntPars(
    k=p.k, c=p.c, dx=p.dx, b=p.b, epsilon=p.epsilon, l=p.n*p.dx)


function wave_step!(wave, pars::IntPars)
    y = wave.y
    y_last = wave.y_last

    l = length(y)
    
    y_copy = Vector{Float32}(undef, l+2)
    y_copy[2:end-1] .= copy(y)
    y_copy[1] = -y[2]
    y_copy[end] = -y[end-1]
    
    # println("Real")
    a = pars.a
    for i in 2:l-1
        ci = i+1

        term1 = a[1] * y_copy[ci] + a[2] * y_last[i] 
        term2 = a[3] * (y_copy[ci+1] + y_copy[ci-1])
        term3 = a[4] * (y_copy[ci+2] + y_copy[ci-2])
        
        # println("$(i): $(term1) | $(term2) | $(term3)")
        y[i] = term1 + term2 + term3
        
        y_last[i] = y_copy[ci]
    end
    
    # k = pars.k
    # y_copy = copy(y)
    # for i in 2:l-1
    #     y[i] = 2 * (1 - k^2) * y_copy[i] - y_last[i] + k^2 * (y_copy[i+1] + y_copy[i-1])
    #     y_last[i] = y_copy[i]
    # end
end

function wave_step_ideal!(wave, pars::IntPars)
    y = wave.y
    y_last = wave.y_last
    
    l = length(y)
    
    k = pars.k
    y_copy = copy(y)
    # println("Ideal")
    for i in 2:l-1
        # y[i] = 2 * (1 - k^2) * y_copy[i] - y_last[i] + k^2 * (y_copy[i+1] + y_copy[i-1])
        
        term1 = 2 * (1 - k^2) * y_copy[i] - y_last[i]
        term2 = k^2 * (y_copy[i+1] + y_copy[i-1])
        
        # println("$(i): $(term1) | $(term2)")
        
        y[i] = term1 + term2
    
        y_last[i] = y_copy[i]
    end
end

function get_speed(k, dt)
    return k / dt     
end

function get_dt(dx, speed)
    return (speed * 1/dx)^-1
end

function watch_point(wave, xo, k, dt, num_periods=1)
    speed = get_speed(k, dt)
    tf = wave.l / speed
    
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