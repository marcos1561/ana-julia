using Plots

function wave_step!(y, y_last, k)
    l = length(y)
    y_copy = copy(y)
    for i in 2:l-1
        y[i] = 2 * (1 - k^2) * y_copy[i] - y_last[i] + k^2 * (y_copy[i+1] + y_copy[i-1])
        y_last[i] = y_copy[i]
    end
end

l = 100
k = 1
num_steps = 200

sigma = 0.05
x0 = l/4

y = Vector{Float32}(undef, l)
x = Vector{Float32}(0:l-1)
# @. y = exp(-sigma*(y - x0)^2)

max_id = Int(l/4)
@. y[1:max_id] = x[1:max_id]/x[max_id]  
@. y[max_id:end] = -1 * (x[max_id:end] - x[max_id]) * 1/(x[end] - x[max_id])  + 1

y_last = copy(y)
y_init = copy(y)

plot(y)
@gif for i in 1:num_steps
    wave_step!(y, y_last, k)
    plot(y_init, label="inicial")
    plot!(y, label="final")
    ylims!(-1, 1)
end



