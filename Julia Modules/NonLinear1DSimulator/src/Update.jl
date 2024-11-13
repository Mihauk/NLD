module Update

export update_partially_move_stuck_waiting_times!, single_upd!

using .Config
using Random

"""
    update_partially_move_stuck_waiting_times!(config::Conf, p::Vector{Bool}, mpl::Float64, mph::Float64)

Update the configuration `config` in place, according to the stochastic rules defined for the model. `p` is a vector of random boolean values for deciding particle movement, `mpl` and `mph` are probabilities of not moving for low and high density triplets respectively.
"""
function update_partially_move_stuck_waiting_times!(config::Conf, p::Vector{Bool}, mpl::Float64, mph::Float64)
    l = length(p)
    st = config.st
    for j in 1:l
        c_f = rand(Bool)
        move_prob_high = rand()
        move_prob_low = rand()
        idx = 3*(j-1)

        if (st[idx+1]==true && st[idx+2]==false && st[idx+3]==false && move_prob_low > mpl)
            @inbounds st[idx+1] = false
            @inbounds st[idx+2] = p[j]
            @inbounds st[idx+3] = !p[j]
        elseif (st[idx+1]==false && st[idx+2]==true && st[idx+3]==false && move_prob_low > mpl)
            @inbounds st[idx+1] = p[j]
            @inbounds st[idx+2] = false
            @inbounds st[idx+3] = !p[j]
        elseif (st[idx+1]==false && st[idx+2]==false && st[idx+3]==true && move_prob_low > mpl)
            @inbounds st[idx+1] = p[j]
            @inbounds st[idx+2] = !p[j]
            @inbounds st[idx+3] = false
        elseif (st[idx+1]==false && st[idx+2]==true && st[idx+3]==true && move_prob_high > mph)
            if c_f == false
                @inbounds st[idx+1] = true
                @inbounds st[idx+2] = false
                @inbounds st[idx+3] = true
            else
                @inbounds st[idx+1] = true
                @inbounds st[idx+2] = true
                @inbounds st[idx+3] = false
            end
        elseif (st[idx+1]==true && st[idx+2]==true && st[idx+3]==false && move_prob_high > mph)
            if c_f == true
                @inbounds st[idx+1] = true
                @inbounds st[idx+2] = false
                @inbounds st[idx+3] = true
            else
                @inbounds st[idx+1] = false
                @inbounds st[idx+2] = true
                @inbounds st[idx+3] = true
            end
        elseif (st[idx+1]==true && st[idx+2]==false && st[idx+3]==true && move_prob_high > mph)
            @inbounds st[idx+1] = p[j]
            @inbounds st[idx+2] = true
            @inbounds st[idx+3] = !p[j]
        end
    end
    return
end

"""
    single_upd!(conf::Conf, off::Int, mpl::Float64, mph::Float64)

Perform a single update of the configuration `conf`, with an offset `off`. The function applies a circular shift, updates the configuration, and shifts it back.
"""
function single_upd!(conf::Conf, off::Int, mpl::Float64, mph::Float64)
    l = length(conf.st)
    rolled_st = circshift(conf.st, -off)
    rolled_conf = Conf(rolled_st)
    p = rand(Bool, div(l, 3))
    update_partially_move_stuck_waiting_times!(rolled_conf, p, mpl, mph)
    conf.st = circshift(rolled_conf.st, off)
    return
end

end # module
