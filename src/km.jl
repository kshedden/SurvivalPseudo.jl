struct SurvivalFunction

    # The time at which each subject's state is first known.
    entry_orig::Vector{Float64}

    # The entry times rounded up to the next value in 'time'.  A subject is considered at risk on
    # the time in 'entry' (this different from R).
    entry::Vector{Float64}

    # The time at which each subject's state is last known, must be strictly larger than 'entry'.
    time::Vector{Float64}

    # If true, the event occurs at the corresponding element of 'time'
    status::Vector{Bool}

    # The sorted, unique values in 'time'
    utime::Vector{Float64}

    # The number of events occuring at each value in 'utime'
    nevent::Vector{Float64}

    # The number of observations that are censored at each value
    # in 'utime'
    ncensor::Vector{Float64}

    # The number of individuals at risk at each value of 'utime'
    nrisk::Vector{Float64}

    # The indices of observations that enter at each element of 'utime'
    enter::Vector{Vector{Int64}}

    # Weights for each observation
    weights::Vector{Float64}

    utime_ix::Vector{UnitRange}

    # The estimated survival probabilities (Kaplan-Meier or product-limit estimates)
    surv::Vector{Float64}

    # The is a permutation that reverses the sorting that is done to the raw data.
    irev::Vector{Int64}
end

function show(io::IO, sf::SurvivalFunction)
    (; time, utime) = sf
    m, n = length(time), length(utime)
    println(io, "SurvivalFunction with $(m) observations and $(n) unique times.")
end

function SurvivalFunction(entry, time, status; weights=nothing)

    if isnothing(weights)
        weights = ones(length(time))
    end

    # Sort by ascending time
    ii = sortperm(time)
    irev = sortperm(ii)
    time = time[ii]
    entry = entry[ii]
    status = status[ii]
    weights = weights[ii]

    if any(entry .>= time)
        error("Entry times must occur before event time")
    end

    # Unique times
    utime = sort(unique(time))

    # Number of events at each time in 'utime'.
    nevent = zeros(length(utime))

    # Number of individuals who are censored at each time in 'utime'
    ncensor = zeros(length(utime))

    # Ranges in 'time' corresponding to each unique event time
    utime_ix = []
    for (i,t) in enumerate(utime)
        ii = searchsorted(time, t)
        push!(utime_ix, ii)
        nevent[i] = sum(weights[ii] .* status[ii])
        ncensor[i] = sum(weights[ii]) - nevent[i]
    end

    # Shift entry times to subsequent event time, so that entry
    # times are a subset of event times
    enter = [zeros(Int, 0) for _ in eachindex(utime)]
    entry_orig = copy(entry)
    for i in eachindex(entry)
        ii = searchsorted(utime, entry[i])
        jj = length(ii) == 0 ? first(ii) : last(ii) + 1
        entry[i] = jj <= length(utime) ? utime[jj] : Inf
        push!(enter[jj], i)
    end

    m = length(utime)
    nrisk = zeros(m)
    for i in eachindex(utime)
        if i == 1
            nrisk[i] = sum(weights[enter[i]])
        else
            nrisk[i] = nrisk[i-1] + sum(weights[enter[i]]) - nevent[i-1] - ncensor[i-1]
        end
    end

    surv = cumprod(1 .- nevent ./ nrisk)

    return SurvivalFunction(entry_orig, entry, time, status, utime, nevent, ncensor, nrisk, enter, weights, utime_ix, surv, irev)
end

function survprob(sf::SurvivalFunction, stime)

    (; utime, surv) = sf

    _, i = findmin(abs.(utime .- stime))
    return surv[i]
end

"""
    pseudo(sf, stime; method)

Compute pseudo-observations for the (log) survival function `sf` at time `stime`.

If method is :IJ, the infintessimal jacknife method is used to obtain pseudo-observations
for the log survival probability.  If method is :mart, the martingale residual method
is used to obtain pseudo-observations for the survival probability.
"""
function pseudo(sf::SurvivalFunction, stime::T; method::Symbol=:IJ) where{T<:Real}

    if method == :IJ
        return pseudo_ij(sf, stime)
    elseif method == :mart
        return pseudo_mart(sf, stime)
    else
        error("Unknown method $(method)")
    end
end

"""
    stderror(sf)

Returns the standard errors of the (possibly transformed) estimated survival
probabilities.

If method == :greenwood, the returned standard errors are for the estimated
survival probabilities.  If method == :log, the returned standard errors
are for the log of the estimated survival probabilities.
"""
function stderror(sf::SurvivalFunction; method=:greenwood)

    (; nrisk, nevent, surv) = sf

    r = nevent ./ (nrisk .* (nrisk - nevent))
    c = sqrt.(cumsum(r))

    if method == :greenwood
        return surv .* c
    elseif method == :log
        return c
    end
end

"""
    confint(sf; [level=0.95])

Returns a confidence band for the survival function.
"""
function confint(sf::SurvivalFunction; level::Real=0.95)
    (; surv) = sf
    s = stderror(sf; method=:log)
    f = quantile(Normal(), 1 - (1 - level) / 2)
    return (lower=surv.*exp.(s - f*s), upper=surv.*exp.(s + f*s))
end

function pseudo_ij(sf::SurvivalFunction, stime)

    (; time, status, entry, utime, nrisk, nevent, utime_ix, surv, irev) = sf

    n = length(time)
    m = length(utime)

    _, js = findmin(abs.(utime .- stime))
    stime = utime[js]
    lest = log(surv[js])

    # dgrad[j] is the partial derivative of (1 - nevent[j]/nrisk[j]) with respect to nevent[j]
    # cgrad[j] is the partial derivative of cumsum (1 - nevent[j]/nrisk[j]) with respect to nrisk[j]
    dgrad = zeros(m)
    cgrad = zeros(m)
    for j in 1:m
        dgrad[j] = 1 / (nevent[j] - nrisk[j])
        cgrad[j] = nevent[j] / (nrisk[j] * (nrisk[j] - nevent[j]))
    end
    cgrad = cumsum(cgrad)

    # grad[i] is the gradient of the estimand with respect to the weight for observation i.
    grad = zeros(n)
    for i in 1:n
        ke = searchsortedfirst(utime, entry[i])
        if ke > js
            continue
        end
        k = searchsortedfirst(utime, time[i])
        if status[i] && k <= js
            grad[i] = dgrad[k]
        end
        k = min(k, js)
        grad[i] += cgrad[k]
        if ke > 1
            grad[i] -= cgrad[ke-1]
        end
    end
    grad = grad[irev]

    # The pseudo-observations
    loo = lest .- grad

    return (pseudo=n*lest .- (n-1)*loo, estimate=lest, grad=grad)
end

function pseudo_mart(sf::SurvivalFunction, stime)

    (; time, status, entry, utime, nrisk, nevent, utime_ix, surv, irev) = sf

    n = length(time)
    m = length(utime)

    ii = searchsorted(utime, stime)
    if last(ii) == 0
        js = first(ii)
    elseif first(ii) == m + 1
        js = last(ii)
    elseif first(ii) == last(ii)
        js = first(ii)
    else
        js = abs(utime[first(ii)] - stime) < abs(utime[last(ii)] - stime) ? first(ii) : last(ii)
    end
    est = surv[js]

    cx = n * cumsum(nevent ./ nrisk.^2)
    ps = zeros(n)

    for i in 1:n

        # Lower limit of integration
        k1 = searchsortedfirst(utime, entry[i])
        if k1 > js
            continue
        end

        # Upper limit of integration
        kt = searchsortedfirst(utime, time[i])
        k2 = min(kt, js)

        f = cx[k2]
        if k1 > 1
            f -= cx[k1-1]
        end
        f = -f
        if status[i] && kt <= js
            f += n / nrisk[kt]
        end

        ps[i] = surv[js] * (1 - f)
    end
    ps = ps[irev]

    return (pseudo=ps, estimate=est, grad=nothing)
end

"""
    meanreslife(sf, tau, ftype)

Pseudo-observations for mean restricted life.  The mean is restricted at tau.
The vector `ftype` indicates which observations experience the failure type
of interest.

References:

Xin Wang, Xiaoming Xue, and Liuquan Sun.  Regression analysis of restricted mean
survival time based on pseudo-observations for competing risks data.
COMMUNICATIONS IN STATISTICS—THEORY AND METHODS
2018, VOL. 47, NO. 22, 5614-5625
"""
function meanreslife(sf::SurvivalFunction, tau, ftype::Vector{Bool}; npt=20)

    (; time, status, entry_orig, entry, utime, nrisk, nevent, utime_ix, surv, irev) = sf

    if length(ftype) != length(time)
        error("The failure type vector 'ftype' must have the same length as 'time'")
    end

    # Integrate over this grid to estimate the MRL
    stimes = collect(range(0, tau, npt))
    mesh = stimes[2] - stimes[1]

    # Reverse survival function.
    csf = SurvivalFunction(entry_orig, time, .!status)

    n = length(time)

    # Evaluate the reverse survival function at each person's event time, for people
    # who are not censored.
    G = zeros(n)
    for i in 1:n
        if status[i]
            ii = searchsorted(utime, time[i])
            @assert first(ii) == last(ii)
            G[i] = csf.surv[first(ii)]
        else
            G[i] = NaN
        end
    end

    qs = zeros(length(stimes))
    for j in eachindex(stimes)
        a = 0.0
        m = 0
        for i in 1:n

            if entry[i] >= stimes[j]
                continue
            end

            m += 1
            if status[i] && ftype[i] && (time[i] > stimes[j])
                a += 1 / G[i]
            end
        end

        if m > 0
            qs[j] = a / m
        end
    end

    # Point estimate of the mean restricted life.
    tgt = sum(qs) * mesh

    ps = [pseudo(csf, t; method=:mart).pseudo for t in stimes]
    ps = hcat(ps...)

    # Jack-knife version of qs
    qsj = zeros(n, length(stimes))
    for j in eachindex(stimes)
        m = 0
        for i in 1:n
            if entry[i] >= stimes[j]
                continue
            end
            m += 1
            if status[i] && ftype[i] && (time[i] > stimes[j])
                # This is an approximation that does not jackknife G.
                qsj[i, j] = 1 / G[i]
            end
        end
        qsj[:, j] = sum(qsj[:, j]) .- qsj[:, j]
        qsj[:, j] ./= (m - 1)
        qsj[:, j] = m*qs[j] .- (m - 1)*qsj[:, j]
    end

    tgt_ps = sum(qsj; dims=2)[:] * mesh

    return (mrl=tgt, mrl_ps=tgt_ps)
end
