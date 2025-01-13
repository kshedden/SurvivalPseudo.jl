using Test
using SurvivalPseudo
using RCall
using StableRNGs
using StatsBase
using Distributions
using Statistics
using UnicodePlots
using FiniteDifferences

@testset "Test survival statistics" begin

    rng = StableRNG(123)

    for m0 in [10, 20, 100]
        for q in [2, 5, 1000]
            for wt in [false, true]
                for _ in 1:10

                    m = rand(rng, Poisson(m0))
                    time = sample(rng, 1.:q, m)
                    status = [rand(rng) < 0.5 for _ in 1:m]
                    entry = rand(rng, m) .* time

                    weights = wt ? 0.1 .* 0.9*rand(rng, m) : ones(m)

                    sf = SurvivalFunction(entry, time, status; weights=weights)

                    @rput time
                    @rput status
                    @rput entry
                    @rput weights

                    R"
                    library(survival)
                    sf = survfit(Surv(entry, time, status) ~ 1, weights=weights)
                    nrisk = sf$n.risk
                    nevent = sf$n.event
                    ncensor = sf$n.censor
                    utime = sf$time
                    rsurv = sf$surv
                    rse = sf$std.err
                    "

                    @rget nrisk
                    @rget nevent
                    @rget utime
                    @rget ncensor
                    @rget rsurv
                    @rget rse

                    @test isapprox(utime, sf.utime)
                    @test isapprox(nevent, sf.nevent)
                    @test isapprox(ncensor, sf.ncensor)
                    @test isapprox(nrisk, sf.nrisk)
                    @test isapprox(rsurv, sf.surv, atol=1e-5, rtol=1e-5)
                    @test isapprox(rse, stderror(sf))
                end
            end
        end
    end
end

@testset "Test survival pseudo-observations via martingale residuals" begin

    rng = StableRNG(123)

    tgt_age = 60

    n = 100
    dtime = rand(rng, Exponential(80), n)
    ctime = rand(rng, Exponential(80), n)
    entry = rand(rng, Exponential(40), n)

    for i in 1:n
        a = sort([entry[i], ctime[i]])
        entry[i] = a[1]
        ctime[i] = a[2]
    end
    ii = dtime .>= entry
    dtime = dtime[ii]
    ctime = ctime[ii]
    entry = entry[ii]
    n = length(dtime)
    time = [min(x, y) for (x, y) in zip(dtime, ctime)]
    status = time .== dtime

    # Get the pseudo-observations and the analytic gradient
    sf = SurvivalFunction(entry, time, status)
    ps, est, _ = SurvivalPseudo.pseudo(sf, tgt_age; method=:mart)

    # Brute-force jack-knife
    n = length(time)
    bf = zeros(n)
    for i in 1:n
        ii = [j for j in 1:n if j != i]
        sf = SurvivalFunction(entry[ii], time[ii], status[ii])
        bf[i] = n*est - (n-1)*SurvivalPseudo.survprob(sf, tgt_age)
    end

    @test cor(ps, bf) > 0.98
    @test abs(mean(ps) - mean(bf)) < 0.05
    @test abs(std(ps) - std(bf)) < 0.25
end

@testset "Test survival pseudo-observations via infinitesimal jackknife" begin

    rng = StableRNG(123)

    tgt_age = 60

    n = 100
    dtime = rand(rng, Exponential(80), n)
    ctime = rand(rng, Exponential(80), n)
    entry = rand(rng, Exponential(40), n)

    for i in 1:n
        a = sort([entry[i], ctime[i]])
        entry[i] = a[1]
        ctime[i] = a[2]
    end
    ii = dtime .>= entry
    dtime = dtime[ii]
    ctime = ctime[ii]
    entry = entry[ii]
    n = length(dtime)
    time = [min(x, y) for (x, y) in zip(dtime, ctime)]
    status = time .== dtime

    # Get the numerical gradient
    f = function(w)
        sf = SurvivalFunction(entry, time, status; weights=w)
        return log(SurvivalPseudo.survprob(sf, tgt_age))
    end
    ngr = grad(central_fdm(5, 1), f, ones(n))[1]

    # Get the pseudo-observations and the analytic gradient
    sf = SurvivalFunction(entry, time, status)
    ps, est, agr = SurvivalPseudo.pseudo(sf, tgt_age; method=:IJ)

    # Check that the analytic gradient and numeric gradient are close
    @test isapprox(ngr, agr, rtol=0.001, atol=0.001)

    # Brute-force jack-knife
    n = length(time)
    bf = zeros(n)
    for i in 1:n
        ii = [j for j in 1:n if j != i]
        sf = SurvivalFunction(entry[ii], time[ii], status[ii])
        bf[i] = n*est - (n-1)*log(SurvivalPseudo.survprob(sf, tgt_age))
    end

    @test cor(ps, bf) > 0.99
    @test abs(mean(ps) - mean(bf)) < 0.05
    @test abs(std(ps) - std(bf)) < 0.2
end
