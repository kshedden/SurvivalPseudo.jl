using RDatasets
using SurvivalPseudo
using UnicodePlots
using Statistics

# Load the lung data
lung = RDatasets.dataset("survival", "lung")

# Create a survival function from the lung data.
time = lung[:, :Time]
status = lung[:, :Status] .== 2
entry = zeros(length(time))
n = length(time)
sf = SurvivalFunction(entry, time, status)

# Compute pseudo-observations at these times.
times = range(10, 600, 10)

function plot_surv(times, psm, pse; log=true, title="")

    if log
        p = lineplot(times, exp.(psm), xlim=(0, 600), ylim=(0, 1), xlabel="Time", ylabel="Survival probability", title=title)
        lineplot!(p, times, exp.(psm - 2*pse))
        lineplot!(p, times, exp.(psm + 2*pse))
    else
        p = lineplot(times, psm, xlim=(0, 600), ylim=(0, 1), xlabel="Time", ylabel="Survival probability", title=title)
        lineplot!(p, times, psm - 2*pse)
        lineplot!(p, times, psm + 2*pse)
    end
    println(p)
end

# Compute the pseudo-observations for log survival probabilities, using infintesimal jackknife.
ps = [pseudo(sf, t).pseudo for t in times]
ps = hcat(ps...)
psm = mean(ps; dims=1)[:]
psd = std(ps; dims=1)[:]
pse = psd / sqrt(n)
plot_surv(times, psm, pse; log=true, title="Infinitesimal jackknife")

# Compute the pseudo-observations for log survival probabilities, using martingale residuals.
ps = [pseudo(sf, t; method=:mart).pseudo for t in times]
ps = hcat(ps...)
psm = mean(ps; dims=1)[:]
psd = std(ps; dims=1)[:]
pse = psd / sqrt(n)
plot_surv(times, psm, pse; log=false, title="Martingale residuals")

# Additive confidence band based on Greenwood's formula.
plot_surv(sf.utime, sf.surv, stderror(sf); log=false, title="Greenwood")

# Transform a confidence band for log survival.
plot_surv(sf.utime, log.(sf.surv), stderror(sf; method=:log); log=true, title="log")
