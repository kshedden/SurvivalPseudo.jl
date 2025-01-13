using SurvivalPseudo
using UnicodePlots
using Statistics
using StatsBase
using CSV
using DataFrames
using Loess

# The myeloma data from the R survival package.
mdat = CSV.read("myeloma.csv.gz", DataFrame)
n = size(mdat, 1)

sf = SurvivalFunction(mdat[:, :entry], mdat[:, :futime], mdat[:, :death])

# Compute pseudo-observations at these times.
times = range(10, 8000, 10)

# Compute the pseudo-observations for log survival probabilities, using infintesimal jackknife.
ps = [pseudo(sf, t).pseudo for t in times]
ps = hcat(ps...)
psm = mean(ps; dims=1)[:]
psd = std(ps; dims=1)[:]

# The standard error from the pseudo-observations
pse = psd / sqrt(n)

# Aymptotic standard error
se = stderror(sf; method=:log)

# Compare the pseudo-observation and asymptotic standard errors.
p = lineplot(sf.utime, se; xlabel="futime", ylabel="SE(log surv)")
lineplot!(p, times, pse)
println(p)

# Obtain mean restricted life pseudo-observations out to day 4000.
ftype = sf.status # There is only one type of failure here.
ss = meanreslife(sf, 365.25*5, ftype)

# Assess how MRL relates to year (it is increasing due to better care)
println(corkendall(ss.mrl_ps, mdat[:, :year]))
println(cor(ss.mrl_ps, mdat[:, :year]))

dm = DataFrame(ps=ss.mrl_ps, year=mdat[:, :year])
ml = loess(dm[:, :year], dm[:, :ps])
xx = range(extrema(dm[:, :year])..., 100)
yy = predict(ml, xx)
p = lineplot(xx, yy)
println(p)
