# Implementation of DifferentialEquations with callbacks

using DifferentialEquations, Plots

f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

using Plots
plot(sol, linewidth=5, title="Solution to the linear ODE with a thick line",
    xaxis="Time (t)", yaxis="u(t) (in μm)", label="My Thick Line!") # legend=false
plot!(sol.t, t -> 0.5 * exp(1.01t), lw=3, ls=:dash, label="True Solution!")

## callbacks
using DifferentialEquations
using Plots

function f(du, u, p, t)
    du[1] = -u[1]
end

u0 = [10.0]
const V = 1
prob = ODEProblem(f, u0, (0.0, 10.0))
sol = solve(prob, Tsit5())
plot(sol)

condition(u, t, integrator) = t == 4
affect!(integrator) = integrator.u[1] += 10
cb = DiscreteCallback(condition, affect!)
sol = solve(prob, Tsit5(), callback=cb)
plot(sol)
sol = solve(prob, Tsit5(), callback=cb, tstops=[4.0])
plot(sol)
dosetimes = [4.0, 8.0]
condition(u, t, integrator) = t ∈ dosetimes
affect!(integrator) = integrator.u[1] += 10
cb = DiscreteCallback(condition, affect!)
sol = solve(prob, Tsit5(), callback=cb, tstops=dosetimes)
plot(sol)
dosetimes = [4.0, 6.0, 8.0]
condition(u, t, integrator) = t ∈ dosetimes && (u[1] < 1.0)
affect!(integrator) = integrator.u[1] += 10integrator.t
cb = DiscreteCallback(condition, affect!)
sol = solve(prob, Tsit5(), callback=cb, tstops=dosetimes)
plot(sol)


# callback for boluses
function f(du, u, p, t)
    du[1] = -u[1]
end

u0 = [10.0]
prob = ODEProblem(f, u0, (0.0, 10.0))
dosetimes = [4.0, 8.0]
doses = [5.0, 12.0]
affect!(integrator) = integrator.u[1] += popfirst!(doses) # can probably do this in a better way
cb = PresetTimeCallback(dosetimes, affect!)
sol = solve(prob, Tsit5(), callback=cb)
plot(sol)

# constant piece-wise infusions ?


## pk3 model
# θ = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
# k10, k12, k13, k21, k31, V1 = θ
# A = [-(k10 + k12 + k13) k12 k13
#     k21 -k21 0.0
#     k31 0.0 -k31]

# B = [1 / V1; 0; 0]
# C = [1 0 0]
# D = 0

k10 = 1.0
k12 = 1.1
k21 = 1.2
k13 = 1.3
k31 = 1.4
V1 = 2.0

function pk3!(dx, x, p, t)
    # x1, x2, x3 = x
    dx[1] = -(k10 + k12 + k13) * x[1] + k12 * x[2] + k13 * x[3] + (1 / V1) * p(t)
    dx[2] = k21 * x[1] - k21 * x[2]
    dx[3] = k31 * x[1] - k31 * x[3]
end


u = t -> 1.0
x0 = [0.0, 0.0, 0.0]
tspan = [0.0, 100.0]
prob = ODEProblem(pk3!, x0, tspan, u)
sol = solve(prob)
plot(sol)

# constant input, works

# changing input u (here, actually a parameter) during the way -works
function pk3!(dx, x, p, t)
    # x1, x2, x3 = x
    dx[1] = -(k10 + k12 + k13) * x[1] + k12 * x[2] + k13 * x[3] + (1 / V1) * p
    dx[2] = k21 * x[1] - k21 * x[2]
    dx[3] = k31 * x[1] - k31 * x[3]
end


p = 1.0 #t -> 1.0
x0 = [0.0, 0.0, 0.0]
tspan = [0.0, 100.0]
prob = ODEProblem(pk3!, x0, tspan, p)
sol = solve(prob)
plot(sol)


infusionrates = [2.0, 4.0]
infusiontimes = [20.0, 40.0]
# affect!(integrator) = integrator.u[1] += popfirst!(bolusdoses) # can probably do this in a better way
# cb = PresetTimeCallback(bolusdosetimes, affect!)
affect!(integrator) = integrator.p += popfirst!(infusionrates) # can probably do this in a better way
cb = PresetTimeCallback(infusiontimes, affect!)

# Solve
sol = solve(prob, callback=cb)
plot(sol)



# now add boluses
function pk3!(dx, x, p, t)
    # x1, x2, x3 = x
    dx[1] = -(k10 + k12 + k13) * x[1] + k12 * x[2] + k13 * x[3] + (1 / V1) * p
    dx[2] = k21 * x[1] - k21 * x[2]
    dx[3] = k31 * x[1] - k31 * x[3]
end

p = 0.0 #t -> 1.0 # starting input u

x0 = [0.0, 0.0, 0.0]
tspan = [0.0, 100.0]
prob = ODEProblem(pk3!, x0, tspan, p)

bolusdoses = [5.0, 12.0]
bolusdosetimes = [30.0, 60.0]

infusionrates = [2.0, 4.0]
infusiontimes = [20.0, 40.0]

affectbolus!(integrator) = integrator.u[1] += popfirst!(bolusdoses) # can probably do this in a better way
# u[1], is that affecting x1? think so
cbbolus = PresetTimeCallback(bolusdosetimes, affectbolus!)

affectinfusion!(integrator) = integrator.p += popfirst!(infusionrates) # can probably do this in a better way
cbinfusion = PresetTimeCallback(infusiontimes, affectinfusion!)

cbs = CallbackSet(cbbolus, cbinfusion)

# Solve
sol = solve(prob, callback=cbs)
plot(sol)

# seems to work

# now add observations
function pk3!(dx, x, p, t)
    # x1, x2, x3 = x
    dx[1] = -(k10 + k12 + k13) * x[1] + k12 * x[2] + k13 * x[3] + (1 / V1) * p
    dx[2] = k21 * x[1] - k21 * x[2]
    dx[3] = k31 * x[1] - k31 * x[3]
end

p = 0.0 #t -> 1.0 # starting input u
x0 = [0.0, 0.0, 0.0]
tspan = [0.0, 100.0]
prob = ODEProblem(pk3!, x0, tspan, p)

bolusdoses = [5.0, 12.0]
bolusdosetimes = [30.0, 60.0]

infusionrates = [2.0, 4.0]
infusiontimes = [20.0, 40.0]

yobs = [10.0,20.0,23.0,60.0,82.0, 83.2, 84.5]

affectbolus!(integrator) = integrator.u[1] += popfirst!(bolusdoses) # can probably do this in a better way
# u[1], is that affecting x1? think so
cbbolus = PresetTimeCallback(bolusdosetimes, affectbolus!)

affectinfusion!(integrator) = integrator.p += popfirst!(infusionrates) # can probably do this in a better way
cbinfusion = PresetTimeCallback(infusiontimes, affectinfusion!)

cbs = CallbackSet(cbbolus, cbinfusion)

# Solve
sol = solve(prob, callback=cbs, saveat=yobs)
plot(sol)

# changing theta values
function pk3!(dx, x, p, t)
    # x1, x2, x3 = x
    k10,k12,k13,k21,k31,V1,w = p
    dx[1] = -(k10 + k12 + k13) * x[1] + k12 * x[2] + k13 * x[3] + (1 / V1) * w
    dx[2] = k21 * x[1] - k21 * x[2]
    dx[3] = k31 * x[1] - k31 * x[3]
end

θ = [1.0,1.1,1.2,1.3,1.4,2.0]
w = 0.0
p = [θ; w]
x0 = [0.0, 0.0, 0.0]
tspan = [0.0, 100.0]
prob = ODEProblem(pk3!, x0, tspan, p)

bolusdoses = [5.0, 12.0]
bolusdosetimes = [30.0, 60.0]

infusionrates = [2.0, 4.0]
infusiontimes = [20.0, 40.0]

yobs = [10.0,20.0,23.0,60.0,82.0, 83.2, 84.5]

affectbolus!(integrator) = integrator.u[1] += popfirst!(bolusdoses) # can probably do this in a better way
# u[1], is that affecting x1? think so
cbbolus = PresetTimeCallback(bolusdosetimes, affectbolus!)

affectinfusion!(integrator) = integrator.p[7] += popfirst!(infusionrates) # can probably do this in a better way
cbinfusion = PresetTimeCallback(infusiontimes, affectinfusion!)

cbs = CallbackSet(cbbolus, cbinfusion)

# Solve
sol = solve(prob, callback=cbs, saveat=yobs)
plot(sol)