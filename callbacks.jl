using DifferentialEquations

function affect!(integrator, t_inv, start_time, n_invaders, n_species)
    println("affect outside")
    if integrator.t >= start_time
        println("affect inside")
        idx = Int(n_species + floor((integrator.t - start_time) / t_inv)) + 1
        if idx <= n_invaders + n_species
            println(idx)
            integrator.u[idx] += 10
        end
    end
end

function create_callbacks(t_inv, start_time, n_invaders, n_species)
    cb = PeriodicCallback(
        integrator -> affect!(integrator, t_inv, start_time, n_invaders, n_species), t_inv
    )
    return cb
end

