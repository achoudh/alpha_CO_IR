# Line width
sσ(Δν) = Δν / (2.0*sqrt(2.0*log(2.0)))

# Gaussian function
gssn(ν, ν0, Δν) = exp(-(ν - ν0)^2 / (2 * sσ(Δν)^2)) / (sσ(Δν) * sqrt(2π))

# define hstatT
function hstatT(ev::Vector{Float64}, eu::Vector{Vector{Float64}}, com_ol::Matrix{Float64})

    h::Matrix{Float64} = zeros(Float64, nmols_ml, nmols_ml)

    for n1::Int64 in 1:nmols_ml

        for n2::Int64 in n1+1:nmols_ml

            rvec12::Vector{Float64} = com_ol[n1,:] - com_ol[n2,:]
            rvec12[1] = rvec12[1] - nx*round(Int, rvec12[1]/(nx))
            rvec12[2] = rvec12[2] - ny*round(Int, rvec12[2]/(ny))
            # rvec12[3] = rvec12[3] - ny*round(Int, rvec12[3]/(nz))

            rvec12 = a0_CO .* rvec12
            r12::Float64 = norm(rvec12)
            en::Vector{Float64} = rvec12/r12

            force::Float64 = (dot(eu[n1][:], eu[n2][:]) - 3.0*dot(en, eu[n1][:])*dot(en, eu[n2][:])) / r12^3
        
            h[n1,n1] += force
            h[n2,n2] += force
            h[n1,n2] = unit1*μ01^2*force
            h[n2,n1] = h[n1,n2]

        end

    end

    for n1::Int64 in 1:nmols_ml
        h[n1,n1] = ev[n1] + unit1*(μ11 - μ00)*μ00*h[n1,n1]
    end

    # return eigenvalues and eigenvectors
    return eigen(h)
       
end

# IR spectra
function ir_spectra(νk::Vector{Float64}, eu::Vector{Vector{Float64}}, com_ol::Matrix{Float64}, Δν)

    # θ = x[1+0*nmols_ml:1*nmols_ml]
    # ϕ = x[1+1*nmols_ml:2*nmols_ml]

    # Unperturbed eigenvalues
    ev::Vector{Float64} = fill(ν0, nmols_ml)
    
    # ev::Vector{Float64} =zeros(nmols_ml)
    # for i in 1:nmols_ml
    #     ev[i] = θ[i] < 0.5*pi ? ν0[1] : ν0[2]
    # end

    # # Orientation of the dipole moments of each vectors
    # eu = zeros(nmols_ml, 3) 
    # for i in 1:nmols_ml
    #     eu[i,:] = [sin(θ[i]) * cos(ϕ[i]), sin(θ[i]) * sin(ϕ[i]), cos(θ[i])]
    # end

    eigenvals, eigenvecs = hstatT(ev, eu, com_ol)
    σ = eigenvals ./ nmols_ml

    μEpda::Vector{Float64} = zeros(nmols_ml)
    μEsda::Vector{Float64} = zeros(nmols_ml)
    μEp::Vector{Float64}   = zeros(nmols_ml)
    μEs::Vector{Float64}   = zeros(nmols_ml)
    pl::Matrix{Float64}    = zeros(nmols_ml,3)
    for m::Int64 in 1:nmols_ml # loop over eigenvecs
        for i::Int64 in 1:nmols_ml # loop over molecules
            pl[m,:] += eigenvecs[i, m]*(μ01 .* eu[i][:])
        end
        # Single domain
        μEp[m] = dot(pl[m,:],ep)^2
        μEs[m] = dot(pl[m,:],es)^2
        # Domain average
        μEpda[m] = (0.5*((pl[m,1])^2 + (pl[m,2])^2)*Tx + (pl[m,3])^2 *Tz)
        μEsda[m] = 0.5*(pl[m,1]^2 + pl[m,2]^2)*Ty
        
    end

    ipda::Vector{Float64} = zeros(size(νk,1))
    isda::Vector{Float64} = zeros(size(νk,1))
    ip::Vector{Float64}   = zeros(size(νk,1))
    is::Vector{Float64}   = zeros(size(νk,1))

    for (iν,ν) in enumerate(νk)
        for m in 1:nmols_ml
            gp = gssn(ν, eigenvals[m], Δν)#1.15
            gs = gssn(ν, eigenvals[m], Δν)
            ipda[iν] += unit2*σ[m]*μEpda[m] * gp
            isda[iν] += unit2*σ[m]*μEsda[m] * gs
            ip[iν] += unit2*σ[m]*μEp[m] * gp
            is[iν] += unit2*σ[m]*μEs[m] * gs
        end
    end

    return ipda, isda, ip, is
end

