# Line width
sσ(Δν) = Δν / (2.0*sqrt(2.0*log(2.0)))

# Gaussian function
gssn(ν, ν0, Δν) = exp(-(ν - ν0)^2 / (2 * sσ(Δν)^2)) / (sσ(Δν) * sqrt(2π))

# define hstatT
function hstatT_centre(ev::Vector{Float64}, eu::Vector{Vector{Float64}}, com_ol::Matrix{Float64}, large_h_matrix_sum)

    h::Matrix{Float64} = zeros(Float64, nmols_centre, nmols_centre)
    count_n2 = 0
    for n1::Int64 in 1:nmols_centre
        
        for n2::Int64 in n1+1:nmols_centre
            rvec12::Vector{Float64} = com_ol[n1,:] - com_ol[n2,:]
            rvec12[1] = rvec12[1] - nx*round(Int, rvec12[1]/(nx)) # periodic boundary conditions in x
            rvec12[2] = rvec12[2] - ny*round(Int, rvec12[2]/(ny)) # periodic boundary conditions in y
            if z_boundary_conditions == true
                rvec12[3] = rvec12[3] - ny*round(Int, rvec12[3]/(nz)) # periodic boundary conditions in z
            end
            
            rvec12 = a0_CO .* rvec12
            r12::Float64 = norm(rvec12)
            en::Vector{Float64} = rvec12/r12
            if Interaction_radius_cutoff == true
                if Interaction_centre_inside_unit_cell == false
                    if norm(com_ol[n1,:]*a0_CO - position_z_centre) >= interaction_cut_off_radius * a0_CO || norm(com_ol[n2,:]*a0_CO - position_z_centre) >= interaction_cut_off_radius * a0_CO # r12 >= interaction_cut_off_radius * a0_CO  || norm(com_ol[n1,:]*a0_CO - position_z_centre) >= interaction_cut_off_radius * a0_CO || norm(com_ol[n2,:]*a0_CO - position_z_centre) >= interaction_cut_off_radius * a0_CO
                        force = 0.0   # force::Float64
                    else
                        force = (dot(eu[n1][:], eu[n2][:]) - 3.0*dot(en, eu[n1][:])*dot(en, eu[n2][:])) / r12^3
                        count_n2 += 1
                        #print("n1:",n1," n2:",n2," ")
                    end
                else # Interaction_centre_inside_unit_cell == true
                    if norm(com_ol[n1,:]*a0_CO - position_z_centre) > interaction_cut_off_radius * a0_CO || norm(com_ol[n2,:]*a0_CO - position_z_centre) > interaction_cut_off_radius * a0_CO # either n1 or n2 are more than interaction_cut_off_radius away from the centre
                        force = 0.0
                    else
                        force = (dot(eu[n1][:], eu[n2][:]) - 3.0*dot(en, eu[n1][:])*dot(en, eu[n2][:])) / r12^3
                        count_n2 += 1
                        #print(" r12:", r12)
                        #print("n1:",n1," n2:",n2," ")
                    end
                end
            else
                force = (dot(eu[n1][:], eu[n2][:]) - 3.0*dot(en, eu[n1][:])*dot(en, eu[n2][:])) / r12^3
                count_n2 += 1
            end

            h[n1,n1] += force
            h[n2,n2] += force
            h[n1,n2] = unit1*μ01^2*force
            h[n2,n1] = h[n1,n2]

        end
    end

    for n1::Int64 in 1:nmols_centre
        h[n1,n1] = large_h_matrix_sum[n1] # ev[n1] + unit1*(μ11 - μ00)*μ00*h[n1,n1]
    end
    println(" count_n2:",count_n2)

    #println(typeof(h))  # Matrix{Float64}
    #println(size(h))    # (copy_size*4*4, copy_size*4*4)
    #println(ndims(h))   # 2
    #println(size(h)[1])

    # return eigenvalues and eigenvectors
    # eigen(h) get eigenvalues and eigenvectors of h
    return eigen(h) # [1D array of len: copy_size*4*4; 2D array of size copy_size*4*4 x copy_size*4*4]
end

# IR spectra
function ir_spectra_centre(νk::Vector{Float64}, eu::Vector{Vector{Float64}}, com_ol::Matrix{Float64}, Δν, nmols_centre::Int64, large_h_matrix_sum)

    # Unperturbed eigenvalues
    ev::Vector{Float64} = fill(ν0, nmols_centre)
    print(nmols_centre)

    eigenvals, eigenvecs = hstatT_centre(ev, eu, com_ol, large_h_matrix_sum) # is with eu=eu_unit_vector_centre,  com_ol=com_ol_centre to get small simulation box h_matrix 

    σ = eigenvals ./ nmols_centre

    μEpda::Vector{Float64} = zeros(nmols_centre)
    μEsda::Vector{Float64} = zeros(nmols_centre)
    μEp::Vector{Float64}   = zeros(nmols_centre)
    μEs::Vector{Float64}   = zeros(nmols_centre)
    pl::Matrix{Float64}    = zeros(nmols_centre,3)
    for m::Int64 in 1:nmols_centre # loop over eigenvecs
        for i::Int64 in 1:nmols_centre # loop over molecules
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
        for m in 1:nmols_centre
            gp = gssn(ν, eigenvals[m], Δν)
            gs = gssn(ν, eigenvals[m], Δν)
            ipda[iν] += unit2*σ[m]*μEpda[m] * gp
            isda[iν] += unit2*σ[m]*μEsda[m] * gs
            ip[iν] += unit2*σ[m]*μEp[m] * gp
            is[iν] += unit2*σ[m]*μEs[m] * gs
        end
    end

    return ipda, isda, ip, is, eigenvecs
end