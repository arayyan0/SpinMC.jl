struct LSWT
end

function findLocalRotation(spin::Vector{Float64})::Matrix{Float64}
    """
    findLocalRotation(spin::Vector{Float64})::Matrix{Float64}

    given spin = (s1, s2, s3) with s1,s2,s3 real, return matrix R such that
                        (s1, s2, s3).T = R (0, 0, 1).T
    then s3 is the quantization axis with s1, s2 spanning the perpendicular plane
    """
    s1, s2, s3 = spin
    sperp = norm(spin[1:2])
    return sperp ≈ 0 ? I : [s3*s1/sperp s3*s2/sperp -sperp; -s2/sperp s1/sperp 0; s1 s2 s3]' 
end