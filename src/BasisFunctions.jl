module BasisFunctions

export all

function phi(x::Real, y::Real,z::Real=0.0)
    """
    Defines the Phi basis functions for a linear triangle.

    Args:
        x: x-coordinate.
        y: y-coordinate.

    Returns:
        A 3x1 matrix (column vector) of the basis function values.
    """
    return [1 - x - y; x; y]
end

function p_matrix(x::Real, y::Real,z::Real=0)
    """
    Defines the P matrix (Phi * Phi').

    Args:
        x: x-coordinate.
        y: y-coordinate.

    Returns:
        A 3x3 matrix.
    """
    phi_val = phi(x, y)
    println(phi_val)
    return phi_val * phi_val'  # Use * for matrix multiplication
end

function grads_matrix(grad_matrix::AbstractMatrix{<:Real}, x::Real=0, y::Real=0, z::Real=0)
    """
    Defines the grads matrix (3x3) from a given matrix.

    Args:
       grad_matrix: 3x3 matrix

    Returns:
        A 3x3 matrix.
    """
    rows, cols = size(grad_matrix)
    if rows != 3 || cols != 3
        error("Input matrix must be 3x3.")
    end
    return grad_matrix
end

function enrich_space(x::Real, y::Real, t::Real, A::Real, B::Real, C::Real)
    """
    Defines the spatial enrichment function. This includes the rearrangment for G * conj(G*)

    Args:
        x: x-coordinate.
        y: y-coordinate.
        A: Parameter A: see matlab code.
        B: Parameter B: : see matlab code.
        C: Parameter C.

    Returns:
        A complex number.
    """
    return exp(1im * (A * x + B * y + C))  # 1im represents the imaginary unit
end

function e_time_mass(x::Real,y::Real,t::Real, w::Vector{Float64}, dt::Real, t0::Real, ww::Vector{Float64})
    """
    Defines the time enrichment function for the mass term.

    Args:
        t: Time.
        w: Parameter w.
        dt: Time step.
        t0: Initial time.
        ww: Parameter ww.

    Returns:
        A complex number.
    """
    return exp(1im * w * dt) * exp(-1im * ww * (t - t0))
end

function enrichment_time(x::Real, y::Real, t::Real, w::Float64)
    """
    Defines the time enrichment function for the ansatz function.

    Args:
        t: Time.
        w: Parameter w. Includes the w - ww

    Returns:
        A complex number.
    """
    return exp(1im * w * t)
end

function e_time_test(x::Real, y::Real, t::Real, ww::Vector{Float64})
    """
    Defines the time enrichment function for the test function.

    Args:
        t: Time.
        ww: Parameter ww.

    Returns:
        A complex number.
    """
    return exp(-1im * ww * t)
end

end 