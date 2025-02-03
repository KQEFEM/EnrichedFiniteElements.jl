module BasisFunctions

export phi, p_matrix, e_space, e_time_mass, e_time_ansatz, e_time_test  

function phi(x::Real, y::Real)
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

function p_matrix(x::Real, y::Real)
    """
    Defines the P matrix (Phi * Phi').

    Args:
        x: x-coordinate.
        y: y-coordinate.

    Returns:
        A 3x3 matrix.
    """
    phi_val = phi(x, y)
    return phi_val * phi_val'  # Use * for matrix multiplication
end

function e_space(x::Real, y::Real, A::Real, B::Real, C::Real)
    """
    Defines the spatial enrichment function.

    Args:
        x: x-coordinate.
        y: y-coordinate.
        A: Parameter A.
        B: Parameter B.
        C: Parameter C.

    Returns:
        A complex number.
    """
    return exp(1im * (A * x + B * y + C))  # 1im represents the imaginary unit
end

function e_time_mass(t::Real, w::Real, dt::Real, t0::Real, ww::Real)
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

function e_time_ansatz(t::Real, w::Real)
    """
    Defines the time enrichment function for the ansatz function.

    Args:
        t: Time.
        w: Parameter w.

    Returns:
        A complex number.
    """
    return exp(1im * w * t)
end

function e_time_test(t::Real, ww::Real)
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