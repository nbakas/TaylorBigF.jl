module TaylorBigF

    include("calc_Vandermonde.jl")
    include("calc_Vandermonde_inv_formula.jl")
    include("calc_Vandermonde_dot.jl")
    include("calc_Vandermonde_dot_dot.jl")
    include("calc_Vandermonde_dot_dot_dot.jl")
    include("calc_Vandermonde_dot_dot_dot_dot.jl")
    include("calc_Vandermonde_int.jl")
    include("calc_Vandermonde_int_int.jl")
    include("calc_Vandermonde_2D.jl")
    include("calc_Vandermonde_xxxx_2D.jl")
    include("calc_Vandermonde_yyyy_2D.jl")
    include("calc_Vandermonde_xxyy_2D.jl")
    include("calc_Vandermonde_xxx_2D.jl")
    include("calc_Vandermonde_yyy_2D.jl")
    include("calc_Vandermonde_xyy_2D.jl")
    include("calc_Vandermonde_yxx_2D.jl")

    export calc_Vandermonde,calc_Vandermonde_inv_formula,calc_Vandermonde_dot,
    calc_Vandermonde_dot_dot,calc_Vandermonde_dot_dot_dot,calc_Vandermonde_dot_dot_dot_dot,
    calc_Vandermonde_int,calc_Vandermonde_int_int,calc_Vandermonde_2D,
    calc_Vandermonde_xxxx_2D,calc_Vandermonde_yyyy_2D,calc_Vandermonde_xxyy_2D,
    calc_Vandermonde_xxx_2D,calc_Vandermonde_yyy_2D,calc_Vandermonde_xyy_2D,calc_Vandermonde_yxx_2D



end