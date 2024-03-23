function inject_λsolution!(
        fname_sci::String, fname_λ::String;
        overwrite::Bool=true
    )

    # Load the wls
    λsols = jldopen(fname_λ) do f
        f["λsols"]
    end

    # Open the sci file
    f = FITSIO.FITS(fname_sci, "r+")

    # Write wls to new HDU
    if overwrite
        while length(f) > 2
            deleteat!(f, length(f))
        end
    end

    # Write the wls
    write(f, λsols)

    # Close
    close(f)

end

function save_λsolution_results(fname::String, result1=nothing, result3=nothing)
    if !isnothing(result1) && !isnothing(result3)
        λsols = merge(result1.λsols, result3.λsols)
    elseif !isnothing(result1)
        λsols = result1.λsols
    else
        λsols = result3.λsols
    end 
    jldsave(fname; λsols, result1, result3)
end

function get_λsolution_from_jld(fname::String, label; pixel_offset=0)
    order, fiber = split(label, ".")
    order = parse(Int, order)
    f = jldopen(fname)
    r = f["result$fiber"]
    close(f)
    xarr = [1:2048;]
    λ = build_λsolution_chebyval2d(xarr .+ pixel_offset, order, r.max_pixel, r.max_order, r.coeffs)[1]
    return λ
end

function get_lfcλsolution_poly2d(
        filename_lfc::String, λ_estimates::AbstractDict{String, <:Any};
        do_orders::AbstractVector{Int}, fiber::Int,
        xrange::Vector{Int}=[20, 2048-20],
        σ_bounds::Vector{<:Real}=[0.2, 4],
        deg_intra_order::Int, deg_inter_order::Int,
        nσ::Real=4,
        fit_window::Union{Vector{<:Real}, Nothing}=nothing,
        orders_full::AbstractVector{Int}=copy(do_orders),
        background_poly_deg::Int=1,
    )

    # Load in the fits file and alias HDU
    file = FITSIO.FITS(filename_lfc)

    # Storage arrays
    pixels = Vector{Float64}[]
    pixels_err = Vector{Float64}[]
    amplitudes = Vector{Float64}[]
    amplitudes_err = Vector{Float64}[]
    integers = Vector{Int}[]
    λs = Vector{Float64}[]
    weights = Vector{Float64}[]
    σs = Vector{Float64}[]
    σs_err = Vector{Float64}[]
    backgrounds = Vector{Float64}[]
    orders = Vector{Float64}[]
    rms = Vector{Float64}[]

    # Data grid
    nx = 2048

    # Loop over orders and get peaks
    for (i, order) ∈ enumerate(do_orders)

        # Label
        label = "$order.$fiber"

        println("Getting LFC peaks for $label")
        
        # Initial wavelength solution
        λ_estimate = λ_estimates[label]
        
        # LFC spec for this order
        lfc_spec = read(file[2], label)[:, 2]
        
        # Get 1D peaks for this order
        result = get_lfc_modes(
            λ_estimate, lfc_spec,  LFC_ν0, LFC_Δν;
            xrange, σ_bounds, background_poly_deg, fit_window
        )

        # Store results
        push!(pixels, result.pixels)
        push!(pixels_err, result.pixels_err)
        push!(orders, fill(order, length(result.pixels)))
        push!(λs, result.λs)
        push!(integers, result.integers)
        push!(amplitudes, result.amplitudes)
        push!(amplitudes_err, result.amplitudes_err)
        #push!(weights, 1 ./ result.pixels_err.^2)
        push!(weights, (result.amplitudes ./ result.rms).^2)
        push!(σs, result.σs)
        push!(backgrounds, result.backgrounds)
        push!(σs_err, result.σs_err)
        push!(rms, result.rms)

    end

    # Close fits file
    close(file)

    # Flatten
    pixels = collect(Iterators.flatten(pixels))
    orders = collect(Iterators.flatten(orders))
    λs = collect(Iterators.flatten(λs))
    weights = collect(Iterators.flatten(weights))
    integers = collect(Iterators.flatten(integers))
    amplitudes_err = collect(Iterators.flatten(amplitudes_err))
    amplitudes = collect(Iterators.flatten(amplitudes))
    backgrounds = collect(Iterators.flatten(backgrounds))
    σs = collect(Iterators.flatten(σs))
    σs_err = collect(Iterators.flatten(σs_err))
    rms = collect(Iterators.flatten(rms))

    # Fit
    println("Fitting LFC modes with 2D Polynomial (deg_intra_order=$deg_intra_order, deg_inter_order=$deg_inter_order)")
    result = EchelleCalibration.fit_peaks_poly2d(
        pixels, orders, λs, weights;
        deg_intra_order, deg_inter_order, nσ
    )

    # Build wls for all orders on full pixel grid
    _λsols = EchelleCalibration.build_λsol2d_poly_full(result.coeffs, 1:2048, orders_full, result.min_pixel, result.max_pixel, result.min_order, result.max_order)
    λsols = OrderedDict{String, Vector{Float64}}("$order.$fiber" => _λsols[i] for (i, order) in enumerate(orders_full))

    # Pack results
    out = (;result..., pixels, pixels_err, orders, λs, weights, integers, amplitudes, amplitudes_err, σs, σs_err, rms, λsols, backgrounds)

    # Return
    return out

end

# LFC info
# 192.1852839999997
const LFC_ν0 = SPEED_OF_LIGHT_MPS / (1559.91370 * 1E-9) # freq of pump line in Hz.
const LFC_Δν = 10.000000174E9 # spacing of peaks [Hz]