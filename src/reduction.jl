
export get_L0_science, get_L0_flats, get_L0_fiberflats, get_L0_darks, get_L0_badpix_mask

const DEFAULT_REDUCE_OPTIONS = Dict(
    "trace_xrange" => [10, 2038],
    "trace_x0" => 1024,
    "trace_poly_deg" => 6,
    "trace_σ_bounds" => [0.1, 3.0],
    "trace_fit_window" => 14,
    "trace_min_spacing" => 12,
    "trace_flux_smooth_width" => 101,
    "profile_deg_x" => 3,
    "profile_deg_y" => 3,
    "profile_knot_spacing_x" => 600,
    "profile_knot_spacing_y" => 1,
    "extract_xrange" => [20, 2048],
    "extract_yrange" => [-5, 5]
)

function reduce(;
        science::Vector{PARVIL0}=PARVIL0[],
        darks::Union{Vector{PARVIL0}, Nothing}=nothing,
        flats::Union{Vector{PARVIL0}, Nothing}=nothing,
        fiberflats::Union{AbstractDict{Int, Vector{PARVIL0}}, Nothing}=nothing,
        extra::Union{AbstractDict{String, Vector{PARVIL0}}, Nothing}=nothing,
        badpix_mask::Union{PARVIL0, Nothing}=nothing,
        sci_output_path::Union{String, Nothing}=nothing,
        cal_output_path::Union{String, Nothing}=nothing,
        traces::Union{AbstractDict{String, <:NamedTuple}, String, Nothing}=nothing,
        profiles::Union{AbstractDict{String, LSQBivariateSpline}, String, Nothing}=nothing,
        options::Union{Dict{String, <:Any}, Nothing}=nothing
    )

    # Options
    options = isnothing(options) ? deepcopy(DEFAULT_REDUCE_OPTIONS) : merge(DEFAULT_REDUCE_OPTIONS, options)

    # Create the output directories
    create_output_dirs(science, sci_output_path, cal_output_path)

    # Generate calibration images
    median_dark = !isnothing(darks) ? gen_median_dark(darks, cal_output_path) : nothing
    median_flat = !isnothing(flats) ? gen_median_flat(flats, cal_output_path; dark=median_dark) : nothing
    median_fiberflats = !isnothing(fiberflats) ? OrderedDict(fiber => gen_median_fiberflat(fiberflats[fiber], cal_output_path; fiber) for fiber in keys(fiberflats)) : nothing
    median_extra = !isnothing(extra) ? OrderedDict(key => gen_median_image(key => extra[key], cal_output_path) for key in keys(extra)) : nothing
    
    # Trace orders
    if isnothing(traces)
        traces, _ = trace(median_fiberflats, cal_output_path, options; dark=median_dark, flat=median_flat)
    elseif traces isa String
        traces = jldopen(traces) do f
            f["traces"]
        end
    end

    # Profiles from fiber flats
    if isnothing(profiles)
        profiles, _ = get_extraction_profiles(median_fiberflats, traces, cal_output_path, options; dark=median_dark, flat=median_flat)
    elseif profiles isa String
        profiles = jldopen(profiles) do f
            f["profiles"]
        end
    end

    # Extract all images
    extract(
        science,
        traces, profiles,
        sci_output_path, cal_output_path,
        options;
        dark=median_dark, flat=median_flat, fiberflats=median_fiberflats, extra=median_extra, badpix_mask=badpix_mask
    )

end


function create_output_dirs(science::Vector{PARVIL0}, sci_output_path::String, cal_output_path::String)
    mkpath(cal_output_path)
    for d ∈ science
        mkpath(sci_output_path * get_object(d))
    end
end


function EchelleReduce.gen_median_dark(darks::Vector{PARVIL0}, output_path::String)
    dark_images = read.(darks)
    median_dark_image = EchelleReduce.gen_median_dark(dark_images)
    median_dark = PARVIL0(get_median_dark_filename(darks, output_path))
    FITSIO.fitswrite(median_dark.filename, median_dark_image, header=deepcopy(read_header(darks[1])))
    return median_dark
end


function EchelleReduce.gen_median_flat(flats::Vector{PARVIL0}, output_path::String; dark::Union{PARVIL0, Nothing}=nothing)
    flat_images = read.(flats)
    dark_image = !isnothing(dark) ? read(dark) : nothing
    median_flat_image = gen_median_flat(flat_images; dark_image, q=0.5)
    median_flat = PARVIL0(get_median_flat_filename(flats, output_path))
    FITSIO.fitswrite(median_flat.filename, median_flat_image, header=deepcopy(read_header(flats[1])))
    return median_flat
end


function gen_median_fiberflat(fiberflats::Vector{PARVIL0}, output_path::String; fiber::Int)
    fiberflat_images = read.(fiberflats)
    median_fiberflat_image = EchelleReduce.median_combine_images(fiberflat_images)
    median_fiberflat = PARVIL0(get_median_fiberflat_filename(fiberflats, output_path; fiber))
    FITSIO.fitswrite(median_fiberflat.filename, median_fiberflat_image, header=deepcopy(read_header(fiberflats[1])))
    return median_fiberflat
end


function gen_median_image(data::Pair{String, Vector{PARVIL0}}, output_path::String)
    images = read.(data.second)
    median_image = EchelleReduce.median_combine_images(images)
    median_data = PARVIL0(output_path * "median_" * data.first * ".fits")
    FITSIO.fitswrite(median_data.filename, median_image, header=deepcopy(read_header(data.second[1])))
    return median_data
end


function trace(
        fiberflats::AbstractDict{Int, PARVIL0},
        output_path::String,
        options::Dict{String, <:Any};
        dark::Union{PARVIL0, Nothing}=nothing,
        flat::Union{PARVIL0, Nothing}=nothing,
    )

    # Store all traces in an ordered dict
    traces = OrderedDict{String, NamedTuple}()

    # Dark/flat
    dark_image = !isnothing(dark) ? read(dark) : nothing
    flat_image = !isnothing(flat) ? read(flat) : nothing

    # Loop ove fiberflats (fibers)
    ts = ""
    for fiber in keys(fiberflats)
        fiberflat = fiberflats[fiber]
        println("Tracing orders for $fiberflat ...")
        fiberflat_image = read(fiberflat)
        #calibrate_image!(fiberflat_image; dark_image, flat_image)
        fiberflat_image[fiberflat_image .< 0] .= NaN
        #clamp!(fiberflat_image, 0, Inf)
        labels = string.(options["trace_orders"]) .* ".$fiber"
        r = trace_gauss(
            fiberflat_image, labels;
            xrange=options["trace_xrange"], yrange=options["trace_yrange"][fiber],
            x0=options["trace_x0"],
            min_spacing=options["trace_min_spacing"],
            deg=options["trace_poly_deg"],
            σ_bounds=options["trace_σ_bounds"],
            fit_window=options["trace_fit_window"],
            flux_smooth_width=options["trace_flux_smooth_width"]
        )

        # Plot image
        filename_out = "$(output_path)tracepos_$(get_timestamp(fiberflat))_fiber$fiber.png"
        plot_tracepos_image(fiberflat_image, getproperty.(values(r), :yc), filename_out, qmax=0.9)

        # Store
        merge!(traces, r)

        # Increase timestamps for merged trace file
        ts *= get_timestamp(fiberflat) * "_"
    end

    # Store in data
    filename_out = "$(output_path)traces_$(ts[1:end-1]).jld"

    # Save
    jldsave(filename_out; traces)

    # Return
    return traces, filename_out

end


function get_extraction_profiles(
        fiberflats::AbstractDict{Int, PARVIL0},
        traces::AbstractDict{String, <:NamedTuple},
        output_path::String,
        options::Dict{String, <:Any};
        dark::Union{PARVIL0, Nothing}=nothing,
        flat::Union{PARVIL0, Nothing}=nothing,
    )

    # Store all profiles in an ordered dict
    profiles = OrderedDict{String, LSQBivariateSpline}()

    # Dark/flat
    dark_image = !isnothing(dark) ? read(dark) : nothing
    flat_image = !isnothing(flat) ? read(flat) : nothing

    # Timestamps
    ts = ""
    for fiberflat in values(fiberflats)
        ts *= get_timestamp(fiberflat) * "_"
    end

    # Loop over traces
    for fiber in keys(fiberflats)
        ff_image = read(fiberflats[fiber])
        calibrate_image!(ff_image; dark_image, flat_image)
        bad = findall(ff_image .<= 0)
        ff_image[bad] .= NaN
        for order in options["trace_orders"]
            trace = traces["$order.$fiber"]
            println("Determining trace profile: $(fiberflats[fiber]), $(trace.label) ...")
            profile = get_profile_spline(
                ff_image, trace.yc, options["xrange"], trace.yrange;
                knot_spacing_x=options["profile_knot_spacing_x"], knot_spacing_y=options["profile_knot_spacing_y"],
                deg_x=options["profile_deg_x"], deg_y=options["profile_deg_y"]
            )
            profiles[trace.label] = profile
        end
    end

    # Save
    filename_out = "$(output_path)profiles_$(ts[1:end-1]).jld"
    jldsave(filename_out; profiles)

    # Return
    return profiles, filename_out

end


function extract(
        data::PARVIL0,
        traces::AbstractDict{String, <:Any},
        profiles::AbstractDict{String, <:Any},
        options::Dict{String, <:Any};
        flat::Union{PARVIL0, Nothing}=nothing,
        dark::Union{PARVIL0, Nothing}=nothing,
        badpix_mask::Union{PARVIL0, Nothing}=nothing
    )

    # Load image (units are ADU/s)
    data_image = read(data)

    # Load cals
    dark_image = !isnothing(dark) ? read(dark) : nothing
    flat_image = !isnothing(flat) ? read(flat) : nothing

    # Pre calibrate (image units are ADU/s, dark units are ADU/s, flat units are normalized ADU/s)
    calibrate_image!(data_image; dark_image, flat_image)

    # Convert to pe from ADU / s
    itime = get_itime(data)
    data_image .*= DETECTOR_GAIN * itime

    # Bad pix masking
    if !isnothing(badpix_mask)
        badpix_mask_image = 1.0 .- read(badpix_mask)
        bad = findall(.~isfinite.(badpix_mask_image) .|| badpix_mask_image .== 0)
        data_image[bad] .= NaN
    end

    # Which traces to extract
    labels_extract = get_extract_labels(data, traces, options)

    # Extract traces -> Type is Union{NamedTuple, Nothing}
    reduced_data = OrderedDict{String, Union{<:NamedTuple, Nothing}}()
    for trace_label in labels_extract
        extractor = get_extractor(data, traces[trace_label], profiles[trace_label], options)
        reduced_data[trace_label] = extract_trace(data.filename, data_image, trace_label, extractor)
    end

    # Return
    return reduced_data

end


function get_extract_labels(data::PARVIL0, traces::OrderedDict{String, <:NamedTuple}, options::Dict{String, Any})
    labels = collect(keys(traces))
    r = split.(labels, '.')
    orders, fibers = parse.(Int, getindex.(r, 1)), parse.(Int, getindex.(r, 2))
    labels_extract = copy(labels)
    orders_extract = copy(orders)
    fibers_extract = copy(fibers)
    if occursin("fiber1", lowercase(data.filename))
        good = findall(fibers_extract .== 1)
        labels_extract, orders_extract, fibers_extract = labels_extract[good], orders_extract[good], fibers_extract[good]
    elseif occursin("fiber3", lowercase(data.filename))
        good = findall(fibers_extract .== 3)
        labels_extract, orders_extract, fibers_extract = labels_extract[good], orders_extract[good], fibers_extract[good]
    end
    good = collect(Iterators.flatten([findall(orders_extract .== order) for order in options["extract_orders"]]))
    labels_extract = labels_extract[good]
    return labels_extract
end


function get_extractor(data::PARVIL0, trace::NamedTuple, profile::Any, options::Dict{String, <:Any})

    # ron
    read_noise = get_read_noise(data)

    # Optimal for sci
    return OptimalExtractor(;
        yc=trace.yc, spl=profile, xrange=options["xrange"], yrange=options["extract_yrange"], max_iterations=20, read_noise=read_noise, badpix_nσ=8, yrange_extract=options["extract_yrange"],
    )

    # Boxcar for emission ... ?
    # return Boxcar(;
    #     xrange::Vector{Int},
    #     yrange::Vector{Float64},
    #     yc::Vector{Float64},
    #     collapse_function::String,
    #     read_noise::Real=0,
    # )

end


function extract(
        science::Vector{PARVIL0},
        traces::AbstractDict{String, <:Any}, profiles::AbstractDict{String, <:Any},
        sci_output_path::String, cal_output_path::String,
        options::Dict{String, <:Any};
        dark::Union{PARVIL0, Nothing}=nothing,
        flat::Union{PARVIL0, Nothing}=nothing,
        fiberflats::Union{AbstractDict{Int, <:Any}, Nothing}=nothing,
        badpix_mask::Union{PARVIL0, Nothing}=nothing,
        extra::Union{AbstractDict{String, <:Any}, Nothing}=nothing,
    )

    # Science
    pmap(science) do sci_data

        # Extract
        reduced_data = extract(sci_data, traces, profiles, options; dark, flat, badpix_mask)

        # Plot Reduced Spectrum in ADU/s
        object = replace(get_object(sci_data), " " => "_")
        output_path = sci_output_path * object * PATHSEP
        filename_out = output_path * basename(sci_data)[1:end-5] * "_reduced.png"
        itime = get_itime(sci_data)
        plot_reduced_spectrum(reduced_data, filename_out; itime, gain=DETECTOR_GAIN)

        # Save Reduced Spectrum in pe
        save_extraction_results(sci_data, output_path, reduced_data)

    end

    # Fiber flats
    pmap(values(fiberflats)) do fiberflat
        
        # Extract all orders and fibers
        reduced_data = extract(fiberflat, traces, profiles, options; dark, flat, badpix_mask)

        # Plot Reduced Spectrum
        filename_out = cal_output_path * basename(fiberflat)[1:end-5] * "_reduced.png"
        itime = get_itime(fiberflat)
        plot_reduced_spectrum(reduced_data, filename_out; itime, gain=DETECTOR_GAIN)

        # Save Reduced Spectrum
        save_extraction_results(fiberflat, cal_output_path, reduced_data)

    end

    # Extract extra
    pmap(collect(keys(extra))) do key

        # Data
        d = extra[key]
        
        # Extract all orders and fibers
        reduced_data = extract(d, traces, profiles, options; dark, flat, badpix_mask)

        # Plot Reduced Spectrum
        filename_out = cal_output_path * basename(d)[1:end-5] * "_reduced.png"
        itime = get_itime(d)
        plot_reduced_spectrum(reduced_data, filename_out; itime, gain=DETECTOR_GAIN)

        # Save Reduced Spectrum
        save_extraction_results(d, cal_output_path, reduced_data)

    end

end

########################


function get_science_key(sci_data::PARVIL0, data::Dict{String, <:Any})
    for key in keys(data["science"])
        if sci_data in data["science"][key]
            return key
        end
    end
end

########################


function save_extraction_results(data::PARVIL0, output_path::String, reduced_data::OrderedDict{String, <:Any})

    # Filename out
    filename_out = output_path * basename(data)[1:end-5] * "_reduced.fits"

    # New fits file
    f = FITS(filename_out, "w")

    # Copy header
    header = deepcopy(read_header(data))

    # Reduced data
    reduced_data_out = OrderedDict{String, Any}()
    for trace_label ∈ keys(reduced_data)
        if !isnothing(reduced_data[trace_label])
            reduced_data_out[trace_label] = [reduced_data[trace_label].spec reduced_data[trace_label].specerr]
        else
            reduced_data_out[trace_label] = fill(NaN, (2048, 2))
        end
    end

    # Write and close
    FITSIO.write(f, reduced_data_out, header=header)
    FITSIO.close(f)

    # Save auxiliary results
    filename_out = output_path * basename(data)[1:end-5] * "_auxiliary.jld"
    
    # Save auxiliary data
    save_auxiliary_data(filename_out, reduced_data)

end

########################################################################

function get_read_noise(data::AnyPARVI)
    header = read_header(data)
    n_reads = header["NGROUPS"]
    ron = DETECTOR_READ_NOISE[n_reads]
    return ron
end

########################################################################

function get_median_fiberflat_filename(fiberflats::Vector{PARVIL0}, output_path::String; fiber::Int)
    ts = get_timestamp(fiberflats[1])
    fname = "$(output_path)median_fiberflat_fiber$(fiber)_$(ts).fits"
    return fname
end


function get_median_flat_filename(flats::Vector{PARVIL0}, output_path::String)
    ts = get_timestamp(flats[1])
    filename = "$(output_path)median_flat_$(ts).fits"
    return filename
end


function get_median_image_filename(data::Vector{PARVIL0}, output_path::String; object::String)
    ts = get_timestamp(data[1])
    filename = "$(output_path)median_$(object)_$(ts).fits"
    return filename
end


function get_median_dark_filename(darks::Vector{PARVIL0}, output_path::String)
    ts = get_timestamp(darks[1])
    itime = get_itime(darks[1])
    filename = "$(output_path)median_dark_$(itime)_$(ts).fits"
    return filename
end

########################################################################

function get_L0_science(path::String, objects::Vector{String}, utdates::Vector{String})
    files = PARVIL0[]
    for o ∈ objects
        for u ∈ utdates
            files = vcat(files, PARVIL0.(glob("*$(o)*_*$(u)*.fits", path * o)))
        end
    end
    return files
end


function get_L0_flats(path::String, timestamps::Vector{String})
    files = PARVIL0[]
    for ts ∈ timestamps
        files = vcat(files, PARVIL0.(glob("*$(ts)*.fits", path)))
    end
    return files
end

function get_L0_darks(path::String, timestamps::Vector{String})
    files = PARVIL0[]
    for ts ∈ timestamps
        files = vcat(files, PARVIL0.(glob("*$(ts)*.fits", path)))
    end
    return files
end


function get_L0_fiberflats(path::String, fibers...)
    fiberflat_files = Dict{Int, Vector{PARVIL0}}()
    for p in fibers
        files = PARVIL0[]
        for ts ∈ p.second
            files = vcat(files, PARVIL0.(glob("*$(ts)*.fits", path)))
        end
        fiberflat_files[p.first] = files
    end
    return fiberflat_files
end