
# Read headers


# Read headers
Echelle.read_header(data::PARVIL0) = read_header(data.filename, 1)
Echelle.read_key(data::PARVIL0, key::Union{Int, String}) = read_key(data.filename, key, hdu=1)

function Echelle.read_key(data::PARVIL1, key::Union{Int, String})
    if occursin("deg0", data.filename)
        return read_key(data.filename, key, hdu=1)
    else
        return read_key(data.filename, key, hdu=2)
    end
end

function Echelle.read_header(data::PARVIL1)
    if occursin("deg0", data.filename)
        return read_header(data.filename, 1)
    else
        return read_header(data.filename, 2)
    end
end


# Header parsing
Echelle.get_itime(data::AnyPARVI) = read_key(data, "EXPTIME")
Echelle.get_object(data::AnyPARVI) = read_key(data, "OBJECT")
Echelle.get_timestamp(data::AnyPARVI) = data.filename[end-18:end-5]
function Echelle.get_utdate(data::AnyPARVI)
    t_start = get_exposure_start_time(data)
    t_start = julian2datetime(t_start)
    y, m, d = year(t_start), month(t_start), day(t_start)
    if m < 10
        m = "0$m"
    end
    if d < 10
        d = "0$d"
    end
    utdate = join((y, m, d))
    return utdate
end

function Echelle.get_sky_coord(data::AnyPARVI)
    header = read_header(data)
    if !isnothing(header["P200RA"]) && !isnothing(header["P200DEC"])
        coord = ICRSCoords(hms2rad(header["P200RA"]), dms2rad(header["P200DEC"]))
    elseif !isnothing(header["RA"]) && !isnothing(header["DEC"])
        coord = ICRSCoords(hms2rad(header["RA"]), dms2rad(header["DEC"]))
    else
        coord = ICRSCoords(NaN, NaN)
    end
    return coord
end

Echelle.get_exposure_start_time(data::AnyPARVI) = parse(Float64, read_key(data, "TIMEI00")) / 1E9 / 86400 + 2440587.5

# Barycenter
Echelle.get_barycentric_corrections(data::AnyPARVI; star_name::String) = Echelle.get_barycentric_corrections(get_exposure_midpoint_time(data); star_name, obs_name=OBSERVATORY, zmeas=0)

# Read in image
function Base.read(data::PARVIL0; hdu::Int=1, trans=false, mask_edges=true)
    image = FITS(data.filename) do f
        Float64.(read(f[hdu]))
    end
    if trans
        image .= collect(transpose(image))
    end
    if mask_edges
        image[1:4, :] .= NaN
        image[2018:end, :] .= NaN
        image[:, 1:4] .= NaN
        image[:, end-3:end] .= NaN
    end
    return image
end

# Read in reduced spectrum and error, wavelength separate
function Base.read(data::PARVIL1; order::Int, fiber::Int=1, xrange::Union{Vector{Int}, Nothing}=nothing, norm::Union{Float64, Nothing}=nothing)
    if occursin("deg0", data.filename)
        oi = 128 - order
        spec, specerr = FITSIO.FITS(data.filename) do f
            d = read(f[fiber + 1])
            return d[:, oi, 1], d[:, oi, 2], d[:, oi, 3]
        end
    else
        spec, specerr = FITS(data.filename) do f
            d = read(f[2], "$order.$fiber")
            return d[:, 1], d[:, 2]
        end
    end
    if !isnothing(xrange)
        mask!(spec, xrange)
        mask!(specerr, xrange)
    end
    if !isnothing(norm)
        v = nanquantile(quantile_filter1d(spec, width=5), norm)
        spec ./= v
        specerr ./= v
    end
    return spec, specerr
end

