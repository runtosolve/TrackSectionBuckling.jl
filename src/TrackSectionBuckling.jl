module TrackSectionBuckling


using CrossSectionGeometry, CUFSM, AISIS100, SectionProperties, OpenSectionBuckling


struct Material

    E
    ν

end

struct Dimensions

    t 

    B   #flange
    H   #web
    r   #inside radius 

end


struct Load
   
    P
    Mxx
    Mzz
    M11
    M22

end

struct Results 

    model 
    Lcr
    Rcr 

end


struct Section

    label
    material
    dimensions
    load
    results

end




function get_section_coordinates(dimensions)

    (;
    t, 

    B,
    H, 
    r  #inside radius 
    ) = dimensions 

    section_dimensions = [B, H, B]
    r = [r+t, r+t]
    n = [3, 3, 3]
    n_r = [3, 3];
    θ = [π, -π/2, 0.0]

    coordinates = CrossSectionGeometry.create_thin_walled_cross_section_geometry(section_dimensions, θ, n, r, n_r, t, centerline = "to left", offset = (section_dimensions[2], section_dimensions[3]))

    X = [coordinates.centerline_node_XY[i][1] for i in eachindex(coordinates.centerline_node_XY)]
    Y = [coordinates.centerline_node_XY[i][2] for i in eachindex(coordinates.centerline_node_XY)]

    coordinates = (X=X, Y=Y)

    return coordinates

end



function get_straight_corner_section_coordinates(dimensions)

    (;
    t, 

    B,
    H, 
    r  #inside radius 
    ) = dimensions 

    section_dimensions = [B, H, B]
    n = [1, 1, 1]
    θ = [π, -π/2, 0.0]

    coordinates = CrossSectionGeometry.create_thin_walled_cross_section_geometry(section_dimensions, θ, n, t, centerline = "to left", offset = (section_dimensions[2], section_dimensions[3]))

    X = [coordinates.centerline_node_XY[i][1] for i in eachindex(coordinates.centerline_node_XY)]
    Y = [coordinates.centerline_node_XY[i][2] for i in eachindex(coordinates.centerline_node_XY)]

    coordinates = (X=X, Y=Y, L=section_dimensions)

    return coordinates

end


function calculate_constrained_buckling_properties(straight_corner_coordinates, dimensions, loads, material, mode_type)

    (;E,
    ν
    ) = material 

    (;P,
    Mxx,
    Mzz,
    M11,
    M22,
    ) = loads

    t = dimensions.t



    ###

    # flat_mesh_size_goal = 0.5
    # corner_mesh_size_goal = π/6
    # centerline_radius = r + t/2


    #coordinates coming in are straight corner X, Y

    coordinates = (straight=straight_corner_coordinates, rounded=straight_corner_coordinates) #rounded not used for now 

    material = (E=E, ν=ν)

    loads = (P=P, Mxx=Mxx, Mzz=Mzz, M11=M11, M22=M22)

    # mode_type = ["D"]
    model = OpenSectionBuckling.properties(coordinates, material, loads, t, mode_type)

    return model 

end


function calculate_section_properties(coordinates, t)

    X = coordinates.X 
    Y = coordinates.Y

    t_all = fill(t, length(X)-1)

    section_properties = SectionProperties.open_thin_walled(X, Y, t_all)

    return section_properties 

end








function calculate_buckling_properties(coordinates, dimensions, loads, material, lengths)

    (;E,
    ν
    ) = material 

    (;P,
    Mxx,
    Mzz,
    M11,
    M22,
    ) = loads

    t = dimensions.t

    constraints = []
    springs = []
    
    neigs = 1

    model = CUFSM.Tools.open_section_analysis(coordinates.X, coordinates.Y, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs)

    eig = 1
    Rcr_curve = CUFSM.Tools.get_load_factor(model, eig)

    # if argmin(Rcr_curve) == length(Rcr_curve)

    Rcr = minimum(Rcr_curve) 
    
    index = argmin(Rcr_curve)
    Lcr = lengths[index]

    results = Results(model, Lcr, Rcr)

    return results 

end

# function calculate_Lcrd(dimensions, material, load_type)

#     (;t, 
#     L, 
#     B,
#     H, 
#     r ) = dimensions 

#     (;E,
#     ν) = material 

#     CorZ = 0
#     θ_top = 90.0
#     #Calculate top flange + lip section properties.
#     CorZ = 0

#     b = B - t
#     d = L - t/2
#     θ = 90.0
#     ho = H
#     μ = ν
#     E = 29500.0
#     G = E / (2 * (1 + μ))
#     kϕ = 0.0
#     Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,xhf,yhf,yof = AISIS100.v16S3.table_2_3_3__1(CorZ,t,b,d,θ)

#     #Calculate the purlin distortional buckling half-wavelength.
#     Lm = 999999999.0

#     if load_type == "P"
#         Lcrd = AISIS100.v16S3.appendix2_2_3_3_1__7(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)
#     elseif load_type == "M"
#         Lcrd = AISIS100.v16S3.appendix2_2_3_3_2__4(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)
#     end

#     return Lcrd

# end



function calculate_Pcrℓ(dimensions, material)

    label = "Pcrℓ"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    #catch indistinct local buckling minimum here 
    mode_type = ["L"]
    
    coordinates = get_straight_corner_section_coordinates(dimensions) #straight corner coordinates 
    model = calculate_constrained_buckling_properties(coordinates, dimensions, load, material, mode_type)

    Lcr_cFSM = model.curve[1][1]

    lengths = range(0.75 * Lcr_cFSM, 1.25 * Lcr_cFSM, 9)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates 

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end




function calculate_Mcrℓ_xx(dimensions, material)

    label = "Mcrℓ_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    mode_type = ["L"]
    
    coordinates = get_straight_corner_section_coordinates(dimensions) #straight corner coordinates 
    model = calculate_constrained_buckling_properties(coordinates, dimensions, load, material, mode_type)

    Lcr_cFSM = model.curve[1][1]

    lengths = range(0.75 * Lcr_cFSM, 1.25 * Lcr_cFSM, 9)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates 

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end




function calculate_Mcrℓ_yy_pos(dimensions, material)

    label = "Mcrℓ_yy_pos"

    load = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    mode_type = ["L"]
    
    coordinates = get_straight_corner_section_coordinates(dimensions) #straight corner coordinates 
    model = calculate_constrained_buckling_properties(coordinates, dimensions, load, material, mode_type)

    Lcr_cFSM = model.curve[1][1]

    lengths = range(0.75 * Lcr_cFSM, 1.25 * Lcr_cFSM, 9)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates 

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Mcrℓ_yy_neg(dimensions, material)

    label = "Mcrℓ_yy_neg"

    load = Load(0.0, 0.0, 1.0, 0.0, 0.0)

    mode_type = ["L"]
    
    coordinates = get_straight_corner_section_coordinates(dimensions) #straight corner coordinates 
    model = calculate_constrained_buckling_properties(coordinates, dimensions, load, material, mode_type)

    Lcr_cFSM = model.curve[1][1]

    lengths = range(0.75 * Lcr_cFSM, 1.25 * Lcr_cFSM, 9)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates 

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)


    return section 

end


function calculate_all_buckling_modes(dimensions, material)


    all_buckling_modes = (Pcrℓ = calculate_Pcrℓ(dimensions, material),
    Mcrℓ_xx = calculate_Mcrℓ_xx(dimensions, material),
    Mcrℓ_yy_pos = calculate_Mcrℓ_yy_pos(dimensions, material),
    Mcrℓ_yy_neg = calculate_Mcrℓ_yy_neg(dimensions, material))


end




end # module TrackSectionBuckling
