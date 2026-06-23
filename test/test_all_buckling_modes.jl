using TrackSectionBuckling


material = TrackSectionBuckling.Material(29500.0, 0.3)

dimensions = TrackSectionBuckling.Dimensions(0.0346, 2.0, 3.625, 0.0765)


all_buckling_modes = TrackSectionBuckling.calculate_all_buckling_modes(dimensions, material)

