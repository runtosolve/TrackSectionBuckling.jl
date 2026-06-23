using TrackSectionBuckling


material = TrackSectionBuckling.Material(29500.0, 0.3)

dimensions = TrackSectionBuckling.Dimensions(0.0346, 2.0, 3.625, 0.0765)


section = TrackSectionBuckling.calculate_Pcrℓ(dimensions, material)

