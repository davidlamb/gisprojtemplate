import math

def geodesic(lat1:float,lon1:float,lat2:float,lon2:float)->float:
    """
    Geodesic distance
    Args:
        lat1 (float)
        lon1 (float)
        lat2 (float)
        lon2 (float)

    returns
        float: distance in meters.
    """
    import math
    R = 20902230.97
    phi1 = lat1 * math.pi / 180
    phi2 = lat2 * math.pi / 180
    deltaphi = (lat2 - lat1) * math.pi / 180
    deltalambda = (lon2 - lon1) * math.pi / 180
    a = math.sin(deltaphi / 2) * math.sin(deltaphi/2) + math.cos(phi1) * math.cos(phi2) * math.sin(deltalambda/2) * math.sin(deltalambda/2)
    c = 2 * math.atan2(math.sqrt(a),math.sqrt(1-a))
    d = R * c
    return d

geodesic(41.259444,-110.999721,41.2735,-111.029)