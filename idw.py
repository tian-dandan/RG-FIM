from scipy.spatial import cKDTree
import numpy as np

def idw_interpolation(x, y, kdTree, known_values, numNeighbors=12, power=2):
    """
    Perform IDW interpolation for a given point.
    
    :param x: x coordinate of the point to interpolate
    :param y: y coordinate of the point to interpolate
    :param known_points: array of known points coordinates
    :param known_values: array of known elevation values
    :param power: power parameter for IDW
    :return: interpolated elevation value
    """
    if len(known_values) == 0:
        return np.nan  # Return NaN if no known points are available

    dist, idx = kdTree.query([x, y], k=numNeighbors)
    if np.any(dist == 0):
        # If the point coincides with a known point, return the known elevation
        return known_values[idx[dist == 0][0]]
    
    weights = 1 / dist**power
    weights /= weights.sum()
    interpolated_value = np.dot(weights, known_values[idx])
    return interpolated_value