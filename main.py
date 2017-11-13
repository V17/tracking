# -*- coding: utf-8 -*-
__author__ = 'Vojtech Vozab'

import networkx as nx
import matplotlib.pyplot as plt
from itertools import izip_longest
from math import sqrt


def read_points(filename, sizex, sizey):
    """
    Reads a list of point offsets from a file and transforms it into a 2-dimensional list of coordinates, where first 
    index points to a time frame and second index points to point coordinates present at that time frame
    :param filename: path to a text file containing point offsets
    :return: 2-dimensional list of coordinates
    """
    with open(filename, "r") as f:
        data = f.readlines()
    data_split = []
    for line in data:
        temp = line.split(":")[1].split(",")[:-1]
        data_split.append([int(x) for x in temp])
    # switch axes of multidimensional list of points so that first axis is time:
    data_split = [[i for i in element if i is not None] for element in list(izip_longest(*data_split))]
    coord_list = []
    for frame in data_split:

        coord_list.append([get_xyz_coords(x, sizex, sizey) for x in frame])
    return coord_list


def get_xyz_coords(offset, sizex, sizey):
    """
    Computes 3-dimensional coordinates from offset
    :param sizex: image size on the x axis
    :param sizey: image size on the y axis
    """
    if offset == -1:
        return -1, -1, -1
    z = offset // (sizex * sizey)
    offset = offset - (z * sizex * sizey)
    y = offset // sizex
    x = offset - (y * sizex)
    return x, y, z


def get_offset_from_xyz(xyz, sizex, sizey):
    x, y, z = xyz
    if xyz == (-1, -1, -1):
        return -1
    else:
        return z * sizex * sizey + y * sizex + x


def build_bipartite_graph(points1, points2):
    """
    Builds a bipartite graph from 2 sets of points and writes the euclidean distance into the edge weights and point
    coordinates into the nodes
    :return: Graph object containing the above
    """
    """points1 = []
    points2 = []
    for point in range(len(points1)):
        if points1_full[point] != (-1, -1, -1) or points2_full[point] != (-1, -1, -1):
            points1.append(points1_full[point])
            points2.append(points2_full[point])"""

    points1 = [x for x in points1 if x != (-1, -1, -1)]
    points2 = [x for x in points2 if x != (-1, -1, -1)]
    size1 = len(points1)
    size2 = len(points2)
    size = max(size1, size2)
    if size1 < size2:
        for i in range(size2-size1):
            points1.append((-1, -1, -1))
    if size1 > size2:
        for i in range(size1-size2):
            points2.append((-1, -1, -1))
    print "matching", len(points1)+len(points2), "points"
    print "set 1:", points1
    print "set 2:", points2
    G = nx.complete_bipartite_graph(size, size)

    for i in range(0, size):
        for j in range(0, size):
            if points1[i] != (-1, -1, -1) and points2[j] != (-1, -1, -1):
                dist = anisotropic_euclid_distance(points1[i], points2[j])
            else:
                dist = 0.0000000001
            G.edge[i][j + size]['weight'] = dist
            G.node[i]['coordinates'] = points1[i]
            G.node[j+size]['coordinates'] = points2[j]
    return G


def find_max_match(bipartite_graph):
    """
    Computes the maximum-weight matching of a bipartite graph and returns the result in the form of coordinate pairs
    :return: list of tuples containing matched coordinate pairs
    """
    H = nx.max_weight_matching(bipartite_graph)
    coord_pairs = []
    for edge in H:
        coord_pairs.append((bipartite_graph.node[edge]['coordinates'], bipartite_graph.node[H[edge]]['coordinates']))
    return coord_pairs[:len(coord_pairs)/2]


def anisotropic_euclid_distance(point1, point2, zratio=0.125):
    xa, ya, za = point1
    xb, yb, zb = point2
    za = za / zratio
    zb = zb / zratio
    try:
        dist = 1/sqrt((xa-xb)**2 + (ya-yb)**2 + (za-zb)**2)
    except ZeroDivisionError:
        dist = 9999999999                                              # Tohle mozna neni tak uplne nejlepsi napad
        print "dist is 0, point coords are", point1, point2
    return dist


def check_points(coord_pairs, points, frame):
    """
    Checks the generated matches with the original matches from the loaded point file
    :param points: original point pairs loaded from file
    :param frame: int, indicates which time frame is checked with its next neighbor
    """
    true_matches = 0
    false_matches = 0
    for pair in coord_pairs:
        coord1 = pair[0]
        coord2 = pair[1]
        point_index = points[frame].index(coord1)
        if points[frame+1][point_index] == coord2:
            true_matches += 1
            print "matched point", get_offset_from_xyz(coord1, sizex, sizey), coord1, "with", get_offset_from_xyz(coord2, sizex, sizey), coord2, "correctly"
        else:
            false_matches += 1
            print "matched point", get_offset_from_xyz(coord1, sizex, sizey), coord1, "with", get_offset_from_xyz(coord2, sizex, sizey), coord2, "incorrectly, correct match was", get_offset_from_xyz(points[frame+1][point_index], sizex, sizey), points[frame+1][point_index]
    print "true matches:", true_matches, ", false matches:", false_matches
    print "---"
    return true_matches, false_matches


filename = "./tip_tracking/ID428/tip_trajectories_pxOffsets.txt"
sizex = 314
sizey = 261
points_list = read_points(filename, sizex, sizey)
total_true_matches = 0
total_false_matches = 0
for i in range(len(points_list)-1):
    print "matching frames", i, ",", i+1
    b_g = build_bipartite_graph(points_list[i], points_list[i+1])
    pairs = find_max_match(b_g)
    matches = check_points(pairs, points_list, i)
    total_true_matches += matches[0]
    total_false_matches += matches[1]

print "Matched", total_true_matches, "correctly and", total_false_matches, "incorrectly in total."

