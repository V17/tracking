# -*- coding: utf-8 -*-
__author__ = 'Vojtech Vozab'

import networkx as nx
from itertools import izip_longest
from math import sqrt


class TipCoords:
    def __init__(self, coords, timeframe):
        self.coords_list = []
        for i in range(timeframe):
            self.coords_list.append((-1, -1, -1))
        self.coords_list.append(coords)

    def add_coordinate(self, coords):
        self.coords_list.append(coords)

    def get_last_coord(self, framedif):
        if len(self.coords_list) <= framedif:
            framedif = len(self.coords_list)
        for i in range(framedif):
            if self.coords_list[len(self.coords_list)-i-1] != (-1, -1, -1):
                return self.coords_list[len(self.coords_list)-i-1]
        return (-1, -1, -1)
        
    def get_coords(self):
        return self.coords_list

    def get_match(self, timeframe):
        if timeframe > len(self.coords_list)-1:
            print "match not yet made"
        else:
            return (self.coords_list[timeframe-1], self.coords_list[timeframe])


class TipSet:
    def __init__(self):
        self.tip_list = []
        self.timeframe = 0

    def add_new_tip(self, coords):
        self.tip_list.append(TipCoords(coords, self.timeframe))

    def set_timeframe(self, timeframe):
        self.timeframe = timeframe

    def get_timeframe(self):
        return self.timeframe

    def increment_timeframe(self):
        self.timeframe += 1
        for tip in self.tip_list:

            if len(tip.get_coords()) < (self.timeframe):
                tip.add_coordinate((-1, -1, -1))

    def add_coords_to_tip(self, old_coords, new_coords, framedif=1):
        for tip in self.tip_list:
            if old_coords == tip.get_last_coord(framedif):
                tip.add_coordinate(new_coords)
                return 1
        print "error: tip that should have existed not found"
        return 0

    def get_last_frame(self, framedif=1):
        coordlist = []
        for tip in self.tip_list:
            coord = tip.get_last_coord(framedif)
            if coord is not None:
                coordlist.append(tip.get_last_coord(framedif))
        return coordlist

    def get_all_points(self):
        for tip in self.tip_list:
            offsets = []
            for coord in tip.get_coords():
                offsets.append(get_offset_from_xyz(coord, sizex, sizey))
            print offsets

    def get_matches_for_frame(self, timeframe):
        matches = []
        if timeframe > self.timeframe or self.timeframe < 2:
            print "no matches made yet for timeframe", timeframe, "current timeframe:", self.timeframe
            return []
        for tip in self.tip_list:
            matches.append(tip.get_match(timeframe))
        #print "successfuly returning matches for timeframes", timeframe, "and", timeframe-1
        return matches



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


def build_bipartite_graph2(matched_points, new_points, maxdist=0.05, framedif=1):
    old_points = matched_points.get_last_frame(framedif)
    new_points = [x for x in new_points if x != (-1, -1, -1)]
    size1 = len(old_points)
    size2 = len(new_points)
    if size1 < size2:
        for i in range(size2-size1):
            old_points.append((-1, -1, -1))
    for i in range(size1):
        new_points.append((-1, -1, -1))
    G = nx.complete_bipartite_graph(len(old_points), len(new_points))

    for i in range(0, len(old_points)):
        for j in range(0, len(new_points)):
            if new_points[j] != (-1, -1, -1) and old_points[i] != (-1, -1, -1):
                dist = anisotropic_euclid_distance(old_points[i], new_points[j])
            elif new_points[j] == (-1, -1, -1) and old_points[i] == (-1, -1, -1):
                dist = 0.00000
            else:
                dist = maxdist
            G.edge[i][j + len(old_points)]['weight'] = dist
            G.node[i]['coordinates'] = old_points[i]
            G.node[j + len(old_points)]['coordinates'] = new_points[j]
    return G

    
def build_bipartite_graph(points1, points2, maxdist=0.05):
    """
    Builds a bipartite graph from 2 sets of points and writes the euclidean distance into the edge weights and point
    coordinates into the nodes
    :return: Graph object containing the above
    """

    points1 = [x for x in points1 if x != (-1, -1, -1)]
    points2 = [x for x in points2 if x != (-1, -1, -1)]
    size1 = len(points1)
    for i in range(len(points1)):
        points2.append((-1, -1, -1))
    size2 = len(points2)

    print "matching", len(points1)+len(points2), "points"
    print "set 1:", points1
    print "set 2:", points2
    G = nx.complete_bipartite_graph(size1, size2)

    for i in range(0, size1):
        for j in range(0, size2):
            if points2[j] != (-1, -1, -1):
                dist = anisotropic_euclid_distance(points1[i], points2[j])
            else:
                dist = maxdist
            G.edge[i][j + size1]['weight'] = dist
            G.node[i]['coordinates'] = points1[i]
            G.node[j+size1]['coordinates'] = points2[j]
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


def add_matched_pairs(coord_pairs, matched_points):
    matched_points.increment_timeframe()
    for pair in coord_pairs:
        coord1 = pair[0]
        coord2 = pair[1]
        if coord1 != (-1, -1, -1):
            if matched_points.add_coords_to_tip(coord1, coord2) == 1:
                continue
            else:
                print "error with coords", coord1, coord2

        elif coord2 != (-1, -1, -1):
            matched_points.add_new_tip(coord2)
            #print "added new tip in coordinates", coord2

        else:
            continue
            #print "match not added", coord1, coord2



def check_points(coord_pairs, points, frame):
    """
    Checks the generated matches with the original matches from the loaded point file
    :param points: original point pairs loaded from file
    :param frame: int, indicates which time frame is checked with its next neighbor
    """
    true_matches = 0
    false_matches = 0
    unmatched = 0
    for pair in coord_pairs:
        coord1 = pair[0]
        coord2 = pair[1]
        point_index = points[frame].index(coord1)
        if coord1 == (-1, -1, -1) or coord2 == (-1, -1, -1):
            unmatched += 1
        else:
            if points[frame+1][point_index] == coord2:
                true_matches += 1
                print "matched point", get_offset_from_xyz(coord1, sizex, sizey), coord1, "with", get_offset_from_xyz(coord2, sizex, sizey), coord2, "correctly"
            else:
                false_matches += 1
                print "matched point", get_offset_from_xyz(coord1, sizex, sizey), coord1, "with", get_offset_from_xyz(coord2, sizex, sizey), coord2, "incorrectly, correct match was", get_offset_from_xyz(points[frame+1][point_index], sizex, sizey), points[frame+1][point_index]
    print "true matches:", true_matches, ", false matches:", false_matches, "unmatched:", unmatched
    print "---"
    return true_matches, false_matches, unmatched


def check_correct_matches(detected_tipset, control_tipset):
    correct_matches = 0
    wrong_matches = 0
    found = 0
    vanished = 0
    for i in range(1, len(detected_tipset.tip_list[0].coords_list)):
        detected = detected_tipset.get_matches_for_frame(i)
        control = control_tipset.get_matches_for_frame(i)
        #print "detected:", detected
        #print "control:", control
        for match1 in detected:
            if match1[0] != (-1, -1, -1) and match1[1] != (-1, -1, -1):
                for match2 in control:
                    if match1[0] == match2[0] and match1[0] != (-1, -1, -1):
                        if match1[1] == match2[1]:
                            print "correctly matched", match1
                            correct_matches += 1
                        else:
                            print "wrongly matched", match1, "should be", match2
                            wrong_matches += 1
            elif match1[0] == (-1, -1, -1) and match1[1] != (-1, -1, -1):
                print "new tip found at", match1[1]
                found += 1
            elif match1[0] != (-1, -1, -1) and match1[1] == (-1, -1, -1):
                print "tip vanished from", match1[0]
                vanished += 1
            else:
                continue #(-1, -1, -1) matched to (-1, -1, -1)
        print "-----"
    print "Total score:", correct_matches, "matched correctly,", wrong_matches, "matched incorrectly,", found, "new tips,", vanished, "vanished tips"


def build_control_tipset(control_points):
    control_tipset = TipSet()
    for coord in range(len(control_points[0])):
        temptip = TipCoords(control_points[0][coord], 0)
        for frame in range(1, len(control_points)):
            temptip.add_coordinate(control_points[frame][coord])
        control_tipset.tip_list.append(temptip)
    control_tipset.set_timeframe(len(control_points[0]))
    return control_tipset

filename = "./tip_tracking/ID426/tip_trajectories_pxOffsets.txt"
sizex = 295
sizey = 317
points_list = read_points(filename, sizex, sizey)
total_true_matches = 0
total_false_matches = 0
total_unmatched = 0
tipset = TipSet()
for coordinate in points_list[0]:
    if coordinate != (-1, -1, -1):
        tipset.add_new_tip(coordinate)


control = build_control_tipset(points_list)

for i in range(len(points_list[0])-1):
    b_g = build_bipartite_graph2(tipset, points_list[i+1])
    pairs = find_max_match(b_g)
    add_matched_pairs(pairs, tipset)
    #tipset.get_all_points()
    print "-----"
check_correct_matches(tipset, control)



