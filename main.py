# -*- coding: utf-8 -*-
__author__ = 'Vojtech Vozab'

import networkx as nx
from itertools import izip_longest
from math import sqrt
from operator import add

class TipSet:
    def __init__(self, sizex, sizey, coordlist, maxdist=0.05):
        """after initialization the timeframe is set to 0, the list of frames is filled and match_list contains
        a 1-long list of 1 coordinate for each point in the first frame"""
        self.frame_list = coordlist
        self.match_list = []
        self.timeframe = 0
        self.sizex = sizex
        self.sizey = sizey
        self.maxdist = maxdist
        for coord in self.frame_list[0]:
            if coord != (-1, -1, -1):
                self.match_list.append([coord])

    def get_max_bipartite_graph(self, framedif=1):
        # Initialize lists of existing points and new points, add dummy points into new points list
        old_points = self.get_last_frame(framedif)
        new_points = []
        for coord in self.frame_list[self.timeframe+1]:
            if coord != (-1, -1, -1):
                new_points.append(coord)
        for i in range(len(old_points)):
            new_points.append((-1, -1, -1))

        # Build a bipartite graph, set weights to inverse of euclidean distance
        G = nx.complete_bipartite_graph(len(old_points), len(new_points))
        for i in range(0, len(old_points)):
            for j in range(0, len(new_points)):
                if new_points[j] != (-1, -1, -1):
                    dist = anisotropic_euclid_distance(old_points[i], new_points[j])
                else:
                    dist = self.maxdist - 0.0001
                G.edge[i][j + len(old_points)]['weight'] = dist
                G.node[i]['coordinates'] = old_points[i]
                G.node[j + len(old_points)]['coordinates'] = new_points[j]

        return G

    def get_max_weight_matches(self, graph):
        # Compute the max weight matching and remove redundant edges (existing edges with inverse orientation)
        H = nx.max_weight_matching(graph)
        unique_keys = []
        # delete non-unique edges where x>y = y<x
        for key in H:
            if key not in unique_keys:
                unique_keys.append(H[key])
        for key in unique_keys:
            del H[key]
        matches = {}
        for edge in H:
            matches[graph.node[edge]['coordinates']] = graph.node[H[edge]]['coordinates']
        return matches

    def add_matches(self, matches_dict, framedif=1):
        # add matched pairs and vanishing points into match_list, find and add new points, increase the timeframe
        new_points = []
        for coord in self.frame_list[self.timeframe+1]:
            if coord != (-1, -1, -1):
                new_points.append(coord)
        for match in self.match_list:
            current_coord = self.get_last_location(match, framedif)
            if current_coord != (-1, -1, -1):
                match.append(matches_dict[current_coord])
                try:
                    if matches_dict[current_coord] != (-1, -1, -1):
                        new_points.remove(matches_dict[current_coord])
                except ValueError:
                    print "something broke, coord not in new_points list:", matches_dict[current_coord]
            else:
                # tip not present in this frame
                match.append((-1, -1, -1))
        for new_point in new_points:
            temp_coords = []
            for i in range(self.timeframe+1):
                temp_coords.append((-1, -1, -1))
            temp_coords.append(new_point)
            self.match_list.append(temp_coords)
        self.timeframe +=1

    def get_last_location(self, match, framedif):
        result = (-1, -1, -1)
        for i in range(framedif):
            if match[self.timeframe-i] != (-1, -1, -1):
                return match[self.timeframe-i]
        return result

    def get_last_frame(self, framedif=1):
        coordlist = []
        for coords in self.match_list:
            for i in range(framedif):
                if coords[self.timeframe-i] != (-1, -1, -1):
                    coordlist.append(coords[self.timeframe-i])
                    break
        return coordlist

    def get_pairs_for_frame(self, timeframe):
        pairs = []
        for coord_sequence in self.match_list:
            pairs.append((coord_sequence[timeframe-1], coord_sequence[timeframe]))
        return pairs


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


def build_control_tipset(coord_list, sizex, sizey):
    control = TipSet(sizex, sizey, coord_list)
    control.match_list = []
    for coord in control.frame_list[0]:
        control.match_list.append([coord])
    for frame in range(1, len(control.frame_list)):
        for point in range(len(control.frame_list[frame])):
            control.match_list[point].append(control.frame_list[frame][point])
    control.timeframe = len(control.match_list[0])-1
    return control


def check_found_and_control(found, control, frame, framedif):  # funguje to, ale je potreba divat se o x kroku zpatky... uf.
    found_pairs = found.get_pairs_for_frame(frame)
    control_pairs = control.get_pairs_for_frame(frame)
    correctly_matched = 0
    incorrectly_matched = 0
    correctly_vanished = 0
    incorrectly_vanished = 0
    correctly_found = 0
    incorrectly_found = 0
    for pair in found_pairs:
        if pair[0] != (-1, -1, -1) and pair[1] != (-1, -1, -1):
            if pair in control_pairs:
                correctly_matched += 1
            else:
                incorrectly_matched += 1
        elif pair[0] == (-1, -1, -1) and pair[1] != (-1, -1, -1):  # tady by se melo checkovat o par frejmu zpatky
            if pair in control_pairs:
                correctly_vanished += 1
            else:
                incorrectly_vanished += 1
        elif pair[0] != (-1, -1, -1) and pair[1] == (-1, -1, -1):
            if pair in control_pairs:
                correctly_found += 1
            else:
                incorrectly_found += 1
        else:  # (-1, -1, -1) matched with (-1, -1, -1) = point doesn't exist in this timeframe
            pass

    return [correctly_matched, incorrectly_matched, correctly_vanished, incorrectly_vanished, correctly_found, incorrectly_found]


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


def anisotropic_euclid_distance(point1, point2, zratio=0.125):
    xa, ya, za = point1
    xb, yb, zb = point2
    za = za / zratio
    zb = zb / zratio
    try:
        dist = 1/sqrt((xa-xb)**2 + (ya-yb)**2 + (za-zb)**2)
    except ZeroDivisionError:
        dist = 9999999999                                              # Tohle mozna neni tak uplne nejlepsi napad
        # print "dist is 0, point coords are", point1, point2
    return dist


def convert_match_to_offset(match, sizex, sizey):
    return get_offset_from_xyz(match[0], sizex, sizey), get_offset_from_xyz(match[1], sizex, sizey)


def apply_and_check(filename, sizex, sizey, frames=0):
    points_list = read_points(filename, sizex, sizey)
    tipset = TipSet(sizex, sizey, points_list, 0.005)  # s moc malou hranici pada, proc?
    for i in range(len(points_list)-1):
        G = tipset.get_max_bipartite_graph(5)
        matches = tipset.get_max_weight_matches(G)
        tipset.add_matches(matches, 5)

    # checking section
    results = [0, 0, 0, 0, 0, 0]
    control = build_control_tipset(points_list, sizex, sizey)
    for frame in range(1, len(points_list)):
        tempresults = check_found_and_control(tipset, control, frame)
        results = map(add, results, tempresults)
    print "Checking file", filename
    print "correctly matched:", results[0], "incorrectly matched:", results[1]
    print "correctly recognized as new:", results[4], "incorrectly recognized as new:", results[5]
    print "correctly recognized as vanishing:", results[2], "incorrectly recognized as vanishing:", results[3]
    print "-------------------"
    for point in tipset.match_list:
        print [get_offset_from_xyz(item, sizex, sizey) for item in point]


apply_and_check("./real_images/real-control-01.txt", 220, 220, 0)

