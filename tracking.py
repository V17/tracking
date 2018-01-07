# -*- coding: utf-8 -*-
__author__ = 'Vojtech Vozab'

import networkx as nx
from itertools import izip_longest
from math import sqrt
from operator import add
import time


class TipSet:
    def __init__(self, sizex, sizey, coordlist, maxdist=0.05):
        """
        After initialization the timeframe is set to 0, the list of frames is filled and match_list contains
        a 1-long list of 1 coordinate for each point in the first frame.
        """
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
        """
        Builds a bipartite graph with inverse of euclidian distance as the edge weights, first set contains points in
        the current frame, second set contains points in the next frame plus a dummy point for each point in the first
        frame. Weights of edges leading to dummy points are set to nearly maximum distance.
        :param framedif: How many timeframes back do we look for potential matches, timeframe=1 means only looking at
        the current one.
        """
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
        """
        Computes the maximum weight matching and removes redundant edges (existing edges with inverse orientation).
        """
        H = nx.max_weight_matching(graph)
        for key in H.keys():
            if key > H[key]:
                del H[key]
        matches = {}
        for edge in H:
            matches[graph.node[edge]['coordinates']] = graph.node[H[edge]]['coordinates']
        return matches

    def add_matches(self, matches_dict, framedif=1):
        """
        Adds matched pairs and vanishing points into match_list, finds and adds new points, increases the timeframe.
        :param matches_dict: Output from get_max_weight_matches method
        """
        new_points = []
        for coord in self.frame_list[self.timeframe+1]:
            if coord != (-1, -1, -1):
                new_points.append(coord)
        for match in self.match_list:
            current_coord = self.get_last_location(match, framedif)
            if current_coord != (-1, -1, -1):
                try:
                    match.append(matches_dict[current_coord])
                except KeyError:
                    print "looking for", current_coord, "in", matches_dict, "failed"
                try:
                    if matches_dict[current_coord] != (-1, -1, -1):
                        new_points.remove(matches_dict[current_coord])
                except ValueError:
                    print "Something broke, coord not in new_points list:", get_offset_from_xyz(matches_dict[current_coord], self.sizex, self.sizey)
                    print "This may happen if two points occupy the same coordinates in a single timeframe."
            else:
                # tip not present in this frame
                match.append((-1, -1, -1))
        for new_point in new_points:
            temp_coords = []
            for i in range(self.timeframe+1):
                temp_coords.append((-1, -1, -1))
            temp_coords.append(new_point)
            self.match_list.append(temp_coords)
        self.timeframe += 1

    def get_last_location(self, match, framedif, index=-1):
        """
        Returns the last known coordinate of a tip. If its coordinate on current timeframe is not (-1, -1, -1), it is
        returned, otherwise the method looks up to framedif-1 frames back.
        :param match: Sequence of coordinates
        :param index: If set, it is used as a starting time index instead of timeframe.
        """
        result = (-1, -1, -1)
        if index == -1:
            index = self.timeframe
        if framedif == -1:
            framedif = index + 1
        for i in range(0, framedif):
            if i <= index and match[index-i] != (-1, -1, -1):
                return match[index-i]
        return result

    def get_last_frame(self, framedif=1):
        """
        Returns a list of all coordinates present either in the current timeframe or looking up to framedif-1 frames
        back.
        """
        coordlist = []
        for coords in self.match_list:
            for i in range(framedif):
                if coords[self.timeframe-i] != (-1, -1, -1):
                    coordlist.append(coords[self.timeframe-i])
                    break
        return coordlist

    def get_complete_frame(self, frame):
        return [self.match_list[i][frame] for i in range(len(self.match_list))]

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


def apply_and_write(filename, sizex, sizey, points_list, xypxsize=125, zratio=0.125, framedif=1, maxdistnm=6250):
    if len(points_list) == 0:
        with open(filename + "-" + str(maxdistnm) + "nm-"+str(int(framedif))+"frames.log", "w") as f:
            f.write("empty list")
        return 0
    maxdistweight = 1.0/(float(maxdistnm)/float(xypxsize))
    tipset = TipSet(sizex, sizey, points_list, maxdistweight)
    for i in range(len(points_list)-1):
        G = tipset.get_max_bipartite_graph(framedif)
        matches = tipset.get_max_weight_matches(G)
        tipset.add_matches(matches, framedif)

    with open(filename + "-" + str(maxdistnm) + "nm-"+ str(int(framedif))+"frames.log", "w") as f:
        for point in tipset.match_list:
            f.write(str([get_offset_from_xyz(item, sizex, sizey) for item in point])+"\n")

    print "Tips matched successfully with maximum distance of", maxdistnm, "and maximum gap of", framedif-1, "frames, " \
          "output was saved as", str(filename + "-" + str(maxdistnm) + "nm-"+str(int(framedif))+"frames.log")
