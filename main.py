# -*- coding: utf-8 -*-
__author__ = 'Vojtech Vozab'

import networkx as nx
from itertools import izip_longest
from math import sqrt


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

    def get_max_bipartite_graph(self):
        # Initialize lists of existing points and new points, add dummy points into new points list
        old_points = self.get_last_frame()
        new_points = []
        for coord in self.frame_list[self.timeframe+1]:
            if coord != (-1, -1, -1):
                new_points.append(coord)
        for i in range(len(old_points)):
            new_points.append((-1, -1, -1))
        # print "old points", len(old_points), old_points
        # print "new points", len(new_points), new_points

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

    def add_matches(self, matches_dict):
        # add matched pairs and vanishing points into match_list, find and add new points, increase the timeframe
        new_points = []
        for coord in self.frame_list[self.timeframe+1]:
            if coord != (-1, -1, -1):
                new_points.append(coord)
        for match in self.match_list:
            if match[self.timeframe] != (-1, -1, -1):
                match.append(matches_dict[match[self.timeframe]])
                try:
                    if matches_dict[match[self.timeframe]] != (-1, -1, -1):
                        new_points.remove(matches_dict[match[self.timeframe]])
                except ValueError:
                    print "something broke, coord not in new_points list:", matches_dict[match[self.timeframe]]
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

    def get_last_frame(self, framedif=1):  # TODO: implement looking back
        coordlist = []
        for coords in self.match_list:
            if coords[self.timeframe] != (-1, -1, -1):
                coordlist.append(coords[self.timeframe])
        return coordlist

    def get_pairs_for_frame(self, timeframe):
        pairs = []
        for coord_sequence in self.match_list:
            pairs.append(([coord_sequence[timeframe-1]], coord_sequence[timeframe]))
        return



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


def check_found_and_control(found, control):
    for frame in range(1, len(found.match_list[0])-1):
        found_pairs = found.get_pairs_for_frame(frame)
        control_pairs = control.get_pairs_for_frame(frame)
        for pair in found_pairs:
            if pair[0] != (-1, -1, -1) and pair[1] != (-1, -1, -1):
                pass




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


def check_points(coord_pairs, points, frame, sizex, sizey):
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


def convert_match_to_offset(match, sizex, sizey):
    return get_offset_from_xyz(match[0], sizex, sizey), get_offset_from_xyz(match[1], sizex, sizey)


def check_correct_matches(detected_tipset, control_tipset):
    correct_matches = 0
    wrong_matches = 0
    incorrectly_found = 0
    incorrectly_vanished = 0
    for i in range(1, len(detected_tipset.tip_list[0].coords_list)):
        detected = detected_tipset.get_matches_for_frame(i)
        control = control_tipset.get_matches_for_frame(i)
        sizex = detected_tipset.sizex
        sizey = detected_tipset.sizey
        print "comparing matches in frames", i-1, "and", i
        for match1 in detected:
            if match1[0] != (-1, -1, -1) and match1[1] != (-1, -1, -1):
                for match2 in control:
                    if match1[0] == match2[0] and match1[0] != (-1, -1, -1):
                        if match1[1] == match2[1]:
                            print "correctly matched", convert_match_to_offset(match1, sizex, sizey)
                            correct_matches += 1
                        else:
                            print "wrongly matched", convert_match_to_offset(match1, sizex, sizey), "should be", convert_match_to_offset(match2, sizex, sizey)
                            wrong_matches += 1
            elif match1[0] == (-1, -1, -1) and match1[1] != (-1, -1, -1):
                for match2 in control:
                    if match2[0] == (-1, -1, -1) and match1[1] == match2[1]:
                        print "new tip found at", get_offset_from_xyz(match1[1], sizex, sizey)
                    if match2[1] == match1[1] and match2[0] != (-1, -1, -1):
                        print "existing tip incorrectly recognized as new at", get_offset_from_xyz(match1[1], sizex, sizey)
                        incorrectly_found += 1

            elif match1[0] != (-1, -1, -1) and match1[1] == (-1, -1, -1):
                for match2 in control:
                    if match2[0] == match1[0] and match2[1] != (-1, -1, -1):
                        print "tip incorrectly recognized as vanishing at", get_offset_from_xyz(match1[0], sizex, sizey)
                        incorrectly_vanished += 1
                    if match2[0] == match1[0] and match2[1] == (-1, -1, -1):
                        print "tip correctly recognized as vanishing at", get_offset_from_xyz(match1[0], sizex, sizey)

            else:
                continue  # (-1, -1, -1) matched to (-1, -1, -1)
        print "-----"
    print "Total score:", correct_matches, "matched correctly,", wrong_matches, "matched incorrectly,", incorrectly_found, "incorrectly labeled as new,", incorrectly_vanished, "incorrectly labeled as vanished"


def apply_and_check(filename, sizex, sizey, frames=0):
    points_list = read_points(filename, sizex, sizey)
    tipset = TipSet(sizex, sizey, points_list)
    for i in range(len(points_list)-1):
        print "timeframe is", tipset.timeframe
        G = tipset.get_max_bipartite_graph()
        matches = tipset.get_max_weight_matches(G)
        tipset.add_matches(matches)
    control = build_control_tipset(points_list, sizex, sizey)
    print control.match_list


apply_and_check("./tip_tracking/ID319/tip_trajectories_pxOffsets.txt", 159, 146, 0)

