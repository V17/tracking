# -*- coding: utf-8 -*-
__author__ = 'Vojtech Vozab'
import tracking
import argparse

if __name__ == '__main__':
    defaults = {'framegap': 1, 'maxdistance': 2000, 'log': True, 'inputfile': None, 'x': None, 'y': None, 'zr': 0.125}
    parser = argparse.ArgumentParser(description='A script for tracking a sequence of tip coordinates')
    parser.add_argument('inputfile', help='Text file containing sequences of tip offsets')
    parser.add_argument('x', type=int, help='y-size of source image')
    parser.add_argument('y', type=int, help='x-size of source image')
    parser.add_argument('--framegap', '-f', type=int, help='Maximum temporal distance between two tip coordinates identified as the same tip.'
                                                           'Default = 1 means that the tip cannot vanish at all, 2 means that it can vanish'
                                                           'for a single frame etc.')
    parser.add_argument('--maxdistance', '-d', type=int, help='Maximum distance between 2 frames in nanometers.')
    parser.add_argument('--log', '-l', type=bool, help='If set to False, output is only printed to console.')
    parser.add_argument('--zratio', '-zr', type=float, help='The ration of x and y voxel dimensions to z dimension in anisotropic images.'
                                                            'Default is set to 0.125')
    namespace = parser.parse_args()
    command_line_args = {k: v for k, v in vars(namespace).items() if v}
    d = defaults.copy()
    d.update(command_line_args)

    pointset = tracking.read_points(d['inputfile'], d['x'], d['y'])
    tracking.apply_and_write(d['inputfile'], d['x'], d['y'], pointset, d['framegap'], d['maxdistance'], d['log'])