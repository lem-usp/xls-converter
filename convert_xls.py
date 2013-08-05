# -*- coding: utf-8 -*-

import xlrd
import string
import ast
from collections import defaultdict
from optparse import OptionParser
import numpy as np

parser = OptionParser()
parser.add_option("-i", "--input", dest="input",
                  help="Input XLS", metavar="FILE")
parser.add_option("-o", "--output", dest="output",
                  help="output CSV")
parser.add_option("-n", "--inds", dest="inds_num",
                  help="Number of individuals")
parser.add_option("-v", "--vistas", dest="vistas",
                  help="Vistas coordinates like (1,'A'), (10,'B')")
parser.add_option("-d", "--dists", dest="distances",
                  help="Distances coordinates like (1,'A'), (10,'B')")

(options, args) = parser.parse_args()

def cell_name_idx(s):
    b = 1
    n = 0
    for i in [string.lowercase.index(string.lower(c)) for c in reversed(s)]:
        n += b*(i+1)
        b *= (string.lowercase.index('z') + 1)

    return n-1

def get_points_names(data, vistas):
    d = defaultdict(list)
    for i, coords in enumerate(vistas):
        line = int(coords[0] + 1)
        col = coords[1]
        p = data.cell(line, col)
        while p.value != '':
            d[data.cell(*vistas[i]).value].append((p.value, (line, col)))
            line += 1
            p = data.cell(line, col)

    return d

def read_all_points(data, vistas, num_mes=1):
    INFO_COL = 1

    points = []
    line = 0

    i = 0

    points_info = get_points_names(data, vistas)
    points_lines_max = max(map(len, points_info.values()))

    parsed = 0

    while parsed < num_mes:
        points.append({})

        # Get Info
        points[-1]['info'] = []
        l = line + 1
        p = data.cell(l, INFO_COL).value

        while p != '':
            points[-1]['info'].append(p)
            l += 1
            p = data.cell(l, INFO_COL).value

        for v in vistas:
            vname = data.cell(*v).value
            points[-1][vname] = {}

            for i, p in enumerate(points_info[vname]):
                points[-1][vname][p[0]] = {}

                # Get X
                points[-1][vname][p[0]]['x'] = []
                points[-1][vname][p[0]]['x'].append(data.cell(p[1][0] + line, p[1][1]+1).value)
                points[-1][vname][p[0]]['x'].append(data.cell(p[1][0] + line, p[1][1]+4).value)

                # Get Y
                points[-1][vname][p[0]]['y'] = []
                points[-1][vname][p[0]]['y'].append(data.cell(p[1][0] + line, p[1][1]+2).value)
                points[-1][vname][p[0]]['y'].append(data.cell(p[1][0] + line, p[1][1]+5).value)

                # Get Z
                points[-1][vname][p[0]]['z'] = []
                points[-1][vname][p[0]]['z'].append(data.cell(p[1][0] + line, p[1][1]+3).value)
                points[-1][vname][p[0]]['z'].append(data.cell(p[1][0] + line, p[1][1]+6).value)

        line += points_lines_max + 3
        parsed += 1

    return points

def get_distance_names(data, columns):
    dists_names = defaultdict(list)
    for c in columns:
        line = c[0] + 1
        col = cell_name_idx(c[1])
        vista = data.cell(line - 2, col).value

        dist = data.cell(line, col).value
        while dist != '':
            dists_names[vista].append(dist)
            line += 1
            dist = data.cell(line, col).value

    for k in dists_names.keys():
        dists_names[k] = map(lambda d: d.split('-'), dists_names[k])
        dists_names[k] = [[p[0].rstrip(), p[1].rstrip()] for p in dists_names[k]]

    return dists_names

def calc_distances(points, dist_names, i):
    ps = points[i]
    distances = {}

    for v in dist_names.keys():
        distances[v] = {}
        for dp in dist_names[v]:
            dist_key = '{}-{}'.format(dp[0], dp[1])

            if dp[1][-1] == 'e' or dp[1][-1] == 'd':
                pk = dp[0] + ' ' + dp[1][-1]
                if ps[v].has_key(pk):
                    if ps[v].has_key(dp[1]):
                        point1_key = pk
                        point2_key = dp[1]
                    else:
                        point1_key = pk
                        point2_key = dp[1][:-2]
                else:
                    point1_key = dp[0]
                    point2_key = dp[1]
            else:
                point1_key = dp[0]
                point2_key = dp[1]

            dist1_x = ps[v][point1_key]['x'][0]- ps[v][point2_key]['x'][0]
            dist1_y = ps[v][point1_key]['y'][0]- ps[v][point2_key]['y'][0]
            dist1_z = ps[v][point1_key]['z'][0]- ps[v][point2_key]['z'][0]

            dist2_x = ps[v][point1_key]['x'][1]- ps[v][point2_key]['x'][1]
            dist2_y = ps[v][point1_key]['y'][1]- ps[v][point2_key]['y'][1]
            dist2_z = ps[v][point1_key]['z'][1]- ps[v][point2_key]['z'][1]

            distances[v][dist_key] = []
            distances[v][dist_key].append(np.sqrt(dist1_x**2 + dist1_y**2 + dist1_z**2))
            distances[v][dist_key].append(np.sqrt(dist2_x**2 + dist2_y**2 + dist2_z**2))

    return distances

def avg_distances(distances):
    avg_distances = {}
    for vista in distances.keys():
        avg_distances[vista] = {}
        for d in distances[vista].keys():
            if d[-1] == 'e':
                if distances[vista].has_key(d[:-1] + 'd'):
                    avg_distances[vista][d[:-1]] = np.average([np.average(distances[vista][d]), np.average(distances[vista][d[:-1] + 'd'])])
                else:
                    avg_distances[vista][d] = np.average(distances[vista][d])
            elif d[-1] == 'd':
                if distances[vista].has_key(d[:-1] + 'e'):
                    avg_distances[vista][d[:-1]] = np.average([np.average(distances[vista][d]), np.average(distances[vista][d[:-1] + 'e'])])
                else:
                    avg_distances[vista][d] = np.average(distances[vista][d])
            else:
                avg_distances[vista][d] = np.average(distances[vista][d])

    return avg_distances

workbook = xlrd.open_workbook(options.input)
data = workbook.sheet_by_name('dados')

vistas_c = ast.literal_eval(options.vistas)
vistas_c = [(v[0], cell_name_idx(v[1])) for v in vistas_c]
dists_c = ast.literal_eval(options.distances)

dist_names = get_distance_names(data, dists_c)
points = read_all_points(data, vistas_c, num_mes=int(options.inds_num))
dists = [calc_distances(points, dist_names, i) for i in xrange(int(options.inds_num))]
avg_dist = [avg_distances(d) for d in dists]

print avg_dist
