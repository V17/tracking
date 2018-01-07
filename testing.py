# -*- coding: utf-8 -*-
__author__ = 'Vojtech Vozab'
from tracking import *
import time

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


def check_found_and_control(found, control, frame, framedif):
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
        elif pair[0] == (-1, -1, -1) and pair[1] != (-1, -1, -1):
            matches_index_found = found.get_complete_frame(frame).index(pair[1])
            matches_index_ctrl = control.get_complete_frame(frame).index(pair[1])
            last_coord = found.get_last_location(found.match_list[matches_index_found], -1, frame-1)
            last_real_coord = control.get_last_location(control.match_list[matches_index_ctrl], -1, frame - 1)
            if last_coord == (-1, -1, -1):
                if last_real_coord == (-1, -1, -1):
                    correctly_found += 1
                else:
                    incorrectly_found += 1
            else:
                if last_coord == last_real_coord:
                    correctly_matched += 1
                else:
                    incorrectly_matched += 1

        elif pair[0] != (-1, -1, -1) and pair[1] == (-1, -1, -1):
            if pair in control_pairs:
                correctly_vanished += 1
            else:
                incorrectly_vanished += 1
        else:  # (-1, -1, -1) matched with (-1, -1, -1) = point doesn't exist in this timeframe
            pass

    return [correctly_matched, incorrectly_matched, correctly_vanished, incorrectly_vanished, correctly_found, incorrectly_found]


def check_found_and_control_f1(found, control, frame, framedif):
    found_pairs2 = found.get_pairs_for_frame(frame)
    control_pairs2 = control.get_pairs_for_frame(frame)
    found_pairs = []
    control_pairs = []
    tp = 0
    tn = 0
    fp = 0
    fn = 0

    for pair in found_pairs2:
        if pair[0] == (-1, -1, -1) and pair[1] != (-1, -1, -1):
            matches_index_found = found.get_complete_frame(frame).index(pair[1])
            last_coord = found.get_last_location(found.match_list[matches_index_found], framedif, frame - 1)
            found_pairs.append((last_coord, pair[1]))
        else:
            found_pairs.append(pair)
    for pair in control_pairs2:
        if pair[0] == (-1, -1, -1) and pair[1] != (-1, -1, -1):
            matches_index_control = control.get_complete_frame(frame).index(pair[1])
            last_coord = control.get_last_location(control.match_list[matches_index_control], -1, frame - 1)
            control_pairs.append((last_coord, pair[1]))
        else:
            control_pairs.append(pair)

    for pair in found_pairs:
        if pair[0] != (-1, -1, -1) and pair[1] != (-1, -1, -1):
            # edge between two existing points
            if pair in control_pairs:
                tp += 1
            else:
                if (pair[0], (-1, -1, -1)) in control_pairs:
                    fp += 1
                else:
                    fp += 1
                    fn += 1

        elif pair[0] == (-1, -1, -1) and pair[1] != (-1, -1, -1):
            # marked as new, edge between existing point and (-1, -1, -1)
            if pair in control_pairs:
                tn += 1
            else:
                pass
                # incorrectly marked as new, no change in results because the false positive edge is handled
                # in one of the other cases

        elif pair[0] != (-1, -1, -1) and pair[1] == (-1, -1, -1):
            # marked as vanishing
            if pair in control_pairs:
                tn += 1
            else:
                fn += 1

        else:
            # (-1, -1, -1) matched with (-1, -1, -1) = point doesn't exist in this timeframe
            pass

    return [tp, tn, fp, fn]


def apply_and_check(filename, sizex, sizey, points_list, xypxsize=125, zratio=0.125, framedif=1, maxdistnm=6250, logging=True, check_n_frames=0):
    #points_list = read_points(filename, sizex, sizey)
    if len(points_list) == 0 and logging:
        logfile = open(filename + "." + str(int(time.time())) + ".log", "w")
        logfile.write("empty list")
        logfile.close()
        return 0
    maxdistpx = 1.0/(float(maxdistnm)/float(xypxsize))

    tipset = TipSet(sizex, sizey, points_list, maxdistpx)
    if check_n_frames == 0 or check_n_frames > len(points_list)-1:
        check_n_frames = len(points_list)-1
    for i in range(check_n_frames):
        G = tipset.get_max_bipartite_graph(framedif)
        matches = tipset.get_max_weight_matches(G)
        tipset.add_matches(matches, framedif)

    # checking section
    results = [0, 0, 0, 0, 0, 0]
    results_f1 = [0., 0., 0., 0.]
    control = build_control_tipset(points_list, sizex, sizey)
    for frame in range(1, check_n_frames+1):
        tempresults = check_found_and_control(tipset, control, frame, framedif)
        tempresults2 = check_found_and_control_f1(tipset, control, frame, framedif)
        results = map(add, results, tempresults)
        results_f1 = map(add, results_f1, tempresults2)

    try:
        precision = results_f1[0]/(results_f1[0]+results_f1[2])
        recall = results_f1[0]/(results_f1[0]+results_f1[3])
        print "precision, recall", precision, recall
        f1score = (2*precision*recall)/(precision+recall)
    except ZeroDivisionError:
        precision = 0
        recall = 0
        f1score = 0

    if logging:
        logfile = open(filename + "." + str(int(time.time())) + ".log", "w")
        logfile.write(
            "Checking file " + filename + " maximum time gap between two tips is " + str(framedif-1) + " frames, " +
            "maximum distance between two points is " + str(maxdistnm) + "nm. \n")
        logfile.write("correctly matched: " + str(results[0]) + " incorrectly matched: " + str(results[1]) + "\n")
        logfile.write(
            "correctly recognized as new: " + str(results[4]) + " incorrectly recognized as new: " + str(
                results[5]) + "\n")
        logfile.write(
            "correctly recognized as vanishing: " + str(results[2]) + " incorrectly recognized as vanishing: " + str(
                results[3]) + "\n")
        logfile.write("precision: "+str(precision)+", recall: "+str(recall)+", F1 score: "+str(f1score)+"\n")
        logfile.write("-------------------" + "\n")

        for point in tipset.match_list:
            logfile.write(str([get_offset_from_xyz(item, sizex, sizey) for item in point])+"\n")
        logfile.close()
    return f1score


def generate_results(filename, sizex, sizey):
    results = [[0]+range(1000, 4001, 125)]
    points_list = read_points(filename, sizex, sizey)
    for framegap in range(1, 10):
        temp = [framegap-1]
        for dist in range(1000, 4001, 125):
            f1score = apply_and_check(filename, sizex, sizey, points_list, framedif=framegap, maxdistnm=dist, logging=False)
            print "computing using framedif and maxdis", framegap, dist, "f1 score", f1score
            temp.append(f1score)
        results.append(temp)
    return results


def print_results(results):
    print "{0:11s}".format("")+" | ",
    for i in range(1, len(results[0])):
        print "{0:5d}".format(results[0][i])+" nm | ",
    print ""
    for line in range(1, len(results)):
        for i in range(len(results[line])):
            if i == 0:
                print "{0:4d}".format(results[line][i])+" frames | ",
            else:
                print "{0:8.4f}".format(results[line][i])+" | ",
        print ""


def write_results(results, filename):
    result_file = open(filename+"-result.txt", 'w')
    result_file.write("{0:11s}".format("")+" | ")
    for i in range(1, len(results[0])):
        result_file.write("{0:5d}".format(results[0][i])+" nm | ")
    result_file.write("\n")
    for line in range(1, len(results)):
        for i in range(len(results[line])):
            if i == 0:
                result_file.write("{0:4d}".format(results[line][i])+" frames | ")
            else:
                result_file.write("{0:8.4f}".format(results[line][i])+" | ")
        result_file.write("\n")
    result_file.close()


def add_2_results(results1, results2):
    for i in range(1, len(results1)):
        for j in range(1, len(results1[0])):
            results1[i][j] = (results2[i][j] + results1[i][j])/2
    return results1


def timer(func):
    def inner_func(*args, **kwargs):
        current_time = time.time()
        result = func(*args, **kwargs)
        print("Results generated in {:.1f}s".format(time.time() - current_time))
        return result
    return inner_func


files_list1 = [["./real_images/real-control-01.txt", 220, 220],
               ["./real_images/real-control-02.txt", 220, 220],
               ["./real_images/real-phosphodefective-01.txt", 300, 400],
               ["./real_images/real-overexpressing-01.txt", 350, 300],
               ["./real_images/real-overexpressing-02.txt", 350, 300],
               ["./real_images/real-overexpressing-03.txt", 350, 300],
               ["./real_images/real-phosphodefective-02.txt", 500, 600],
               ["./real_images/real-phosphodefective-03.txt", 400, 270]]

files_list2 = [["./tip_tracking/ID319/tip_trajectories_pxOffsets.txt", 159, 146],
               ["./tip_tracking/ID320/tip_trajectories_pxOffsets.txt", 157, 148],
               ["./tip_tracking/ID321/tip_trajectories_pxOffsets.txt", 160, 148],
               #["./tip_tracking/ID426/tip_trajectories_pxOffsets.txt", 295, 317],
               #["./tip_tracking/ID427/tip_trajectories_pxOffsets.txt", 271, 275],
               ["./tip_tracking/ID428/tip_trajectories_pxOffsets.txt", 314, 261]]

files_list = files_list2
generate_results_timed = timer(generate_results)
results = generate_results_timed(files_list[0][0], files_list[0][1], files_list[0][2])
write_results(results, files_list[0][0])
"""
for i in range(1, len(files_list)):
    results2 = generate_results(files_list[i][0], files_list[i][1], files_list[i][2])
    write_results(results, files_list[i][0])
    results = add_2_results(results, results2)

print_results(results)
write_results(results, "./real_images/total")


# u ID426 to je duplikaci offsetu 641721
# u ID427 je duplikace 1993618, 318888
"""
# print apply_and_check(files_list2[3][0], files_list2[3][1], files_list2[3][2], framedif=1, maxdistnm=2000, logging=True, check_n_frames=0)
#print apply_and_check(filename, size_x, size_y, framedif=2, maxdistnm=1000, logging=True)
