from tracking import *

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

files_list = files_list1

"""results = generate_results(files_list[0][0], files_list[0][1], files_list[0][2])
write_results(results, files_list[0][0])

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
