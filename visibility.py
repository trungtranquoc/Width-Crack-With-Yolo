import numpy as np
import matplotlib.pyplot as plt
import functools
import math

largest_range = 1000000000
delta = 0.0000005

def onSegment(p1, p2, p):
    """
    Check if the point p lies in the axis-aligned segment with endpoints p1, p2
    """
    a = (p1[0] - p[0], p1[1] - p[1])
    b = (p2[0] - p[0], p2[1] - p[1])    
    x1 = a[0] * b[1] 
    x2 = a[1] * b[0] 
    return ((x1 >= x2 - delta and x1 <= x2 + delta) and p[0] >= min(p1[0] - delta,p2[0] - delta) and p[0] <= max(p1[0]+delta,p2[0]+delta)) and p[1] >= min(p1[1]-delta,p2[1]-delta) and p[1] <= max(p1[1]+delta,p2[1]+delta)

#function calculate intersection
def line_intersection(p1, p2, q1, q2):
    """
    Find the intersection between the line containing p1, p2 and the line containing
    q1, q2
    """
    xdiff = (float(p1[0] - p2[0]), float(q1[0] - q2[0]))          #directed vector
    ydiff = (float(p1[1] - p2[1]), float(q1[1] - q2[1]))
    
    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)

    if(div == 0):
        return (largest_range, largest_range)

    d = (det(p1, p2), det(q1, q2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return (x, y)

#Distance between 2 points:
def distance(a,b):
    xdiff = a[0] - b[0]
    ydiff = a[1] - b[1]
    return abs(math.sqrt(xdiff * xdiff + ydiff * ydiff))

#Check whether 2 points are in same side of ray from another point:
def check_side(p,a,b):
    x_a = p[0] - a[0]
    x_b = p[0] - b[0]
    x_p = x_a + x_b
    return x_p >= max(x_a, x_b) or x_p <= min(x_a,x_b)

#check angle of 2 line:
def angle_of_ray(p,q1,q2, incase = 0):
    xdiff = (q1[0] - p[0], q1[1] - p[1])
    ydiff = (q2[0] - p[0], q2[1] - p[1])

    if abs(q1[0] - q2[0]) < delta and abs(q1[1] - q2[1]) < delta:
        print("Error in case: ", q1, " and ", q2, " with p = ", p, " in case ", incase)

    angle = (math.asin((xdiff[0] * ydiff[1] - xdiff[1] * ydiff[0])/math.sqrt((xdiff[0] * xdiff[0] + xdiff[1] * xdiff[1]) * (ydiff[0] * ydiff[0] + ydiff[1] * ydiff[1]))))
    
    diff = xdiff[0] + ydiff[0]
    if(diff > min(xdiff[0], ydiff[0]) and diff < max(xdiff[0], ydiff[0])):
        angle = math.pi - angle
    while (angle < 0):
        angle += 2*math.pi
    return angle

#check whether 2 segments intersect
def check_intersect(p1, p2, q1, q2):
    x = line_intersection(p1, p2, q1, q2)
    if x[0] == largest_range:
        return False
    else:
        return onSegment(p1,p2,x) and onSegment(q1,q2,x)

#function setting up visibility polygon from point p
def visibility_polygon(p, pol):
    visibility_set = []
    number_of_edge = len(pol)
    min_distance = largest_range
    p1 = (p[0]+1,p[1])
    p2 = (largest_range, p[1])
    test_set = []
    
    #Tìm giá trị đầu tiên bắt đầu duyệt bằng cách lấy giao điểm của tia song song Oy kẻ từ p và đa giác.
    for track in range (number_of_edge):
        intersect = line_intersection(p,p1,pol[track], pol[(track + 1) % number_of_edge])
        if intersect[0] == largest_range:
            continue
        distance_from_p = distance(p,intersect)
        if onSegment(pol[track], pol[(track + 1) % number_of_edge], intersect) and check_side(p,intersect,p1):
            if min_distance == distance_from_p:
                if angle_of_ray(pol[track], pol[(track + 1) % number_of_edge], pol[(track - 1) % number_of_edge], 1) <= math.pi:
                    init_vertex_1 = pol[track]
                    init_vertex_2 = pol[(track + 1) % number_of_edge]
                    index = track
            if min_distance > distance_from_p:
                init_point = intersect
                min_distance = distance_from_p
                init_vertex_1 = pol[track]
                init_vertex_2 = pol[(track + 1) % number_of_edge]
                index = track
    visibility_set.append(init_point)
    test_set.append(init_point)

    #Check wise of graph, then, in each while step, we will plus index of point to check_wise:
    #Case 1: Start checking from point init_vertex_2:
    if init_vertex_2[1] > p[1] or init_vertex_1[1] < p[1]:
        visibility_set.append(init_vertex_2)
        test_set.append(init_vertex_2)
        init_index = (index+1) % number_of_edge
        check_wise = 1
    #Case 2: Start checking from point init_vertex_1:
    else:
        test_set.append(init_vertex_1)
        visibility_set.append(init_vertex_1)
        init_index = index
        check_wise = -1

    #General Step: 
    track = (init_index + check_wise) % number_of_edge

    
    # Trong thuật toán H.Elgindy, ta duyệt các qua từng cạnh với độ phức tạp O(n) và chia các trường hợp để xử lí
    # Kết quả lúc duyệt được lưu vào trong visibility_set, nó hoat động như stack. Nghĩa là ở mỗi bước, ta push phần tử mới vào hoặc pop các phần tử cũ ra.

    is_visible = True
    max_angle = 0
    # Dùng để lưu lại điểm đang xét
    # test_point_1 = [0,0]
    # test_point_2 = [0,0]

    while track != init_index:
        vis_len = len(visibility_set)
        # Trường hợp duyệt qua các cạnh nhỏ hơn góc 360 độ
        if is_visible:

            if angle_of_ray(p, p1, visibility_set[vis_len-1], 2) > max_angle:
                max_angle = angle_of_ray(p, p1, visibility_set[vis_len-1])
            # if abs(visibility_set[vis_len-1][0] - pol[track][0]) < delta and abs(visibility_set[vis_len-1][0] - pol[track][0]) < delta:
            #     track = (track + check_wise) % number_of_edge
            #     continue
            angle_diff = angle_of_ray(p, visibility_set[vis_len-1], pol[track], 3)

            # Trường hợp 1: Các cạnh xoay theo chiều ngược kim đồng hồ, ta xét 2 trường hợp nhỏ
            if angle_diff <= math.pi:
                # Trường hợp 1.1: Điểm tiếp theo nằm ở khoảng góc lớn hơn 360 độ nên ta set is_visible = 0 và chuyển sang xét các ở trường hợp lớn ở dưới
                if pol[track][1] > p[1] and visibility_set[vis_len-1][1] <= p[1] and max_angle > math.pi:
                    intersect_1 = line_intersection(p, p1, pol[track], visibility_set[vis_len-1])
                    visibility_set.append(intersect_1)
                    prev_point = pol[track]
                    is_visible = False

                # Trường hợp 1.2: Cạnh thực sự được nhìn thấy bởi p nên ta add vào tập visibility_set
                else:
                    visibility_set.append(pol[track])

            # Trường hợp 2: Các cạnh xoay theo chiều ngược kim đồng hồ, ta xét 2 trường hợp nhỏ
            else:
                # Trường hợp 2.1: Cạnh này visible từ p và nằm phía sau các cạnh trong visibility_set
                if (visibility_set[vis_len-1][0] == visibility_set[vis_len-2][0] and visibility_set[vis_len-1][1] == visibility_set[vis_len-2][1]) or angle_of_ray(visibility_set[vis_len-1], visibility_set[vis_len-2], pol[track], 4) <= math.pi:
                    # Ta duyệt đến cạnh tiếp theo visible từ p
                    while angle_of_ray(p, visibility_set[vis_len-1], pol[track], 5) > math.pi:
                        last_point = pol[track]
                        track = (track + check_wise) % number_of_edge

                    # Ta lấy giao điểm của tia từ p đi qua điểm của cùng của visibility_set với đoạn thảng nói trên. Sau đó push điểm này vào trong visibility_set
                    intersect_point = line_intersection(p, visibility_set[vis_len-1], pol[track], last_point)
                    visibility_set.append(intersect_point)
                    track = (track - check_wise) % number_of_edge

                # Trường hợp 2.2: Cạnh này visible từ p và nằm chắn lên các cạnh trước trong visibility_set
                else:
                    first_removed_point = visibility_set[vis_len - 1]
                    last_removed_point = visibility_set[vis_len-1]
                    is_hidden = False
                    visibility_set.pop()
                    track_index = vis_len-2
                    angle_different = angle_of_ray(p ,visibility_set[track_index], pol[track], 6)

                    # Ta vừa duyệt vừa xóa phần tử trong visibility_set
                    while angle_different > math.pi and angle_different < 2*math.pi - delta and is_hidden == False:
                        last_removed_point = visibility_set[track_index]       #last_point is point next to track_index
                        visibility_set.pop()
                        track_index -= 1
                        angle_different = angle_of_ray(p ,visibility_set[track_index], pol[track], 7)
                        is_hidden = check_intersect(first_removed_point, pol[track], visibility_set[track_index], last_removed_point)

                    # Trường hợp 2.2.1: Việc pop phần tử dừng lại khi một cạnh trong visibility_set mà cắt cạnh đang xét.
                    # Lưu ý: Cạnh trong visibility_set không phải cạnh trong polygon ban đầu và chắc chắn nó thẳng hàng với điểm p. 
                    if is_hidden:
                        intersect_2 = line_intersection(first_removed_point, pol[track], visibility_set[track_index], last_removed_point)
                        prev_point = pol[track]
                        track = (track + check_wise) % number_of_edge
                        
                        # Trong trường hợp này, ta tiến hành duyệt tiếp cạnh trong polygon đến khi có cạnh tiếp theo visible từ p.
                        while check_intersect(intersect_2, visibility_set[track_index], prev_point, pol[track]) == False:
                            prev_point = pol[track]
                            track = (track + check_wise) % number_of_edge
                        intersect_3 = line_intersection(intersect_2, visibility_set[track_index], prev_point, pol[track])
                        visibility_set.append(intersect_3)
                        track = (track - check_wise) % number_of_edge 

                    # Trường hợp 2.2.2: Việc pop phần tử dừng lại khi đụng một cạnh trong visibility_set không bị chắn hoàn toàn bởi cạnh đang xét
                    if is_hidden == False:
                        intersect_point = line_intersection(p, pol[track], visibility_set[track_index], last_removed_point)     
                        # Ta push giao điểm từ p đến cạnh vừa xét và push điểm cuối cùng của cạnh đang xét.                    
                        visibility_set.append(intersect_point)
                        visibility_set.append(pol[track])
                        # if abs(pol[track][0] - 438) < delta and abs(pol[track][1] - 237) < delta:
                            # test_point_1 = intersect_point
                            # test_point_2 = pol[(track + check_wise) % number_of_edge]
                            # test_point_3 = pol[track]

        # Trường hợp duyệt qua các cạnh góc lớn hơn 360 độ và bị khuất bởi các cạnh đầu tiên.
        else:
            if check_intersect(p, visibility_set[vis_len-1], pol[track], prev_point):
                p2 = pol[track]
                track = (track - check_wise) % number_of_edge
                is_visible = True
            else:
                prev_point = pol[track]
        if abs(pol[track][0] - pol[(track + check_wise) % number_of_edge][0]) < delta and abs(pol[track][1] - pol[(track + check_wise) % number_of_edge][1]) < delta and track != init_index:
            track = (track + check_wise) % number_of_edge
        track = (track + check_wise) % number_of_edge

    #Delete 2 adjacent point with same coordinate:
    vis_len = len(visibility_set)
    visibility_set_shorten = []
    for track in range(vis_len):
        len_shorten = len(visibility_set_shorten)
        if len_shorten == 0:
            visibility_set_shorten.append(visibility_set[track])
        elif visibility_set[track][0] != visibility_set_shorten[len_shorten-1][0] or visibility_set[track][1] != visibility_set_shorten[len_shorten-1][1]:
            if track < vis_len-1:
                visibility_set_shorten.append(visibility_set[track])
            elif visibility_set[track][0] != visibility_set_shorten[0][0] or visibility_set[track][1] != visibility_set_shorten[0][1]:
                visibility_set_shorten.append(visibility_set[track])

    return visibility_set_shorten
