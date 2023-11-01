import numpy as np
import matplotlib.pyplot as plt
import functools
import math

largest_range = 1000000000;
delta = 0.0000005;

def onSegment(p1, p2, p):
    """
    Check if the point p lies in the axis-aligned segment with endpoints p1, p2
    """
    a = (p1[0] - p[0], p1[1] - p[1]);
    b = (p2[0] - p[0], p2[1] - p[1]);    
    x1 = a[0] * b[1] ;
    x2 = a[1] * b[0] ;
    return ((x1 >= x2 - delta and x1 <= x2 + delta) and p[0] >= min(p1[0] - delta,p2[0] - delta) and p[0] <= max(p1[0]+delta,p2[0]+delta)) and p[1] >= min(p1[1]-delta,p2[1]-delta) and p[1] <= max(p1[1]+delta,p2[1]+delta);

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
        return (largest_range, largest_range);

    d = (det(p1, p2), det(q1, q2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return (x, y);

#Distance between 2 points:
def distance(a,b):
    xdiff = a[0] - b[0];
    ydiff = a[1] - b[1];
    return abs(math.sqrt(xdiff * xdiff + ydiff * ydiff));

#Check whether 2 points are in same side of ray from another point:
def check_side(p,a,b):
    x_a = p[0] - a[0];
    x_b = p[0] - b[0];
    x_p = x_a + x_b;
    return x_p >= max(x_a, x_b) or x_p <= min(x_a,x_b)

#check angle of 2 line:
def angle_of_ray(p,q1,q2):
    xdiff = (q1[0] - p[0], q1[1] - p[1]);
    ydiff = (q2[0] - p[0], q2[1] - p[1]);

    angle = (math.asin((xdiff[0] * ydiff[1] - xdiff[1] * ydiff[0])/math.sqrt((xdiff[0] * xdiff[0] + xdiff[1] * xdiff[1]) * (ydiff[0] * ydiff[0] + ydiff[1] * ydiff[1]))));
    
    diff = xdiff[0] + ydiff[0];
    if(diff > min(xdiff[0], ydiff[0]) and diff < max(xdiff[0], ydiff[0])):
        angle = math.pi - angle;
    while (angle < 0):
        angle += 2*math.pi;
    return angle;

#check whether 2 segments intersect
def check_intersect(p1, p2, q1, q2):
    x = line_intersection(p1, p2, q1, q2);
    if x[0] == largest_range:
        return False;
    else:
        return onSegment(p1,p2,x) and onSegment(q1,q2,x);

#function setting up visibility polygon from point p
def visibility_polygon(p, pol):
    visibility_set = [];
    number_of_edge = len(pol);
    min_distance = largest_range;
    p1 = (p[0]+1,p[1]);
    p2 = (largest_range, p[1]);
    test_set = [];
    
    #Search first visibility point to start
    for track in range (number_of_edge):
        intersect = line_intersection(p,p1,pol[track], pol[(track + 1) % number_of_edge]);
        if intersect[0] == largest_range:
            continue;
        distance_from_p = distance(p,intersect);
        if onSegment(pol[track], pol[(track + 1) % number_of_edge], intersect) and check_side(p,intersect,p1):
            if min_distance == distance_from_p:
                if angle_of_ray(pol[track], pol[(track + 1) % number_of_edge], pol[(track - 1) % number_of_edge]) <= math.pi:
                    init_vertex_1 = pol[track];
                    init_vertex_2 = pol[(track + 1) % number_of_edge];
                    index = track;
            if min_distance > distance_from_p:
                init_point = intersect;
                min_distance = distance_from_p;
                init_vertex_1 = pol[track];
                init_vertex_2 = pol[(track + 1) % number_of_edge];
                index = track;
    visibility_set.append(init_point);
    test_set.append(init_point);

    #Check wise of graph, then, in each while step, we will plus index of point to check_wise:
    #Case 1: Start checking from point init_vertex_2:
    if init_vertex_2[1] > p[1] or init_vertex_1[1] < p[1]:
        visibility_set.append(init_vertex_2);
        test_set.append(init_vertex_2);
        init_index = (index+1) % number_of_edge;
        set_wise = 1;
        check_wise = 1;
    #Case 2: Start checking from point init_vertex_1:
    else:
        test_set.append(init_vertex_1);
        visibility_set.append(init_vertex_1);
        init_index = index;
        set_wise = -1;
        check_wise = -1;

    #General Step: 
    track = (init_index + check_wise) % number_of_edge;
    is_visible = True;
    max_angle = 0;
    while track != init_index:
        vis_len = len(visibility_set);
        if is_visible:
            if angle_of_ray(p, p1, visibility_set[vis_len-1]) > max_angle:
                max_angle = angle_of_ray(p, p1, visibility_set[vis_len-1]);
            angle_diff = angle_of_ray(p, visibility_set[vis_len-1], pol[track]);
            if angle_diff <= math.pi:
                if pol[track][1] > p[1] and visibility_set[vis_len-1][1] <= p[1] and max_angle > math.pi:
                    intersect_1 = line_intersection(p, p1, pol[track], visibility_set[vis_len-1]);
                    visibility_set.append(intersect_1);
                    prev_point = pol[track];
                    is_visible = False;
                    # check_case.append(pol[track]);
                else:
                    visibility_set.append(pol[track]);
                # visibility_set.append(pol[track]);

            else:
                #Case 3:
                if (visibility_set[vis_len-1][0] == visibility_set[vis_len-2][0] and visibility_set[vis_len-1][1] == visibility_set[vis_len-2][1]) or angle_of_ray(visibility_set[vis_len-1], visibility_set[vis_len-2], pol[track]) <= math.pi:
                    while angle_of_ray(p, visibility_set[vis_len-1], pol[track]) > math.pi:
                        last_point = pol[track];
                        track = (track + check_wise) % number_of_edge;
                    intersect_point = line_intersection(p, visibility_set[vis_len-1], pol[track], last_point);
                    visibility_set.append(intersect_point);
                    track = (track - check_wise) % number_of_edge;
                else:
                    first_removed_point = visibility_set[vis_len - 1];
                    last_removed_point = visibility_set[vis_len-1];
                    is_hidden = False;
                    visibility_set.pop();
                    track_index = vis_len-2;
                    angle_different = angle_of_ray(p ,visibility_set[track_index], pol[track]);
                    while angle_different > math.pi and angle_different < 2*math.pi - delta and is_hidden == False:
                        last_removed_point = visibility_set[track_index];       #last_point is point next to track_index
                        visibility_set.pop();
                        track_index -= 1;
                        angle_different = angle_of_ray(p ,visibility_set[track_index], pol[track])
                        is_hidden = check_intersect(first_removed_point, pol[track], visibility_set[track_index], last_removed_point);

                    if is_hidden:
                        intersect_2 = line_intersection(first_removed_point, pol[track], visibility_set[track_index], last_removed_point);
                        prev_point = pol[track];
                        track = (track + check_wise) % number_of_edge;
                        while check_intersect(intersect_2, visibility_set[track_index], prev_point, pol[track]) == False:
                            prev_point = pol[track];
                            track = (track + check_wise) % number_of_edge;
                        intersect_3 = line_intersection(intersect_2, visibility_set[track_index], prev_point, pol[track]);
                        visibility_set.append(intersect_3);
                        track = (track - check_wise) % number_of_edge; 

                    if is_hidden == False:
                        intersect_point = line_intersection(p, pol[track], visibility_set[track_index], last_removed_point);                         
                        visibility_set.append(intersect_point);
                        visibility_set.append(pol[track]);
                    
        else:
            if check_intersect(p, visibility_set[vis_len-1], pol[track], prev_point):
                p2 = pol[track];
                track = (track - check_wise) % number_of_edge;
                is_visible = True;
            else:
                prev_point = pol[track];
        track = (track + check_wise) % number_of_edge;

    track_complete = track;
    #Delete 2 adjacent point with same coordinate:
    vis_len = len(visibility_set);
    visibility_set_shorten = [];
    for track in range(vis_len):
        len_shorten = len(visibility_set_shorten);
        if len_shorten == 0:
            visibility_set_shorten.append(visibility_set[track]);
        elif visibility_set[track][0] != visibility_set_shorten[len_shorten-1][0] or visibility_set[track][1] != visibility_set_shorten[len_shorten-1][1]:
            if track < vis_len-1:
                visibility_set_shorten.append(visibility_set[track]);
            elif visibility_set[track][0] != visibility_set_shorten[0][0] or visibility_set[track][1] != visibility_set_shorten[0][1]:
                visibility_set_shorten.append(visibility_set[track]);

    return min_distance, init_point, check_wise, visibility_set_shorten;
if __name__ == "__main__":
    # pol = [(0,0), (7,0), (7,6), (6,5), (6,5.5), (5,5), (7,7), (9,6), (7,12), (8,7), (7,9), (-2,8), (0,8), (0,7), (5,7), (-2,5), (1,4), (-2,4), (-3,5), (-3,-3)];
    # p = (3,4);
    pol = [(-1,1), (1,1), (2,-1), (4,-2), (6,1), (6,4), (5,3), (4,1), (3,0), (2,2), (0,4), (-3,4), (-3,3), (-4,2), (-4,-10), (-2,-1)];
    p = (-3, 1);

    min_distance, init_point, check_wise, visibility_set = visibility_polygon(p ,pol);

    print(visibility_set);
    visibility_set.append(visibility_set[0]);
    pol.append(pol[0]);

    x, y = zip(*pol)
    a, b = zip(*visibility_set)
    fig = plt.figure()
    plt.plot(p[0], p[1], marker='o')
    plt.plot(x, y, color='black')
    plt.plot(a, b, color='red')
    plt.show()