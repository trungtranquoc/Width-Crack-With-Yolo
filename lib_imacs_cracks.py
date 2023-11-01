import numpy as np
import cv2
from skimage.morphology import skeletonize
from shapely.geometry import Point , LineString,  Polygon, MultiPoint, mapping, shape
from scipy.spatial.distance import directed_hausdorff

def hh_skeletonize(img, points):
    Ny = img.shape[0]
    Nx = img.shape[1]
    img_skel = np.zeros((Ny, Nx, 3), np.uint8)

    cv2.fillPoly(img_skel, pts=[points], color=(255, 255, 255))
    img_skel = img_skel.astype(bool)
    skeleton = skeletonize(img_skel).astype(np.uint8)#,method='lee'
    skeleton_points = np.argwhere(skeleton > 0)
    skeleton_points[:, [0, 1]] = skeleton_points[:, [1, 0]]

    return skeleton_points

def hh_get_segment_angle(p, points, theta):
    poly = Polygon(points)
    poly_bounds = poly.bounds
    rho = np.sqrt((poly_bounds[0]-poly_bounds[2])**2+(poly_bounds[1]-poly_bounds[3])**2)
    q1x = p[0] + rho*np.cos(theta)
    q1y = p[1] + rho*np.sin(theta)
    q2x = p[0] + rho*np.cos(theta+np.pi)
    q2y = p[1] + rho*np.sin(theta+np.pi)
    q1, q2 = Point(q1x,q1y), Point(q2x, q2y)
    p = Point(p[0],p[1])

    line = LineString([q2,q1])
    p_intersections = line.intersection(poly)

    if p_intersections.geom_type=='LineString':
        line = p_intersections
    elif p_intersections.geom_type=='GeometryCollection':
        line = LineString([(0,0),(0,0)])
        dis = 0
        for kk in p_intersections.geoms:
            if kk.geom_type != 'Point':
                iline = kk
                dis_new = iline.distance(p)
                if iline.distance(p) < dis:
                    dis=dis_new
                    line = iline
    else:
        ppp = mapping(p_intersections)['coordinates']
        
        point1 = Point(ppp[0][0])
        point2 = Point(ppp[0][1])
        line = LineString([point1,point2])
        dis = line.distance(p)
        for i in range(1,len(ppp)):
            point1 = Point(ppp[i][0])
            point2 = Point(ppp[i][1])
            iline = LineString([point1,point2])
            dis_new = iline.distance(p)
            if iline.distance(p) < dis:
                dis=dis_new
                line = iline

    def LINESTRING_point(line):
        line = mapping(line)
        line = line['coordinates']
        xx = ((line[0][0]),(line[0][1]))
        yy = ((line[1][0]),(line[1][1]))
        return [yy,xx]
    return LINESTRING_point(line),line.length
        

class hh_width_point_approx(object):
    def __init__(self, p, points, N):
        self.p = p
        self.points = points
        self.N = N

        angle_list = [2 * np.pi * i / 2 ** self.N for i in range(2 ** (self.N - 1) + 1)]
        line_min, line_length_min = hh_get_segment_angle(self.p, self.points, angle_list[0])

        angle_list.pop(0)
        self.all = []
        segment_min = None
        for theta in angle_list:
            line, line_length = hh_get_segment_angle(self.p, self.points, theta)
            if line==None:
                print(line_length)
            self.all = self.all + [np.array(line)]
            if line_length < line_length_min:
                segment_min, line_length_min = line, line_length
        self.min = np.array(segment_min)