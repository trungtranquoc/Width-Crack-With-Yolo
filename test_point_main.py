"""
This code is copyrighted by Institute of Mathematical and Computational Sciences - IMACS, 
    Ho Chi Minh City University of Technology (HCMUT).  
Contact: Prof. Phan Thanh An thanhan@hcmut.edu.vn
"""
from ultralytics import YOLO
import cv2
import numpy as np
import matplotlib.pyplot as plt
from lib_imacs_cracks import *
from width_crack import *

#Nhập thông tin
path_import = 'fig/'#Thư mục chứa ảnh
path_export = 'fig-width/'#Thư mục trả kết quả

# img_name = 'vd2.jpg'
img_name = '7.png'

img_path = path_import+img_name
original_image = cv2.imread(img_path)
k_width = 0.011754334835464791

#Nhận diện vết nứt S
model = YOLO("crack_width_YOLOv8.pt")
cracks = model.predict(source=img_path, conf=0.4)[0]

crack = cracks.masks.xy[0]
crack = np.concatenate((crack, np.array([crack[0]])))
crack = crack.astype(np.int32)

#Lấy tập chuẩn G trong S
standard_set = hh_skeletonize(original_image,crack)

#Tính độ rộng tại hữu hạn các điểm trong G
standard_set_length=len(standard_set)-1
Num_p=int(standard_set_length*0.1)
Angle_p=4
standard_width_lenght=0
# standard_width = None
minU = None
minV = None

error_set = []


# p = standard_set[int(standard_set_length*k)]
p = [523, 245]

# Cập nhật code chỗ này
p_width, min_u, min_v, visibility_set, test_set = width_crack(p, crack)
######

u = [min_u[0], min_v[0]]
v = [min_u[1], min_v[1]]

if p_width>standard_width_lenght:
    standard_width_lenght = p_width
    minU = u
    minV = v
    #####

visibility_set.append(visibility_set[0])
a, b = zip(*visibility_set)

# Vẽ đoạn thẳng qua u, v
plt.plot(u, v, color='yellow', linewidth=0.5)

# Vẽ điểm p của visibility polygon
plt.plot(p[0], p[1], marker='o')

# plot visibility_polygon
plt.plot(a, b, color='blue', linewidth=0.8)

# plot 2 points in M(P,p)
plt.plot(min_u[0], min_u[1], marker='o')
plt.plot(min_v[0], min_v[1], marker='o')

# Plot 2 edge that contain u and v:
u1 = [test_set[0][0][0], test_set[0][1][0]]
v1 = [test_set[0][0][1], test_set[0][1][1]]
u2 = [test_set[1][0][0], test_set[1][1][0]]
v2 = [test_set[1][0][1], test_set[1][1][1]]
plt.plot(u1, v1, color='red', linewidth=1.5)
plt.plot(u2, v2, color='red', linewidth=1.5)

# Cập nhật code chỗ này
# plt.plot(standard_width[:, 0], standard_width[:, 1], color='red', linewidth=1, label=str(round(standard_width_lenght*k_width,3))+' mm')
plt.plot(minU, minV, color='green', linewidth=1, label=str(round(standard_width_lenght*k_width,3))+' mm')
#####

plt.scatter(standard_set[:, 0], standard_set[:, 1], s=0.05, color='yellow')
# plt.fill(crack[:, 0], crack[:, 1], color='yellow', alpha=0.2)
plt.plot(crack[:, 0], crack[:, 1], color='red', linewidth=0.5)

print("standard_width_lenght = ", standard_width_lenght)
plt.show()