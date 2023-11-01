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

img_name = '7.png'#Tên ảnh

img_path = path_import+img_name
original_image = cv2.imread(img_path)
k_width = 0.011754334835464791

#Nhận diện vết nứt S
model = YOLO("crack_width_YOLOv8.pt")
cracks = model.predict(source=img_path, conf=0.4)[0]

for crack in cracks.masks.xy:
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

    print("crack: ", crack)

    for k in [1/Num_p*i for i in range(Num_p+1)]:
        try:
            p = standard_set[int(standard_set_length*k)]

            # Cập nhật code chỗ này
            # p_segments = hh_width_point_approx(p, crack, Angle_p)
            p_width, min_u, min_v = width_crack(p, crack)
            # p_width = np.sqrt((p_segments.min[0][0]-p_segments.min[1][0])**2+(p_segments.min[0][1]-p_segments.min[1][1])**2)
            ######

            u = [min_u[0], min_v[0]]
            v = [min_u[1], min_v[1]]

            if p_width>standard_width_lenght:
                standard_width_lenght = p_width
                # Cập nhật code chỗ này
                # standard_width = p_segments.min
                minU = u
                minV = v
                #####

            # for seg in p_segments.all:
            #     plt.plot(seg[:, 0], seg[:, 1], color='green', linewidth=0.5)
            
            # plt.plot(p_segments.min[:, 0], p_segments.min[:, 1], color='blue', linewidth=0.2)
            plt.plot(u, v, color='yellow', linewidth=0.2)

            # plt.scatter(p[0], p[1], s=2, c='orange', zorder=2)
        except:
            pass
    
    # print('standard_width = ', standard_width[:, 0])

    # Cập nhật code chỗ này
    # plt.plot(standard_width[:, 0], standard_width[:, 1], color='red', linewidth=1, label=str(round(standard_width_lenght*k_width,3))+' mm')
    plt.plot(minU, minV, color='green', linewidth=1, label=str(round(standard_width_lenght*k_width,3))+' mm')
    #####

    plt.scatter(standard_set[:, 0], standard_set[:, 1], s=0.05, color='yellow')
    # plt.fill(crack[:, 0], crack[:, 1], color='yellow', alpha=0.2)
    plt.plot(crack[:, 0], crack[:, 1], color='red', linewidth=1)
    
plt.imshow(cv2.cvtColor(original_image, cv2.COLOR_BGR2RGB))
plt.axis('off')
plt.legend(loc="upper right",fontsize=15)
plt.savefig('fig-width/'+img_name)
plt.show()