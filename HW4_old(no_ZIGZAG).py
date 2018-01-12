#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from pylab import *

from Tkinter import *
from PIL import ImageTk, Image
import tkMessageBox
import tkFileDialog 
from ttk import Frame, Button, Label, Style

import numpy
import scipy.fftpack

import pickle

import colorsys

#======  搜尋dataset中的圖片數量  =============
import sys

fileList = []
fileSize = 0
folderCount = 0
rootdir = './dataset'

for root, subFolders, files in os.walk(rootdir):
    folderCount += len(subFolders)
    for file in files:
        f = os.path.join(root,file)
        fileSize = fileSize + os.path.getsize(f)
        #print(f)
        fileList.append(f)

print("Total Size is {0} bytes".format(fileSize))
print("Total Files ", len(fileList))
print("Total Folders ", folderCount)


images_quantities = len(fileList)
#=================================================


#images_quantities = 10
all_images_dictionary = {}

#with open('HW3_dictData_kmean500.pickle', 'rb') as handle:
#    all_images_dictionary = pickle.load(handle)




class Tkinter_GUI(Frame):
  
    def __init__(self, parent):
        Frame.__init__(self, parent)   
         
        self.parent = parent
        
        self.initUI()
    
        
    def initUI(self):
        self.parent.title("HW4") 
        self.pack(fill=BOTH, expand=1)

        Button(self, text = "Select File", command = openFile).grid(row=0, column=0, pady=5)
        self.fileName = StringVar()
        Label(self, textvariable=self.fileName).grid(row=0, column=1, columnspan=2, pady=5, sticky=W)
        self.queryImage = []
        self.queryImage.append(Label(self))
        self.queryImage[0].grid(row=0, column=3, pady=5, sticky=W)

        Label(self, text = "Select Mode: ").grid(row=1, column=0, pady=5)
        mode = StringVar(self)
        mode.set("ColorHistogram")
        om = OptionMenu(self, mode, "ColorHistogram", "ColorLayout", "AverageRGB")
        om.grid(row=1, column=1, pady=5, sticky=W)
        
        mode2 = StringVar(self)
        mode2.set("8")
        om2 = OptionMenu(self, mode2, "8", "20", "32", "40", "80")
        om2.grid(row=1, column=2, pady=5, sticky=W)
        
        mode3 = StringVar(self)
        mode3.set("RGB")
        om3 = OptionMenu(self, mode3, "RGB", "HSV")
        om3.grid(row=1, column=3, pady=5, sticky=W)

        Button(self, text = "SEARCH", command = lambda: startSearching(self.fileName.get(),mode.get(),mode2.get(),mode3.get())).grid(row=3, column=0, pady=5)

        
        #用來將top10名的圖片貼到window上
        self.images = []
        self.images.append(Label(self))
        self.images[0].grid(row=5, column=4, padx=0, pady=0)
        #for i in range(number_of_grids_n*number_of_grids_n):
        #    self.images.append(Label(self))
        #    self.images[i].grid(row=(i/number_of_grids_n)+4, column=i%number_of_grids_n+4, padx=0, pady=0)
            
        #用來將top10名的文字貼到window上
        #self.rank_labels = []
        #for i in range(10):
        #    self.rank_labels.append(Label(self))
        #    self.rank_labels[i].grid(row=(i/5)*2+5, column=i%5, padx=5, pady=10)
            


def openFile ():
    fileName = tkFileDialog.askopenfilename(initialdir = "./dataset")
    app.fileName.set(fileName)
    
    image = Image.open(fileName)
    image = image.resize((50, 50), Image.ANTIALIAS) #The (80, 80) is (height, width)
    photo = ImageTk.PhotoImage(image)
    app.queryImage[0].configure(image=photo)
    app.queryImage[0].image = photo

global is_there_pickle
is_there_pickle = False

number_of_grids_n = 20

def startSearching (fileName, mode, mode2, mode3):
    #print "Your Code Here."
    global all_images_dictionary
    all_images_dictionary.clear()
    global is_there_pickle
    is_there_pickle = False
    
    #read n from menu2
    if mode2 == "8":
        number_of_grids_n = 8
    elif mode2 == "20":
        number_of_grids_n = 20
    elif mode2 == "32":
        number_of_grids_n = 32
    elif mode2 == "40":
        number_of_grids_n = 40
    elif mode2 == "80":
        number_of_grids_n = 80
        
        
        
    if mode == "ColorHistogram":
        print "Q1."
        Query_fileName = fileName
        
        Query_im = Image.open(Query_fileName).resize((640, 480)).convert('RGB') #讀入圖片
        Query_pixel = Query_im.load()
        Query_width, Query_height = Query_im.size
        
        
        #確認是否已經有pickle檔
        pickle_path = 'HW4_ColorHistogram_' + str(number_of_grids_n) + 'x' + str(number_of_grids_n) + '.pickle'
        if os.path.exists(pickle_path):
                is_there_pickle = True
                with open(pickle_path, 'rb') as handle:
                    all_images_dictionary = pickle.load(handle)
        
        #讀入所有jpg的Color Histogram(只有第一次按Search時需要)
        if (is_there_pickle == False):
            for i in range(0, images_quantities):
                temp_str = str(i)
                temp_str_fillZero = temp_str.zfill(4)
                i_th_fileName = "./dataset/ukbench0"+ temp_str_fillZero +".jpg"
                
                if i not in all_images_dictionary:
                    all_images_dictionary[i] = {}
                    all_images_dictionary[i]["fileName"] = i_th_fileName
                    
                all_images_dictionary[i]["color_histogram_R_list"] = list()
                all_images_dictionary[i]["color_histogram_G_list"] = list()
                all_images_dictionary[i]["color_histogram_B_list"] = list()
                
                for j in range(0, 256):
                    all_images_dictionary[i]["color_histogram_R_list"].append(0)
                    all_images_dictionary[i]["color_histogram_G_list"].append(0)
                    all_images_dictionary[i]["color_histogram_B_list"].append(0)
                
                #HSV
                all_images_dictionary[i]["color_histogram_H_list"] = list()
                all_images_dictionary[i]["color_histogram_S_list"] = list()
                all_images_dictionary[i]["color_histogram_V_list"] = list()
                
                for j in range(0, 360):
                    all_images_dictionary[i]["color_histogram_H_list"].append(0)
                    all_images_dictionary[i]["color_histogram_S_list"].append(0)
                    all_images_dictionary[i]["color_histogram_V_list"].append(0)
                
                
                #讀入Color Histogram    
                im = Image.open(i_th_fileName).convert('RGB').resize((Query_width/number_of_grids_n, Query_height/number_of_grids_n)) #讀入圖片
                pixel = im.load()
                width, height = im.size

                for x in xrange(width):    
                    for y in xrange(height):
                        R, G, B = pixel[x,y]
                        all_images_dictionary[i]["color_histogram_R_list"][R]+=1
                        all_images_dictionary[i]["color_histogram_G_list"][G]+=1
                        all_images_dictionary[i]["color_histogram_B_list"][B]+=1
                        H, S, V = colorsys.rgb_to_hsv(R/255., G/255., B/255.)
                        all_images_dictionary[i]["color_histogram_H_list"][int(H*360)]+=1
                        #all_images_dictionary[i]["color_histogram_H_list"][(int(H*360))/3]+=1
                        all_images_dictionary[i]["color_histogram_S_list"][int(S*100)]+=1
                        all_images_dictionary[i]["color_histogram_V_list"][int(V*100)]+=1
                        
                        
            #輸出pickle
            output_pickle_path = 'HW4_ColorHistogram_' + str(number_of_grids_n) + 'x' + str(number_of_grids_n) + '.pickle'
            with open(output_pickle_path, 'wb') as handle:
                pickle.dump(all_images_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
        #-------------------------------------------
        
        #------- mosaic  ---------------------------
        
        mosaic_images_dictionary = {}
        
        im = Image.open(Query_fileName).resize((640, 480)).convert('RGB') #讀入圖片
        pixel = im.load()
        width, height = im.size

        for x in xrange(width):    
            for y in xrange(height):
                if(((width%number_of_grids_n) == 0) and ((height%number_of_grids_n) == 0)):
                    curr_block_key = "block" + str((x/(width/number_of_grids_n))).zfill(4)+ str((y/(height/number_of_grids_n))).zfill(4)
                elif(((width%number_of_grids_n) == 0) and ((height%number_of_grids_n) != 0)):
                    curr_block_key = "block" + str((x/(width/number_of_grids_n))).zfill(4)+ str((y/((height/number_of_grids_n)+1))).zfill(4)
                elif(((width%number_of_grids_n) != 0) and ((height%number_of_grids_n) == 0)):
                    curr_block_key = "block" + str((x/((width/number_of_grids_n)+1))).zfill(4)+ str((y/(height/number_of_grids_n))).zfill(4)
                else:
                    curr_block_key = "block" + str((x/((width/number_of_grids_n)+1))).zfill(4)+ str((y/((height/number_of_grids_n)+1))).zfill(4)
                
                if curr_block_key not in mosaic_images_dictionary:
                    mosaic_images_dictionary[curr_block_key] = {}
                    mosaic_images_dictionary[curr_block_key]["color_histogram_R_list"] = list()
                    mosaic_images_dictionary[curr_block_key]["color_histogram_G_list"] = list()
                    mosaic_images_dictionary[curr_block_key]["color_histogram_B_list"] = list()
                    
                    for j in range(0, 256):
                        mosaic_images_dictionary[curr_block_key]["color_histogram_R_list"].append(0)
                        mosaic_images_dictionary[curr_block_key]["color_histogram_G_list"].append(0)
                        mosaic_images_dictionary[curr_block_key]["color_histogram_B_list"].append(0)
                        
                    #HSV
                    mosaic_images_dictionary[curr_block_key]["color_histogram_H_list"] = list()
                    mosaic_images_dictionary[curr_block_key]["color_histogram_S_list"] = list()
                    mosaic_images_dictionary[curr_block_key]["color_histogram_V_list"] = list()
                    
                    for j in range(0, 360):
                        mosaic_images_dictionary[curr_block_key]["color_histogram_H_list"].append(0)
                        mosaic_images_dictionary[curr_block_key]["color_histogram_S_list"].append(0)
                        mosaic_images_dictionary[curr_block_key]["color_histogram_V_list"].append(0)
                    
                
                R, G, B = pixel[x,y]
                mosaic_images_dictionary[curr_block_key]["color_histogram_R_list"][R]+=1
                mosaic_images_dictionary[curr_block_key]["color_histogram_G_list"][G]+=1
                mosaic_images_dictionary[curr_block_key]["color_histogram_B_list"][B]+=1
                
                H, S, V = colorsys.rgb_to_hsv(R/255., G/255., B/255.)
                mosaic_images_dictionary[curr_block_key]["color_histogram_H_list"][int(H*360)]+=1
                #mosaic_images_dictionary[curr_block_key]["color_histogram_H_list"][(int(H*360))/3]+=1
                mosaic_images_dictionary[curr_block_key]["color_histogram_S_list"][int(S*100)]+=1
                mosaic_images_dictionary[curr_block_key]["color_histogram_V_list"][int(V*100)]+=1
        #-------------------------------------------
        print "calculate_distance"
        if mode3 == "RGB":
            distance_QueryWithImages = list()
            for x in xrange(number_of_grids_n):    
                for y in xrange(number_of_grids_n):
                    curr_block_key = "block" + str(x).zfill(4)+ str(y).zfill(4)
                    for i in range(0, images_quantities):
                        distance_QueryWithImages.append(list())
                        ith_distance = 0
                        ith_R_distance = 0
                        ith_G_distance = 0
                        ith_B_distance = 0
                        for j in range(0, 256):
                            ith_R_distance += ((mosaic_images_dictionary[curr_block_key]["color_histogram_R_list"][j] - all_images_dictionary[i]["color_histogram_R_list"][j])**2)
                            ith_G_distance += ((mosaic_images_dictionary[curr_block_key]["color_histogram_G_list"][j] - all_images_dictionary[i]["color_histogram_G_list"][j])**2)
                            ith_B_distance += ((mosaic_images_dictionary[curr_block_key]["color_histogram_B_list"][j] - all_images_dictionary[i]["color_histogram_B_list"][j])**2)
                         
                        #print ith_distance
                        ith_R_distance = ith_R_distance ** (0.5)
                        ith_G_distance = ith_G_distance ** (0.5)
                        ith_B_distance = ith_B_distance ** (0.5)
                        ith_distance = ith_R_distance + ith_G_distance + ith_B_distance
                        distance_QueryWithImages[i].append(i)
                        distance_QueryWithImages[i].append(ith_distance)
                        #print curr_block_key,i,ith_distance
                        
                    distance_QueryWithImages = sorted(distance_QueryWithImages,key=lambda l:l[1], reverse=False) #讓list照distance由小到大排序
                    mosaic_images_dictionary[curr_block_key]["closest_fileName"] = all_images_dictionary[distance_QueryWithImages[0][0]]["fileName"]
                    del distance_QueryWithImages[:]
                    
        else:
            distance_QueryWithImages = list()
            for x in xrange(number_of_grids_n):    
                for y in xrange(number_of_grids_n):
                    curr_block_key = "block" + str(x).zfill(4)+ str(y).zfill(4)
                    for i in range(0, images_quantities):
                        distance_QueryWithImages.append(list())
                        ith_distance = 0
                        ith_H_distance = 0
                        ith_S_distance = 0
                        ith_V_distance = 0
                        for j in range(0, 360):
                        #for j in range(0, 36):
                            ith_H_distance += ((mosaic_images_dictionary[curr_block_key]["color_histogram_H_list"][j] - all_images_dictionary[i]["color_histogram_H_list"][j])**2)
                            ith_S_distance += ((mosaic_images_dictionary[curr_block_key]["color_histogram_S_list"][j] - all_images_dictionary[i]["color_histogram_S_list"][j])**2)
                            ith_V_distance += ((mosaic_images_dictionary[curr_block_key]["color_histogram_V_list"][j] - all_images_dictionary[i]["color_histogram_V_list"][j])**2)
 
                        #print ith_distance
                        ith_distance = ith_H_distance + ith_S_distance + ith_V_distance
                        distance_QueryWithImages[i].append(i)
                        distance_QueryWithImages[i].append(ith_distance)
                        #print curr_block_key,i,ith_distance
                        
                    distance_QueryWithImages = sorted(distance_QueryWithImages,key=lambda l:l[1], reverse=False) #讓list照distance由小到大排序
                    mosaic_images_dictionary[curr_block_key]["closest_fileName"] = all_images_dictionary[distance_QueryWithImages[0][0]]["fileName"]
                    del distance_QueryWithImages[:]
                    
        #-------------------------------------------
        print "save_result"
        result_im = Image.open(Query_fileName).resize((640, 480)).convert('RGB') #讀入圖片
        result_pixel = result_im.load()
        result_width, result_height = result_im.size
        
        for x in xrange(number_of_grids_n):    
            for y in xrange(number_of_grids_n):
                curr_block_key = "block" + str(x).zfill(4)+ str(y).zfill(4)
        
                im_tile = Image.open(mosaic_images_dictionary[curr_block_key]["closest_fileName"])
                im_tile = im_tile.resize((Query_width/number_of_grids_n, Query_height/number_of_grids_n))
                pixel_tile = im_tile.load()
                width_tile, height_tile = im_tile.size
                
                for x_tile in xrange(width_tile):    
                    for y_tile in xrange(height_tile):
                        R, G, B = pixel_tile[x_tile,y_tile]
                        result_pixel[x*(Query_width/number_of_grids_n)+x_tile,y*(Query_height/number_of_grids_n)+y_tile] = (R, G, B)
                        
                        
        photo = ImageTk.PhotoImage(result_im)
        app.images[0].configure(image=photo)
        app.images[0].image = photo
        """
        for x in xrange(number_of_grids_n):    
            for y in xrange(number_of_grids_n):
                curr_block_key = "block" + str(x).zfill(4)+ str(y).zfill(4)
                #print curr_block_key
        #app.rank_labels[i].configure(text="Rank " + str(i) + " is number " + str(distance_QueryWithImages[i][0]) + ", distance is " + str(int(distance_QueryWithImages[i][1])))
                photo = ImageTk.PhotoImage(Image.open(mosaic_images_dictionary[curr_block_key]["closest_fileName"]).resize((width/number_of_grids_n, height/number_of_grids_n)))
                app.images[x*number_of_grids_n + y].configure(image=photo)
                app.images[x*number_of_grids_n + y].image = photo
        """

                           
    elif mode == "ColorLayout":
        print "Q2."
        
        if number_of_grids_n > 40: #超過20的話目前會有問題?
            print "error grids_n"
            return
            
        Q2_dict = {} 
        Query_fileName = fileName
        
        Query_im = Image.open(Query_fileName).resize((640, 480)).convert('RGB') #讀入圖片
        Query_pixel = Query_im.load()
        Query_width, Query_height = Query_im.size
        
        
        #確認是否已經有pickle檔
        pickle_path = 'HW4_ColorLayout_' + str(number_of_grids_n) + 'x' + str(number_of_grids_n) + '.pickle'
        if os.path.exists(pickle_path):
                is_there_pickle = True
                with open(pickle_path, 'rb') as handle:
                    all_images_dictionary = pickle.load(handle)
        
        #計算所有jpg的DCT(只有第一次按Search時需要)
        if (is_there_pickle == False):
            for i in range(0, images_quantities):
                Q2_dict.clear()#Q2_dict用來計算每次輸入的image的DCT_Y, DCT_Cb, DCT_Cr，每次匯入新圖時清空
                temp_str = str(i)
                temp_str_fillZero = temp_str.zfill(4)
                i_th_fileName = "./dataset/ukbench0"+ temp_str_fillZero +".jpg"
                
                if i not in all_images_dictionary:
                    all_images_dictionary[i] = {}
                    all_images_dictionary[i]["fileName"] = i_th_fileName
                    
                all_images_dictionary[i]["DCT_Y_list"] = list()
                all_images_dictionary[i]["DCT_Cb_list"] = list()
                all_images_dictionary[i]["DCT_Cr_list"] = list()

                
                #開始計算DCT 
                im = Image.open(i_th_fileName).convert('RGB').resize((Query_width/number_of_grids_n, Query_height/number_of_grids_n)) #讀入圖片
                pixel = im.load()
                width, height = im.size
                
                for x in xrange(width):    
                    for y in xrange(height):
                        if(((width%8) == 0) and ((height%8) == 0)):
                            curr_block_key = "block" + str((x/(width/8)))+ str((y/(height/8)))
                        elif(((width%8) == 0) and ((height%8) != 0)):
                            if((y/((height/8)+1)) < (height%8)):
                                curr_block_key = "block" + str((x/(width/8)))+ str((y/((height/8)+1)))
                            else:
                                curr_block_key = "block" + str((x/(width/8)))+ str((height%8)+((y-(((height/8)+1)*(height%8)))/(height/8)))
                        elif(((width%8) != 0) and ((height%8) == 0)):
                            if((x/((width/8)+1)) < (width%8)):
                                curr_block_key = "block" + str((x/((width/8)+1)))+ str((y/(height/8)))
                            else:
                                curr_block_key = "block" + str((width%8)+((x-(((width/8)+1)*(width%8)))/(width/8)))+ str((y/(height/8)))
                        else:
                            if((x/((width/8)+1)) < (width%8)):
                                curr_block_key = "block" + str((x/((width/8)+1)))
                            else:
                                curr_block_key = "block" + str((width%8)+((x-(((width/8)+1)*(width%8)))/(width/8)))
                            
                            if((y/((height/8)+1)) < (height%8)):
                                curr_block_key = curr_block_key + str((y/((height/8)+1)))
                            else:
                                curr_block_key = curr_block_key + str((height%8)+((y-(((height/8)+1)*(height%8)))/(height/8)))
                            
                            
                        #print width, height, x,y,curr_block_key
                        if curr_block_key not in Q2_dict:
                            Q2_dict[curr_block_key] = list()
                            Q2_dict[curr_block_key].append(0) #Q2_dict[curr_block_key][0] 儲存pixels個數
                            Q2_dict[curr_block_key].append(0) #Q2_dict[curr_block_key][1] Representative Color R
                            Q2_dict[curr_block_key].append(0) #Q2_dict[curr_block_key][2] Representative Color G
                            Q2_dict[curr_block_key].append(0) #Q2_dict[curr_block_key][3] Representative Color B
                        
                        R, G, B = pixel[x,y]
                        Q2_dict[curr_block_key][0] += 1
                        Q2_dict[curr_block_key][1] += R
                        Q2_dict[curr_block_key][2] += G
                        Q2_dict[curr_block_key][3] += B
                    
                    
                #前一個雙層for只有做加總，這邊將加總完的R,G,B除以該block的總pixels數量        
                for j in Q2_dict:
                    Q2_dict[j][1] = Q2_dict[j][1]/Q2_dict[j][0]
                    Q2_dict[j][2] = Q2_dict[j][2]/Q2_dict[j][0]
                    Q2_dict[j][3] = Q2_dict[j][3]/Q2_dict[j][0]
                    
                  
                 
                #步驟二：把算完的平均值存進8x8的image
                im_8x8 = Image.open(fileName)
                im_8x8 = im_8x8.resize((8, 8))
                pixel_8x8 = im_8x8.load()
                width_8x8, height_8x8 = im_8x8.size
                
                for x in xrange(width_8x8):    
                    for y in xrange(height_8x8):
                        curr_block_key = "block" + str(x)+ str(y)
                        pixel_8x8[x,y] = (Q2_dict[curr_block_key][1],Q2_dict[curr_block_key][2],Q2_dict[curr_block_key][3])
                
                
                
                #步驟三：把8x8的RGB轉成8x8的YCbCr
                im_ycbcr = im_8x8.convert('YCbCr')
                pixel_ycbcr = im_ycbcr.load()
                width_ycbcr, height_ycbcr = im_ycbcr.size
                
                #用來print出8x8的im_ycbcr的每個pixel的值
                #for x in xrange(width_ycbcr):    
                #    for y in xrange(height_ycbcr):
                #        Y, Cb, Cr = pixel_ycbcr[x,y]
                #        print Y, Cb, Cr
                        
                
                
                #步驟四：對8x8的YCbCr image(im_ycbcr)做DCT運算，Y、Cb、Cr分別運算
                dctSize = im_ycbcr.size[0]
                
                pixels_Y = numpy.array(im_ycbcr.getdata(0), dtype=numpy.float).reshape((dctSize, dctSize))
                pixels_Cb = numpy.array(im_ycbcr.getdata(1), dtype=numpy.float).reshape((dctSize, dctSize))
                pixels_Cr = numpy.array(im_ycbcr.getdata(2), dtype=numpy.float).reshape((dctSize, dctSize))


                # perform 2-dimensional DCT (discrete cosine transform):
                DCT_Y = scipy.fftpack.dct(scipy.fftpack.dct(pixels_Y.T, norm="ortho").T, norm="ortho")
                DCT_Cb = scipy.fftpack.dct(scipy.fftpack.dct(pixels_Cb.T, norm="ortho").T, norm="ortho")
                DCT_Cr = scipy.fftpack.dct(scipy.fftpack.dct(pixels_Cr.T, norm="ortho").T, norm="ortho")
                
                for x in xrange(width_ycbcr):    
                    for y in xrange(height_ycbcr):
                        all_images_dictionary[i]["DCT_Y_list"].append(DCT_Y[x][y])
                        all_images_dictionary[i]["DCT_Cb_list"].append(DCT_Cb[x][y])
                        all_images_dictionary[i]["DCT_Cr_list"].append(DCT_Cr[x][y])
            
            
            #輸出pickle
            output_pickle_path = 'HW4_ColorLayout_' + str(number_of_grids_n) + 'x' + str(number_of_grids_n) + '.pickle'
            with open(output_pickle_path, 'wb') as handle:
                pickle.dump(all_images_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
        #-------------------------------------------
        
        #------- mosaic  ---------------------------
        
        mosaic_images_dictionary = {}
        
        for block_x in xrange(number_of_grids_n):    
            for block_y in xrange(number_of_grids_n):
                Q2_dict.clear()#Q2_dict用來計算每次輸入的image的DCT_Y, DCT_Cb, DCT_Cr，每次匯入新圖時清空
                mosaic_curr_block_key = "block" + str(block_x).zfill(4)+ str(block_y).zfill(4)
                
                if mosaic_curr_block_key not in mosaic_images_dictionary:
                    mosaic_images_dictionary[mosaic_curr_block_key] = {}
                    
                mosaic_images_dictionary[mosaic_curr_block_key]["DCT_Y_list"] = list()
                mosaic_images_dictionary[mosaic_curr_block_key]["DCT_Cb_list"] = list()
                mosaic_images_dictionary[mosaic_curr_block_key]["DCT_Cr_list"] = list()
                    
                
                #開始計算DCT 
                im = Image.open(Query_fileName).convert('RGB').resize((Query_width/number_of_grids_n, Query_height/number_of_grids_n)) #讀入圖片
                pixel = im.load()
                width, height = im.size
                
                for x in xrange(width):    
                    for y in xrange(height):
                        R, G, B = Query_pixel[block_x*(Query_width/number_of_grids_n)+x,block_y*(Query_height/number_of_grids_n)+y]
                        pixel[x, y] = (R, G, B)
                   

                for x in xrange(width):    
                    for y in xrange(height):
                        if(((width%8) == 0) and ((height%8) == 0)):
                            curr_block_key = "block" + str((x/(width/8)))+ str((y/(height/8)))
                        elif(((width%8) == 0) and ((height%8) != 0)):
                            if((y/((height/8)+1)) < (height%8)):
                                curr_block_key = "block" + str((x/(width/8)))+ str((y/((height/8)+1)))
                            else:
                                curr_block_key = "block" + str((x/(width/8)))+ str((height%8)+((y-(((height/8)+1)*(height%8)))/(height/8)))
                        elif(((width%8) != 0) and ((height%8) == 0)):
                            if((x/((width/8)+1)) < (width%8)):
                                curr_block_key = "block" + str((x/((width/8)+1)))+ str((y/(height/8)))
                            else:
                                curr_block_key = "block" + str((width%8)+((x-(((width/8)+1)*(width%8)))/(width/8)))+ str((y/(height/8)))
                        else:
                            if((x/((width/8)+1)) < (width%8)):
                                curr_block_key = "block" + str((x/((width/8)+1)))
                            else:
                                curr_block_key = "block" + str((width%8)+((x-(((width/8)+1)*(width%8)))/(width/8)))
                            
                            if((y/((height/8)+1)) < (height%8)):
                                curr_block_key = curr_block_key + str((y/((height/8)+1)))
                            else:
                                curr_block_key = curr_block_key + str((height%8)+((y-(((height/8)+1)*(height%8)))/(height/8)))
                        
                        if curr_block_key not in Q2_dict:
                            Q2_dict[curr_block_key] = list()
                            Q2_dict[curr_block_key].append(0) #Q2_dict[curr_block_key][0] 儲存pixels個數
                            Q2_dict[curr_block_key].append(0) #Q2_dict[curr_block_key][1] Representative Color R
                            Q2_dict[curr_block_key].append(0) #Q2_dict[curr_block_key][2] Representative Color G
                            Q2_dict[curr_block_key].append(0) #Q2_dict[curr_block_key][3] Representative Color B
                        
                        R, G, B = pixel[x,y]
                        Q2_dict[curr_block_key][0] += 1
                        Q2_dict[curr_block_key][1] += R
                        Q2_dict[curr_block_key][2] += G
                        Q2_dict[curr_block_key][3] += B
                    
                    
                #前一個雙層for只有做加總，這邊將加總完的R,G,B除以該block的總pixels數量        
                for j in Q2_dict:
                    Q2_dict[j][1] = Q2_dict[j][1]/Q2_dict[j][0]
                    Q2_dict[j][2] = Q2_dict[j][2]/Q2_dict[j][0]
                    Q2_dict[j][3] = Q2_dict[j][3]/Q2_dict[j][0]
                    
                  
                 
                #步驟二：把算完的平均值存進8x8的image
                im_8x8 = Image.open(fileName)
                im_8x8 = im_8x8.resize((8, 8))
                pixel_8x8 = im_8x8.load()
                width_8x8, height_8x8 = im_8x8.size
                
                for x in xrange(width_8x8):    
                    for y in xrange(height_8x8):
                        curr_block_key = "block" + str(x)+ str(y)
                        pixel_8x8[x,y] = (Q2_dict[curr_block_key][1],Q2_dict[curr_block_key][2],Q2_dict[curr_block_key][3])
                
                
                
                #步驟三：把8x8的RGB轉成8x8的YCbCr
                im_ycbcr = im_8x8.convert('YCbCr')
                pixel_ycbcr = im_ycbcr.load()
                width_ycbcr, height_ycbcr = im_ycbcr.size
                
                #用來print出8x8的im_ycbcr的每個pixel的值
                #for x in xrange(width_ycbcr):    
                #    for y in xrange(height_ycbcr):
                #        Y, Cb, Cr = pixel_ycbcr[x,y]
                #        print Y, Cb, Cr
                        
                
                
                #步驟四：對8x8的YCbCr image(im_ycbcr)做DCT運算，Y、Cb、Cr分別運算
                dctSize = im_ycbcr.size[0]
                
                pixels_Y = numpy.array(im_ycbcr.getdata(0), dtype=numpy.float).reshape((dctSize, dctSize))
                pixels_Cb = numpy.array(im_ycbcr.getdata(1), dtype=numpy.float).reshape((dctSize, dctSize))
                pixels_Cr = numpy.array(im_ycbcr.getdata(2), dtype=numpy.float).reshape((dctSize, dctSize))


                # perform 2-dimensional DCT (discrete cosine transform):
                DCT_Y = scipy.fftpack.dct(scipy.fftpack.dct(pixels_Y.T, norm="ortho").T, norm="ortho")
                DCT_Cb = scipy.fftpack.dct(scipy.fftpack.dct(pixels_Cb.T, norm="ortho").T, norm="ortho")
                DCT_Cr = scipy.fftpack.dct(scipy.fftpack.dct(pixels_Cr.T, norm="ortho").T, norm="ortho")
                
                for x in xrange(width_ycbcr):    
                    for y in xrange(height_ycbcr):
                        mosaic_images_dictionary[mosaic_curr_block_key]["DCT_Y_list"].append(DCT_Y[x][y])
                        mosaic_images_dictionary[mosaic_curr_block_key]["DCT_Cb_list"].append(DCT_Cb[x][y])
                        mosaic_images_dictionary[mosaic_curr_block_key]["DCT_Cr_list"].append(DCT_Cr[x][y])
                
        
        #-------------------------------------------
        #計算Query跟所有image個別的distance(利用Eucildean)
        distance_QueryWithImages = list()
        
        for x in xrange(number_of_grids_n):    
            for y in xrange(number_of_grids_n):
                curr_block_key = "block" + str(x).zfill(4)+ str(y).zfill(4)
                for i in range(0, images_quantities):
                    distance_QueryWithImages.append(list())
                    ith_distance = 0
                    ith_Y_distance = 0
                    ith_Cb_distance = 0
                    ith_Cr_distance = 0
                    for j in range(0, 64):
                        ith_Y_distance += ((mosaic_images_dictionary[curr_block_key]["DCT_Y_list"][j] - all_images_dictionary[i]["DCT_Y_list"][j])**2)
                        ith_Cb_distance += ((mosaic_images_dictionary[curr_block_key]["DCT_Cb_list"][j] - all_images_dictionary[i]["DCT_Cb_list"][j])**2)
                        ith_Cr_distance += ((mosaic_images_dictionary[curr_block_key]["DCT_Cr_list"][j] - all_images_dictionary[i]["DCT_Cr_list"][j])**2)
                     
                    #print ith_distance
                    ith_Y_distance = ith_Y_distance ** (0.5)
                    ith_Cb_distance = ith_Cb_distance ** (0.5)
                    ith_Cr_distance = ith_Cr_distance ** (0.5)
                    ith_distance = ith_Y_distance + ith_Cb_distance + ith_Cr_distance
                    distance_QueryWithImages[i].append(i)
                    distance_QueryWithImages[i].append(ith_distance)
                    #print curr_block_key,i,ith_distance
                    
                distance_QueryWithImages = sorted(distance_QueryWithImages,key=lambda l:l[1], reverse=False) #讓list照distance由小到大排序
                mosaic_images_dictionary[curr_block_key]["closest_fileName"] = all_images_dictionary[distance_QueryWithImages[0][0]]["fileName"]
                del distance_QueryWithImages[:]
                    
        #----------------------------------------------------------------------------------------
        
        result_im = Image.open(Query_fileName).resize((640, 480)).convert('RGB') #讀入圖片
        result_pixel = result_im.load()
        result_width, result_height = result_im.size
        
        for x in xrange(number_of_grids_n):    
            for y in xrange(number_of_grids_n):
                curr_block_key = "block" + str(x).zfill(4)+ str(y).zfill(4)
        
                im_tile = Image.open(mosaic_images_dictionary[curr_block_key]["closest_fileName"])
                im_tile = im_tile.resize((Query_width/number_of_grids_n, Query_height/number_of_grids_n))
                pixel_tile = im_tile.load()
                width_tile, height_tile = im_tile.size
                
                for x_tile in xrange(width_tile):    
                    for y_tile in xrange(height_tile):
                        R, G, B = pixel_tile[x_tile,y_tile]
                        result_pixel[x*(Query_width/number_of_grids_n)+x_tile,y*(Query_height/number_of_grids_n)+y_tile] = (R, G, B)
                        
                        
        photo = ImageTk.PhotoImage(result_im)
        app.images[0].configure(image=photo)
        app.images[0].image = photo
             
        
    elif mode == "AverageRGB":
        print "Q3."
        
    else:
        print "ELSE."


if __name__ == '__main__':
    root = Tk()
    size = 220, 220

    app = Tkinter_GUI(root)
    #root.geometry("1024x720")
    root.geometry("1280x820")
    root.mainloop()
