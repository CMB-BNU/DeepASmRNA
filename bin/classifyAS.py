# TensorFlow and tf.keras
import tensorflow as tf
from tensorflow import keras
# Helper libraries
import numpy as np
#import matplotlib.pyplot as plt
import os
import sys 
print(tf.__version__)
#os.environ["CUDA_VISIBLE_DEVICES"] = "0"
#import time
import random
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Arguments for use AttentionCNN model predict Alternative Splice event")
    parser.add_argument(dest='input', help='input sequence of AS event, a csv file contains 4 columns, the 1st columns for seq id, 2nd for upstream seq, 3rd for alternative seq, 4th for downstream seq. 2nd and 4th should in length 50', type=argparse.FileType('r'))
    parser.add_argument('-m', help='model for specie, choose from [arabidopsis, human, rice, fine_tune]', default="human", choices=['arabidopsis', 'human', 'rice', 'fine_tune'])
    parser.add_argument('-o', help='output file name', default="AttCnn.result.txt")
    parser.add_argument('-ft', help='input sequence and label of fine tune dataset, a csv file contains 4 columns, the 1st columns for label, in ["SE", "RI", "A3", "A5"], 2nd for upstream seq, 3rd for alternative seq, 4th for downstream seq. 2nd and 4th should in length 50', default="none")
    args = parser.parse_args()

    M = args.m
    data_file = args.input.name
    if args.ft != "none":
        print("start fine tune training!")
        import fine_tune
        fine_tune.fine_tune(args.ft, M)
        M = "fine_tune"
    max_len = 50

    def label_trans(array, array_dict):
        new_array = []
        for i in array:
            new_array.append(array_dict[i])
        return new_array
    
    ### input 转化成one—hot，这里把二倍体转化成4列，理论上损失了phase信息，用八列的one-hot可以保留phase信息
    def to_one_hot_coding(string):
        tmp_dict = [0,0,0,0]
        for i in string:
            if i == 'A' or i == 'a':
                yield [1,0,0,0,0]
            elif i == "T" or i == 't':
                yield [0,1,0,0,0]
            elif i == "C" or i == "c":
                yield [0,0,1,0,0]
            elif i == "G" or i == "g":
                yield [0,0,0,1,0]
            elif i == ",":
                yield [0,0,0,0,1]
            else:
                yield [0,0,0,0,0]

    labels_list = {'SE': 0, 'A3': 1, 'A5': 2, 'RI': 3} 
    reverse_labels_list = {}
    for i, j in labels_list.items():
        reverse_labels_list[j] = i

    with open(data_file, 'r') as F:
        print("loading data from {} >>>".format(data_file))
        lines = F.readlines()
    print("data loaded from {} !".format(data_file))

    tmp_id = []
    tmp_data = []

    for line in lines:
        tmp = line.strip('\n').split(',')
        tmp_id.append(tmp[0])
        if len(tmp[2]) >= max_len:
            tmp_start = tmp[2][:max_len]
            tmp_end = tmp[2][-max_len:]
        else:
            tmp_start = tmp[2]+"0"*(max_len - len(tmp[2]))
            tmp_end = "0"*(max_len - len(tmp[2]))+tmp[2]  
        if len(tmp[1]) > max_len:
            tmp[1] = tmp[1][-max_len:]
        elif len(tmp[1]) < max_len:
            #print("yes")
            tmp[1] = "0"*(max_len-len(tmp[1])) + tmp[1]
        if len(tmp[3]) > max_len:
            tmp[3] = tmp[3][:max_len]
        elif len(tmp[3])<max_len:
            tmp[3] = tmp[3] + "0"*(max_len-len(tmp[3]))

        tmp_data.append(','.join((tmp[1],tmp_start))+','.join((tmp_end,tmp[3])))

    data0 = []
    for i in tmp_data:
        tmp = []
        for j in to_one_hot_coding(i):
            tmp.append(j)
        data0.append(tmp)
    tmp_data = np.array(data0).astype("float32")
    print(tmp_data.shape)
    
    test_set = tmp_data

    test_ds = tf.data.Dataset.from_tensor_slices(test_set).batch(100)

    class SelfAttention(keras.Model):
        def __init__(self, filters, kernel_size, activation="relu", padding="valid"):
            super(SelfAttention, self).__init__()
            self.filters = filters
            self.kernel_size = kernel_size
            self.activation = activation
            self.padding = padding
            self.conv_k = keras.layers.Conv1D(filters=self.filters, kernel_size=self.kernel_size, activation=self.activation, padding = self.padding)
            self.conv_q = keras.layers.Conv1D(filters=self.filters, kernel_size=self.kernel_size, activation=self.activation, padding = self.padding)
            self.conv_v = keras.layers.Conv1D(filters=self.filters, kernel_size=self.kernel_size, activation=self.activation, padding = self.padding)
            self.att1 = keras.layers.Attention()
            self.concat = keras.layers.Concatenate()
        def call(self, x):
            x_k = self.conv_k(x)
            x_q = self.conv_q(x)
            x_v = self.conv_v(x)
            y = self.att1([x_q, x_v, x_k])
            x = self.concat([x_q, y])
            return x
            
    class AttentionCNN(keras.Model):
        def __init__(self):
            super(AttentionCNN, self).__init__()
            self.SelfAtt1 = SelfAttention(filters=30, kernel_size=12)
            self.SelfAtt2 = SelfAttention(filters=100,kernel_size=12)
            self.pool = keras.layers.MaxPool1D(pool_size=2, strides=2)
            self.conv = keras.layers.Conv1D(filters=200, kernel_size=42)
            self.flatten = keras.layers.Flatten()
            self.d1 = keras.layers.Dense(64, activation='relu')
            self.dropout = keras.layers.Dropout(0.6)
            self.d2 = keras.layers.Dense(4, activation='softmax')
    
        def call(self, x):
            x = self.SelfAtt1(x)
            x = self.pool(x)
            x = self.SelfAtt2(x)
            x = self.pool(x)
            x = self.conv(x)
            x = self.flatten(x)
            x = self.d1(x)
            x = self.dropout(x)
            return self.d2(x)
    
    model = AttentionCNN()
    model.load_weights(filepath="{}/final_parameters/CnnAtt_{}_final".format(sys.path[0], M))


    ### model test
    predict_label = []
    true_label = []
    for images in test_ds:
        predict_label.extend(model(images).numpy())
    predict_label = np.array(predict_label)

    print(predict_label.shape)
    print(len(tmp_id))

    with open(args.o , "w") as fw:
        for i in range(len(tmp_id)):
            #fw.write("{}\t{}\n".format(tmp_id[i], ))
            fw.write(tmp_id[i])
            fw.write("\t")
            j = np.argmax(predict_label[i])
            fw.write(reverse_labels_list[j])
            fw.write("\t")
            fw.write(str(predict_label[i][j]))
            fw.write("\n")






