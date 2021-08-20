# TensorFlow and tf.keras
import tensorflow as tf
from tensorflow import keras

# Helper libraries
import numpy as np
import os 
import sys
print(tf.__version__)
import time
import random

def fine_tune(data_file, model_wieght):
    ### loading data
    #### read the csv file and return two list with data and lable

    max_len = 50
    EPOCHS = 20
    
    ### 对数据 转化成numpy多维数组
    #### 随机打乱数据，避免label的集中分布，（测试集中全是同一种标签）
    ### label 转化成id
    #random.seed(20201014)
    
    def filter_lines(lines, labels_list):
        new_lines = []
        for line in lines:
            if line[:2] in labels_list.keys():
                new_lines.append(line)
        return new_lines
    
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
    
    labels = {}
    data = {}
    
    labels_list = {'SE': 0, 'A3': 1, 'A5': 2, 'RI': 3} 
    reverse_labels_list = {}
    for i, j in labels_list.items():
        reverse_labels_list[j] = i
    try:
        F = open(data_file, 'r')
    except IOError:
        raise("Can not find file{}".format(data_file))
    else:
        print("loading data from {} >>>".format(data_file))
        lines = F.readlines()
        F.close()
        print("data loaded from {} !".format(data_file))

    lines = filter_lines(lines, labels_list)
    tmp_labels = []
    tmp_data = []
    for line in lines:
        tmp = line.strip('\n').split(',')
        tmp_labels.append(tmp[0])
        if len(tmp[2]) >= max_len:
            tmp_start = tmp[2][:max_len]
            tmp_end = tmp[2][-max_len:]
        else:
            tmp_start = tmp[2]+"0"*(max_len - len(tmp[2]))
            tmp_end = "0"*(max_len - len(tmp[2]))+tmp[2]
            
        tmp_data.append(','.join((tmp[1],tmp_start))+','.join((tmp_end,tmp[3])))
    tmp_labels = label_trans(tmp_labels, labels_list)
    #### 如果使用softmax + cross_entrophy label_type为 "int64"
    tmp_labels = np.array(tmp_labels).astype('int64')
    data0 = []
    for i in tmp_data:
        tmp = []
        for j in to_one_hot_coding(i):
            tmp.append(j)
        data0.append(tmp)
    tmp_data = np.array(data0).astype("float32")
    print(tmp_data.shape)
    fine_tune_set = tmp_data
    fine_tune_labels = tmp_labels
    fine_tune_ds = tf.data.Dataset.from_tensor_slices((fine_tune_set, fine_tune_labels)).batch(100)
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
    print("{}/final_parameters/CnnAtt_{}_final".format(sys.path[0], model_wieght))
    model.load_weights(filepath="{}/final_parameters/CnnAtt_{}_final".format(sys.path[0], model_wieght))
    fine_tune_at = 4
    for layer in model.layers[:fine_tune_at]:
        layer.trainable =  False
    loss_object = tf.keras.losses.SparseCategoricalCrossentropy()
    optimizer = tf.keras.optimizers.Adam()
    train_loss = tf.keras.metrics.Mean(name='train_loss')
    train_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(name='train_accuracy')
    ### model train
    #@tf.function
    def train_step(images, labels):
        with tf.GradientTape() as tape:
            predictions = model(images)
            loss = loss_object(labels, predictions)
        gradients = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))
    
        train_loss(loss)
        train_accuracy(labels, predictions)

    now = time.time()
    for epoch in range(EPOCHS):
        # 在下一个epoch开始时，重置评估指标
        train_loss.reset_states()
        train_accuracy.reset_states()
        for images, labels in fine_tune_ds:
            train_step(images, labels)
    
    
        template = 'Epoch {}, Loss: {}, Accuracy: {}'
        print (template.format(epoch+1,
                               train_loss.result(),
                               train_accuracy.result()*100))
    now = time.time() - now
    print(now)
    model.save_weights(filepath="{}/final_parameters/CnnAtt_fine_tune_final".format(sys.path[0]), overwrite=True, save_format=None)
    print("parameter have saved!")
            
