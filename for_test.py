# from __future__ import print_function
# import keras
# from keras.datasets import mnist
# from keras.models import Sequential
# from keras.layers import Dense, Dropout, Flatten
# from keras.layers import Conv2D, MaxPooling2D
# from keras import backend as K
#
#
# batch_size = 128
# num_classes = 10
# epochs = 12
#
# # input image dimensions
# img_rows, img_cols = 28, 28
#
# # the data, split between train and test sets
# (x_train, y_train), (x_test, y_test) = mnist.load_data()
#
# if K.image_data_format() == 'channels_first':
#     x_train = x_train.reshape(x_train.shape[0], 1, img_rows, img_cols)
#     x_test = x_test.reshape(x_test.shape[0], 1, img_rows, img_cols)
#     input_shape = (1, img_rows, img_cols)
# else:
#     x_train = x_train.reshape(x_train.shape[0], img_rows, img_cols, 1)
#     x_test = x_test.reshape(x_test.shape[0], img_rows, img_cols, 1)
#     input_shape = (img_rows, img_cols, 1)
#
# x_train = x_train.astype('float32')
# x_test = x_test.astype('float32')
# x_train /= 255
# x_test /= 255
# print('x_train shape:', x_train.shape)
# print(x_train.shape[0], 'train samples')
# print(x_test.shape[0], 'test samples')
#
# # convert class vectors to binary class matrices
# y_train = keras.utils.to_categorical(y_train, num_classes)
# y_test = keras.utils.to_categorical(y_test, num_classes)
#
# model = Sequential()
# model.add(Conv2D(32, kernel_size=(3, 3),
#                  activation='relu',
#                  input_shape=input_shape))
# model.add(Conv2D(64, (3, 3), activation='relu'))
# model.add(MaxPooling2D(pool_size=(2, 2)))
# model.add(Dropout(0.25))
# model.add(Flatten())
# model.add(Dense(128, activation='relu'))
# model.add(Dropout(0.5))
# model.add(Dense(num_classes, activation='softmax'))
#
# model.compile(loss=keras.losses.categorical_crossentropy,
#               optimizer=keras.optimizers.Adadelta(),
#               metrics=['accuracy'])
#
# model.fit(x_train, y_train,
#           batch_size=batch_size,
#           epochs=epochs,
#           verbose=1,
#           validation_data=(x_test, y_test))
# score = model.evaluate(x_test, y_test, verbose=0)
# print('Test loss:', score[0])
# print('Test accuracy:', score[1])


import sqlite3
from Database import Database
import pandas as pd

# df = pd.read_csv("/home/mingyu/Downloads/Galaxy37.txt", sep='\t')
# df.to_sql('Galaxy37', sqlite3.connect("/home/mingyu/Downloads/Galaxy37.db"), index=None)

con = sqlite3.connect("/home/mingyu/Downloads/Galaxy37.db")
df = pd.read_sql("SELECT * FROM 'Galaxy37_org'", con)
target_id = df['target_id'].str.split('|', expand=True)
target_id.columns = ['reference_id', 'ref_gene_id', 'dummy_id', 'dummy_gene_id', 'transcript_id', 'gene_id', 'length', 'label', 'dummy']
target_id = target_id[['reference_id', 'ref_gene_id', 'transcript_id', 'gene_id']]

df_res = pd.concat([target_id, df.drop('target_id', axis=1)], axis=1)
df_res.to_sql('Galaxy37', con, index=None, if_exists='replace')
