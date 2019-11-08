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

con = sqlite3.connect('/home/mingyu/Bioinformatics/database/Fantom/v5/cluster/result_all_clusters.db')
con_out = sqlite3.connect('/home/mingyu/Bioinformatics/database/Fantom/v5/cluster/result_all_clusters_spt.db')
df = pd.read_sql("SELECT * FROM 'all_clusters'", con)
for chr, df_chr in df.groupby('chromosome'):
    for str, df_str in df_chr.groupby('strand'):
        df_str.drop(['chromosome', 'strand'], axis=1).to_sql('_'.join([chr, str]), con_out, index=None)

# loc = df['location'].str[1:-1].str.split(',', expand=True)
# loc.columns = ['start', 'end']
# df = pd.concat([df.drop('location', axis=1), loc], axis=1)
# df = df[['chromosome', 'start', 'end', 'strand', 'gene_name', 'transcript_id', 'attribute']]
# df.to_sql('all_clusters', con, if_exists='replace', index=None)
