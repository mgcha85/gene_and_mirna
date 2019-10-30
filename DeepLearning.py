import numpy as np
import os
# from keras.models import Sequential, load_model
# from keras.layers import Dense, Dropout
# from keras.optimizers import SGD
# from keras.layers import Conv2D, MaxPooling2D, BatchNormalization
# from keras.layers import Flatten


class DeepLearning:
    def sequence_LSTM(self, x_train, y_train, x_test, y_test, fname):
        from keras.models import Sequential
        from keras.layers import LSTM, Dense

        input_shape = x_train.shape[1:]
        num_cls = y_train.shape[1]

        # expected input data shape: (batch_size, timesteps, data_dim)
        model = Sequential()
        model.add(LSTM(32, return_sequences=True, input_shape=input_shape))
        model.add(LSTM(32, return_sequences=True))
        model.add(LSTM(32))
        model.add(Dense(num_cls, activation='softmax'))

        model.compile(loss='categorical_crossentropy',
                      optimizer='rmsprop',
                      metrics=['accuracy'])

        model.fit(x_train, y_train, batch_size=32, epochs=300, validation_data=(x_test, y_test))
        model.save(fname)

    def vgg_like_convnet(self, x_train, y_train, x_test, y_test, fname, layer_sizes, filter_sizes):
        input_shape = x_train.shape[1:]
        num_cls = y_train.shape[1]

        model = Sequential()
        for i, fsize in enumerate(layer_sizes):
            if i == 0:
                model.add(Conv2D(fsize, filter_sizes[i], activation='relu', input_shape=input_shape))
            else:
                model.add(Conv2D(fsize, filter_sizes[i], activation='relu'))
            model.add(BatchNormalization())
            model.add(MaxPooling2D(pool_size=(1, 2)))
            model.add(Dropout(0.25))

        model.add(Flatten())
        model.add(Dense(layer_sizes[-1] * 2, activation='relu'))
        model.add(Dropout(0.5))
        model.add(Dense(num_cls, activation='softmax'))

        sgd = SGD(lr=1e-2, decay=1e-6, momentum=0.9, nesterov=True)
        model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])

        model.fit(x_train, y_train, batch_size=16, epochs=300)
        score = model.evaluate(x_test, y_test, batch_size=16)
        print(score)
        model.save(fname)

    def split_data_set(self, x_trn, y_trn, ratio=0.2, pred=0):
        N = len(x_trn)
        np.random.seed(0)
        ridx = np.random.randint(0, N, N)
        x_trn = x_trn[ridx]
        y_trn = y_trn[ridx]

        N = len(x_trn)
        M = int(N * (1 - ratio))
        if pred > 0:
            M2 = int(N * (1 - 2 * ratio))
            return x_trn[:M2], y_trn[:M2], x_trn[M2: M], y_trn[M2: M], x_trn[M:], y_trn[M:]
        else:
            M = int(N * (1 - ratio))
            return x_trn[:M], y_trn[:M], x_trn[M:], y_trn[M:]

    def prediction(self, x_pred, y_pred, fname):
        from keras.models import load_model
        model = load_model(fname)

        predictions = model.predict(x_pred).argmax(axis=1)
        gtruth = y_pred.argmax(axis=1)

        results = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0}
        for p, g in zip(predictions, gtruth):
            if p == 0 and g == 0:
                results['tp'] += 1
            elif p == 0 and g == 1:
                results['tn'] += 1
            elif p == 1 and g == 0:
                results['fp'] += 1
            else:
                results['fn'] += 1

        print(results)
