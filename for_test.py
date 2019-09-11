from sklearn import linear_model
import numpy as np
from sklearn.linear_model import Lasso


clf = linear_model.Lasso(alpha=0.1)
A = np.array([[0, 0], [1, 1], [2, 2]])
B = np.array([0, 1, 2])
clf.fit(A, B)

Lasso(alpha=0.1, copy_X=True, fit_intercept=True, max_iter=1000,
   normalize=False, positive=False, precompute=False, random_state=None,
   selection='cyclic', tol=0.0001, warm_start=False)

print(clf.coef_)
print(clf.intercept_)
