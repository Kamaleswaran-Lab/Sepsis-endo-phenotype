import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import normalize
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE

meta_train = pd.read_csv("./trainlabel.csv", index_col=0) 
data_train = pd.read_csv("./trainmatrix.csv", index_col=0) #normalized matrix from original endo/phenotype paper
data_train = data_train.T

data_pred = pd.read_csv("/psudobulk_deseq_norm.csv", index_col=0) # scRNA-seq to psudobulk matrix, then use DEseq to normalize the data
data_pred = data_pred.T


##fine-tune SVM
X_train, X_test, y_train, y_test = train_test_split(data_train1, meta_train["label"], test_size=0.3, random_state=400)

param_grid = {
    'C': [0.01, 0.1, 1, 10],
    'kernel': ['linear', 'rbf'],
    'gamma': ['scale', 'auto']
}

svc = SVC()
grid = GridSearchCV(svc, param_grid, cv=5, scoring='accuracy')
grid.fit(X_train, y_train)

print("Best parameters:", grid.best_params_)
print("Best cross-val accuracy:", grid.best_score_)

best_model = grid.best_estimator_
y_val_pred = best_model.predict(X_test)

print("Validation Classification Report:")
print(classification_report(y_test, y_val_pred))


#predict on our own data
all_feature = data_train.columns
data_pred = data_pred[all_feature]

y_pred = best_model.predict(data_pred)
df = pd.DataFrame(y_pred, columns=['predict'])
df.index = data_pred.index

df.to_csv('./label_assignment.csv', index=True)




