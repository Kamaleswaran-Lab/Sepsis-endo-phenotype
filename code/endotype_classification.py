import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import normalize
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE

meta_train = pd.read_csv("./trainlabel.csv", index_col=0) 
data_train = pd.read_csv("./trainmatrix.csv", index_col=0) #normalized matrix
data_train = data_train.T

data_pred = pd.read_csv("/psudobulk_deseq_norm.csv", index_col=0) # scRNA-seq to psudobulk matrix, then use DEseq to normalize the data
data_pred = data_pred.T

X_train, X_test, y_train, y_test = train_test_split(data_train1, meta_train["label"], test_size=0.3, random_state=400)
svm_model = SVC(kernel='linear', random_state=42)
svm_model.fit(X_train, y_train)
y_pred = svm_model.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)

all_feature = data_train.columns
data_pred = data_pred[all_feature]

y_pred = svm_model.predict(data_pred)
df = pd.DataFrame(y_pred, columns=['predict'])
df.index = data_pred.index

df.to_csv('./label_assignment.csv', index=True)




