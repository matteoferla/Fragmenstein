import os
import numpy as np
import pandas as pd
import scipy
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LassoCV
from sklearn.linear_model import LinearRegression, Ridge, LogisticRegression, Lasso
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import make_scorer
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
from sklearn.utils import parallel_backend
from sklearn.metrics import roc_auc_score

from fragmenstein.scoring.scoring_config import FEATURES_DIR
from matplotlib import pyplot as plt

# https://gist.github.com/fabianp/2020955

RANDOM_SEED = 131
DO_PLOT = True


fname_prepared_dataset = os.path.join( FEATURES_DIR, "features.csv")

data= pd.read_csv(fname_prepared_dataset)
data.dropna(axis=0, inplace=True)
print(data.shape)

# matteo_scoring= data["comRMSD"]+ data["∆∆G"] / 5 + (data["N_unconstrained_atoms"] - data["N_constrained_atoms"] / 2) / 5
# print( scipy.stats.pearsonr(matteo_scoring,  data.loc[:, "label"]))
# plt.scatter( matteo_scoring,  data.loc[:, "label"]); plt.show()


x= data.iloc[:, 2:]
y= data.loc[:, "label"]

print( list(data.columns))
print( np.sum(y), len(y)- np.sum(y))

# if DO_PLOT:
#     for variable in x.columns:
#         # data_clean= data[data["∆∆G"]<30]
#         data.plot( x=variable, y="label", kind="scatter")
#         plt.show()

data["label"] = y
print(data.head())

if DO_PLOT:
    for variable in x.columns:
        data.boxplot( variable, by= "label", showfliers=False)
        plt.show()



# scoring_function= 'r2' #make_scorer(balanced_accuracy_score) # 'accuracy' # lambda estimator, x,y: estimator.predict_proba(x)[:,-1]
#
# cl= LassoCV(  cv=LeaveOneOut(), max_iter= int(1e3))
# cl.fit(x, y)
# print(cl.alpha_)
# y_pred = cross_val_predict(Lasso(alpha=cl.alpha_, max_iter=1000), x, y, cv=LeaveOneOut())
# print(scipy.stats.pearsonr(y, y_pred))


scoring_function= 'accuracy' #make_scorer(balanced_accuracy_score) # 'accuracy' # lambda estimator, x,y: estimator.predict_proba(x)[:,-1]

cl= LogisticRegressionCV(Cs= np.linspace(1e-7,1e-1, 10),  cv=LeaveOneOut(), scoring=  scoring_function, max_iter= int(1e3), class_weight="balanced")
cl.fit(x, y)

y_pred = cross_val_predict(LogisticRegression(C=cl.C_[0], class_weight="balanced", max_iter=1000), x, y, cv=LeaveOneOut(), method="predict_proba")[:,1]
print(roc_auc_score(y, y_pred))
print(confusion_matrix(y, (y_pred>0.5).astype(np.int32) ) )



clf = RandomForestClassifier(n_jobs=1) #Initialize with whatever parameters you want to

param_grid = {
                 'n_estimators': [5, 10, 15, 20, 100],
                 'max_depth': [2, 5, 7, 9],
                 'class_weight': ["balanced"]
             }

with parallel_backend( 'multiprocessing' ):
    grid_clf = GridSearchCV(clf, param_grid, cv=LeaveOneOut(), scoring= scoring_function, n_jobs= 4)
    grid_clf.fit(x, y)

print( grid_clf.best_params_, grid_clf.best_score_)

y_pred = cross_val_predict(RandomForestClassifier(** grid_clf.best_params_), x, y, cv=LeaveOneOut(), method="predict_proba")[:,1]
print(roc_auc_score(y, y_pred))
print(confusion_matrix(y, (y_pred>0.5).astype(np.int32) ) )