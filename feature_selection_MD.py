import numpy as np
import pandas as pd

feature_index = pd.read_csv("File_list/feature_index_select.csv")

def select_features(all_features):
    print('Feature selection...')
    new_feature = []
    original_data = pd.DataFrame(all_features)
    for i in list(feature_index.iloc[0, 1:]):
        new_feature.append(original_data[int(i)])  # 选择特征
    features = np.array(new_feature).T  # 在转化为矩阵，feature即为最终选出的最优特征子集
    return features