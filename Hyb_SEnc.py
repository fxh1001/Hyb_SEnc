# coding=UTF-8
from flask import Flask,render_template,request
import feature_extraction
import feature_selection_MD
import read_fasta_sequences
import numpy as np
import pandas as pd
import joblib
import feature_selection
from sklearn.preprocessing import MinMaxScaler
import sys,re
from check_sequences import *
from gevent import pywsgi

# initialize flask application
app = Flask(__name__)

@app.route('/')
def main():
   return render_template('services2.html')


@app.route('/submit',methods=["POST","GET"])
def predict():
    if request.method == "POST":
        sequencedata = request.form.get("sequence")
        select_classifier = request.form.get("classifier")
        fastas, sequence_name = read_fasta_sequences.read_protein_sequences(sequencedata)
        check=check_sequences(sequencedata, fastas)
        if check == False:
            return render_template("services2.html")

        else:
            encodings = feature_extraction.get_features(fastas)
            if select_classifier == 'Hyb_SEnc_RD':
                model = joblib.load('./Model/sclf_RD_special.pkl')
                newdataset = feature_selection.select_features(encodings)
            else:
                model = joblib.load('./Model/sclf_model.pkl')
                newdataset = feature_selection_MD.select_features(encodings)
            X = newdataset
            scaler = MinMaxScaler()
            scaler.fit(X)
            X = scaler.transform(X)
            y_pred_prob = model.predict_proba(X)
            df_out = pd.DataFrame(np.zeros((y_pred_prob.shape[0], 3)),
                                  columns=["Sequence_name", "Prediction", "probability"])
            y_pred = model.predict(X)
            for i in range(y_pred.shape[0]):
                df_out.iloc[i, 0] = str(sequence_name[i])
                if y_pred[i] == 1:
                    df_out.iloc[i, 1] = "AntiTb peptide"
                    df_out.iloc[i, 2] = "%.2f%%" % (y_pred_prob[i, 1] * 100)
                if y_pred[i] == 0:
                    df_out.iloc[i, 1] = "non-AntiTb peptide"
                    df_out.iloc[i, 2] = "%.2f%%" % (y_pred_prob[i, 0] * 100)

            results = df_out.to_numpy()
            return render_template('result.html', Fasta=sequencedata, result=results, Predicting="Output results ",Classifier = select_classifier)
    else:
        return render_template("services.html")

if __name__=="__main__":
    app.run(debug=True,host="119.13.83.151",port=5000)


