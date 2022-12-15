import argparse
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, cross_val_predict, GridSearchCV

def read_arg():
    parser = argparse.ArgumentParser()
    parser.description = f'Training script with random forest model.'
    parser.add_argument('-i', '--input', required=True, help='Data frame for model training, csv format.')
    args = parser.parse_args()
    print(parser.description)
    return args

def training(file):
	# load data
	train = pd.read_csv(file)
	x_train = train[["s__[Eubacterium]_hallii",
		"s__Bacteroides_coprophilus",
		"s__Bacteroides_dorei",
		"s__Bacteroides_faecis",
		"s__Bacteroides_massiliensis",
		"s__Bacteroides_sp._CAG:98",
		"s__Bacteroides_stercoris",
		"s__Bacteroides_stercoris_CAG:120",
		"s__Bifidobacterium_longum",
		"s__Bilophila_wadsworthia",
		"s__Blautia_obeum",
		"s__Blautia_sp._CAG:237",
		"s__Blautia_sp._GD8",
		"s__Blautia_sp._KLE_1732",
		"s__Blautia_sp._Marseille-P2398",
		"s__butyrate-producing_bacterium_SS3/4",
		"s__Clostridium_sp._CAG:217",
		"s__Clostridium_sp._CAG:302",
		"s__Clostridium_sp._CAG:417",
		"s__Clostridium_sp._CAG:75",
		"s__Dorea_longicatena",
		"s__Dorea_sp._CAG:105",
		"s__Enterobacter_cloacae",
		"s__Enterobacter_sp._GN02315",
		"s__Eubacterium_hallii_CAG:12",
		"s__Eubacterium_sp._CAG:146",
		"s__Eubacterium_sp._CAG:156",
		"s__Firmicutes_bacterium_CAG:110",
		"s__Firmicutes_bacterium_CAG:227",
		"s__Firmicutes_bacterium_CAG:41",
		"s__Klebsiella_pneumoniae",
		"s__Parasutterella_excrementihominis",
		"s__Ruminococcus_gnavus_CAG:126",
		"s__Ruminococcus_sp._5_1_39BFAA",
		"s__Ruminococcus_sp._CAG:17",
		"s__Ruminococcus_sp._CAG:9",
		"s__Shigella_sonnei"]]
	y_train = train[["group"]]
	# training
	min_samples_leaf = range(3,20,1)
	min_samples_split = range(3,20,1)
	random_state= range(21,100,1)
	tune_list = []
	for a in min_samples_leaf:
	    for b in min_samples_split:
	        for c in random_state:
	            rfc = RandomForestClassifier(
	            	n_estimators = 500, 
	            	max_depth = 5, 
	            	min_samples_leaf = a,
	            	min_samples_split = b, 
	            	random_state = c
	            	)
	            rfc.fit(x_train, y_train.values.ravel())
	            scores = cross_val_score(rfc, x_train, y_train.values.ravel(), cv=5)
	            rfc_y_predict = rfc.predict(x_train)
	            y_predprob = rfc.predict_proba(x_train)[:, 1]
	            tune_list.append([a, b, c])
	            print("min_samples_split: %d" % a)
	            print("min_samples_leaf: %d" % b)
	            print("random_state: %d" % c)
	            print("accuracy score (CV): %f" % scores.mean())
	            print("accuracy score (All): %f" % rfc.score(x_train, y_train))
	print("Finish !")

def main():
	arg = read_arg()
	training(arg.input)

if __name__ == '__main__':
    main()