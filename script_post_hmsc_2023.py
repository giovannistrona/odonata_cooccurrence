import shap
import csv
from sklearn.ensemble import RandomForestClassifier,RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV #GridSearchCV
from numpy import array,corrcoef,mean
from random import shuffle,sample
import numpy as np
from scipy.spatial.distance import mahalanobis
import matplotlib.pyplot as plt

cooc = [i for i in csv.reader(open('supported_co-occurrence.csv','r'))]
tr = [i for i in csv.reader(open('traits.csv','r'))]
#tr = [tr[0]]+[i for i in tr if i[0] in cooc[0]]
nas = [sum([i[j]=='NA' for i in tr]) for j in range(23)]
for i in range(len(tr[0])):
	print (tr[0][i],nas[i])



###limit dataset
#(0, 'GenusSpecies')
#(5, 'Flight_season_start')
#(6, 'Flight_season_end')
#(7, 'body_colors')
#(8, 'body_colortypes')
#(9, 'body_patterns')
#(11, 'flight_mode')
#(13, 'aquatic_habitats')
#(14, 'habitat_openness')
#(15, 'body_lengths')
#(16, 'hindwing_lengths')
#(17, 'has_wing_pigment')


tr_ok = [[i[j] for j in [0,5,6,7,8,9,11,13,14,15,16,17]] for i in tr]
tr_ok = [i for i in tr_ok if 'NA' not in i]
for i in range(1,len(tr_ok)):
	if tr_ok[i][11]=='yes':
		tr_ok[i][11] = 1
	else:
		tr_ok[i][11] = 0



all_body_col = sorted(list(set([i[3] for i in tr_ok[1:]])))
all_body_col_t = sorted(list(set([i[4] for i in tr_ok[1:]])))
all_body_col_p = sorted(list(set([i[5] for i in tr_ok[1:]])))
all_fl_mods = sorted(list(set([i[6] for i in tr_ok[1:]])))
all_aq_hab = sorted(list(set([i[7] for i in tr_ok[1:]])))
all_hab_op = sorted(list(set([i[8] for i in tr_ok[1:]])))

to_expand_ids = [3,4,5,6,7,8]
to_expand_data = [all_body_col,
		  all_body_col_t,
		  all_body_col_p,
		  all_fl_mods,
		  all_aq_hab,
		  all_hab_op]


to_expand_lab = ['body_col_',
		'body_col_t_',
		'body_col_p_',
		'fl_mods_',
		'aq_hab_',
		'all_hab_op_']



tr_exp = []
row = []
for j in range(len(tr_ok[0])):
	if j in to_expand_ids:
		k = to_expand_ids.index(j)
		row += [to_expand_lab[k]+z for z in to_expand_data[k]]
	else:
		row.append(tr_ok[0][j])


tr_exp.append(row)
for i in tr_ok[1:]:
	row = []
	for j in range(len(i)):
		if j in to_expand_ids:
			k = to_expand_ids.index(j)
			row += [1 if i[j]==z else 0 for z in to_expand_data[k]]
		else:
			row.append(i[j])
	tr_exp.append(row)


spp = sorted(list(set([i[0] for i in tr_exp])&set([j[0] for j in cooc])))
tr_mean = []
for i in spp:
	rep = []
	for j in tr_exp:
		if j[0]==i:
			rep.append(j)
	row = [rep[0][j] for j in range(3)]
	for j in range(3,24):
		row.append(max([float(k[j]) for k in rep]))
	for j in [24,25]:
		row.append(mean([float(k[j]) for k in rep]))
	row.append(1*(mean([float(k[26]) for k in rep])>=0.5))
	tr_mean.append(row)


tr_dict = dict([[i[0],i[1:]] for i in tr_mean])
###make the dataset for the machine learning

cooc = [i for i in csv.reader(open('mean_co-occurrence.csv','r'))]
spp_cooc = cooc[0]


from random import random
#####niche overlap mahalanobis
data = array([i[3:] for i in tr_mean if i[0] in spp])
for i in range(len(data)):
	for j in range(len(data[0])):
		data[i][j]+=random()/10**10


v = np.linalg.inv(np.cov(data.T))
out = open('cooc_vs_dist.csv','w')
out.write('sp1,sp2,cooc,dist,'+','.join(['tr_'+str(i) for i in range(24)])+'\n')
for i in range(1,len(cooc)):
	for j in range(1,len(cooc)):
		if i>j:
			c = float(cooc[i][j])
			try:
				sp1 = tr_dict[spp_cooc[i-1]][2:]
				sp2 = tr_dict[spp_cooc[j-1]][2:]
				d = mahalanobis(sp1, sp2, v)
				ind_d = []
				for trait in range(len(sp1)):
					tr1,tr2 = sp1[trait],sp2[trait]
					ind_d.append(abs(tr1-tr2))
				out.write(','.join(map(str,[spp_cooc[i],spp_cooc[j],c,d]+ind_d))+'\n')
			except:
				pass


out.close()

from datetime import datetime
from collections import namedtuple
###RANDOM FOREST MODELS
####complete dataset for random forest regressor
cooc = [i for i in csv.reader(open('mean_co-occurrence.csv','r'))]
yx = []
for i in range(1,len(cooc)):
	for j in range(1,len(cooc)):
		if i>j:
			c = float(cooc[i][j])
			try:
				sp1 = tr_dict[cooc[i][0]]
				sp2 = tr_dict[cooc[j][0]]
				#ov = overlap(sp1[0],sp1[1],sp2[0],sp2[1])
				yx.append([c]+sp1[2:]+sp2[2:])
			except:
				pass


out = open('complete_dataset_rf.csv','w')
out.write('cooc,'+','.join(['sp1_'+i for i in tr_exp[0][3:]]+['sp2_'+i for i in tr_exp[0][3:]])+'\n')
for i in yx:
	out.write(','.join(map(str,i))+'\n')


out.close()
