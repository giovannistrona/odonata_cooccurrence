import shap
import csv
from sklearn.ensemble import RandomForestClassifier,RandomForestRegressor
from sklearn import preprocessing
#from sklearn.inspection import permutation_importance
from sklearn.tree import export_graphviz# Export as dot file
from datetime import datetime
from collections import namedtuple
from numpy import mean
from numpy import arange,array,corrcoef,inf,linspace,ma,mean,ones,percentile,trapz,where,zeros
from random import shuffle,sample
from subprocess import call

cooc = [i for i in csv.reader(open('supported_co-occurrence.csv','r'))]
tr = [i for i in csv.reader(open('traits.csv','r'))]

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
		row.append(max([k[j] for k in rep]))
	for j in [24,25]:
		row.append(mean([float(k[j]) for k in rep]))
	row.append(1*(mean([float(k[26]) for k in rep])>=0.5))
	tr_mean.append(row)




tr_dict = dict([[i[0],i[1:]] for i in tr_mean])
###make the dataset for the machine learning
spp_cooc = cooc[0]
cooc = [i for i in csv.reader(open('mean_co-occurrence.csv','r'))]

#####niche overlap mahalanobis

from scipy.spatial.distance import mahalanobis
import numpy as np
data = array([i[3:] for i in tr_mean])
v = np.linalg.inv(np.cov(data.T))


####complete dataset for random forest regressor
out = open('cooc_vs_dist.csv','w')
out.write('sp1,sp2,cooc,dist\n')
for i in range(1,len(cooc)):
	for j in range(1,len(cooc)):
		if i>j:
			c = float(cooc[i][j])
			try:
				sp1 = tr_dict[spp_cooc[i]][2:]
				sp2 = tr_dict[spp_cooc[j]][2:]
				d = mahalanobis(sp1, sp2, v)
				out.write(','.join(map(str,[spp_cooc[i],spp_cooc[j],c,d]))+'\n')
			except:
				pass



out.close()


####complete dataset for random forest regressor
yx = []
for i in range(1,len(cooc)):
	for j in range(1,len(cooc)):
		if i>j:
			c = float(cooc[i][j])
			try:
				sp1 = tr_dict[spp_cooc[i]]
				sp2 = tr_dict[spp_cooc[j]]
				yx.append([c]+sp1[2:]+sp2[2:])
			except:
				pass

out = open('complete_dataset_rf.csv','w')
out.write('cooc,'+','.join(['sp1_'+i for i in tr_exp[0][3:]]+['sp2_'+i for i in tr_exp[0][3:]])+'\n')

for i in yx:
	out.write(','.join(map(str,i))+'\n')


out.close()



rf = RandomForestRegressor(bootstrap=True,
               max_depth=100, max_features='auto', max_leaf_nodes=None,
               min_impurity_decrease=0.0,
               min_samples_leaf=1, min_samples_split=2,
               min_weight_fraction_leaf=0.0, n_estimators=1000, n_jobs=4,
               oob_score=True, random_state=0, verbose=0, warm_start=False)



out = open('rf_both_res.csv','w')
out_all = open('rf_both_res_xy.csv','w')
out.write('r2\n')
for rep in range(100):
	shuffle(yx)
	split = int(round(len(yx)*0.7))
	train = yx[:split]
	test = yx[split:]
	train_y = [i[0] for i in train]
	train_x = [i[1:] for i in train]
	test_y = [i[0] for i in test]
	test_x = [i[1:] for i in test]
	rf.fit(train_x,train_y)
	pred = rf.predict(test_x)
	r2 = corrcoef(pred,test_y)[0][1]**2
	for i in range(len(pred)):
		out_all.write(','.join(map(str,[pred[i],test_y[i]]))+'\n')
	out.write(str(r2)+'\n')
	print ('both',rep,r2)


out.close()
out_all.close()



###oob score on complete dataset
train_y = [i[0] for i in yx]
train_x = [i[1:] for i in yx]
rf.fit(train_x,train_y)
var_imp_comp = rf.feature_importances_
out_oob = open('oob_scores.csv','w')
out_oob.write('complete,'+str(rf.oob_score_)+'\n')


####barplots
explainer = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(array(train_x[0]))
x_var_names = ['sp1_'+i for i in tr_exp[0][3:]]+['sp2_'+i for i in tr_exp[0][3:]]

shap.bar_plot(explainer.shap_values(array(train_x[100])),
              feature_names=x_var_names,
              max_display=len(x_var_names))



##negative co-occurrence
cooc = [i for i in csv.reader(open('supported_co-occurrence.csv','r'))]
yx = []
for i in range(1,len(cooc)):
	for j in range(1,len(cooc)):
		if i>j:
			c = float(cooc[i][j])
			c1 = 0
			if c<-0.5:# or c>0.7:
				c1 = 1
			try:
				sp1 = tr_dict[spp_cooc[i]]
				sp2 = tr_dict[spp_cooc[j]]
				yx.append([c1]+sp1[2:]+sp2[2:])
			except:
				pass


#random forest
###generate balanced presence-absence dataset

rf = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
               max_depth=100, max_features='auto', max_leaf_nodes=None,
               min_impurity_decrease=0.0,
               min_samples_leaf=1, min_samples_split=2,
               min_weight_fraction_leaf=0.0, n_estimators=1000, n_jobs=4,
               oob_score=True, random_state=0, verbose=0, warm_start=False)


yx_pres = [i for i in yx if i[0]==1]
yx_abs = [i for i in yx if i[0]==0]

out = open('rf_neg_res.csv','w')
out.write('r2,sens,spec,eI,eII,tss\n')
for rep in range(100):
	yx = yx_pres+sample(yx_abs,len(yx_pres))
	shuffle(yx)
	split = int(round(len(yx)*0.7))
	train = yx[:split]
	test = yx[split:]
	train_y = [i[0] for i in train]
	train_x = [i[1:] for i in train]
	test_y_cont = [i[0] for i in test]
	test_y_bin = [i[0] for i in test]
	test_x = [i[1:] for i in test]
	rf.fit(train_x,train_y)
	pred = rf.predict_proba(test_x)
	r2 = corrcoef([i[1] for i in pred],test_y_bin)[0][1]**2
	pred_ = []
	for i in pred:
		pred_.append(1*(i[0]<=0.5))
	test_y = test_y_bin[:]
	a,b,c,d = 0.0,0.0,0.0,0.0
	for i in range(len(test_y)):
		if test_y[i]>0.5 and pred_[i]>0.5: #predicted & observed
			a+=1.0
		elif test_y[i]==0 and pred_[i]>0: #predicted not observed
			b+=1.0
		elif test_y[i]>0 and pred_[i]==0: #observed not predicted
			c+=1.0
		elif test_y[i]==0 and pred_[i]==0: #not observed not predicted
			d+=1.0
	sens = a/float(len([i for i in test_y if i>0])) #sensitivity
	spec = d/float(len([i for i in test_y if i==0])) #specificity
	eI = c/float(len([i for i in test_y if i>0])) # type II error rate
	eII = b/float(len([i for i in test_y if i==0])) # type I error rate
	tss = ((a*d)-(b*c))/((a+c)*(b+d))
	out.write(','.join(map(str,[r2,sens,spec,eI,eII,tss]))+'\n')
	print ('neg',rep,r2,tss)


out.close()


###oob score on complete dataset
train_y = [i[0] for i in yx]
train_x = [i[1:] for i in yx]
rf.fit(train_x,train_y)
var_imp_neg = rf.feature_importances_

out_oob.write('negative,'+str(rf.oob_score_)+'\n')



################Positive cooc
yx = []
for i in range(1,len(cooc)):
	for j in range(1,len(cooc)):
		if i>j:
			c = float(cooc[i][j])
			c1 = 0
			if c>0.5:
				c1 = 1
			try:
				sp1 = tr_dict[spp_cooc[i]]
				sp2 = tr_dict[spp_cooc[j]]
				yx.append([c1]+sp1[2:]+sp2[2:])
			except:
				pass


#random forest
###generate balanced presence-absence dataset
yx_pres = [i for i in yx if i[0]==1]
yx_abs = [i for i in yx if i[0]==0]

out = open('rf_pos_res.csv','w')
out.write('r2,sens,spec,eI,eII,tss\n')

for rep in range(100):
	yx = yx_pres+sample(yx_abs,len(yx_pres))
	shuffle(yx)
	split = int(round(len(yx)*0.7))
	train = yx[:split]
	test = yx[split:]
	train_y = [i[0] for i in train]
	train_x = [i[1:] for i in train]
	test_y_cont = [i[0] for i in test]
	test_y_bin = [i[0] for i in test]
	test_x = [i[1:] for i in test]
	rf.fit(train_x,train_y)
	pred = rf.predict_proba(test_x)
	r2 = corrcoef([i[1] for i in pred],test_y_bin)[0][1]**2
	pred_ = []
	for i in pred:
		pred_.append(1*(i[0]<=0.5))
	test_y = test_y_bin[:]
	a,b,c,d = 0.0,0.0,0.0,0.0
	for i in range(len(test_y)):
		if test_y[i]>0.5 and pred_[i]>0.5: #predicted & observed
			a+=1.0
		elif test_y[i]==0 and pred_[i]>0: #predicted not observed
			b+=1.0
		elif test_y[i]>0 and pred_[i]==0: #observed not predicted
			c+=1.0
		elif test_y[i]==0 and pred_[i]==0: #not observed not predicted
			d+=1.0
	sens = a/float(len([i for i in test_y if i>0])) #sensitivity
	spec = d/float(len([i for i in test_y if i==0])) #specificity
	eI = c/float(len([i for i in test_y if i>0])) # type II error rate
	eII = b/float(len([i for i in test_y if i==0])) # type I error rate
	tss = ((a*d)-(b*c))/((a+c)*(b+d))
	out.write(','.join(map(str,[r2,sens,spec,eI,eII,tss]))+'\n')
	print ('pos',rep,r2,tss)


out.close()



###oob score on complete dataset
train_y = [i[0] for i in yx]
train_x = [i[1:] for i in yx]
rf.fit(train_x,train_y)
var_imp_pos = rf.feature_importances_

out_oob.write('positive,'+str(rf.oob_score_)+'\n')

out_oob.close()


##################Variable importance (mean impurity decrease)

x_var_names = ['sp1_'+i for i in tr_exp[0][3:]]+['sp2_'+i for i in tr_exp[0][3:]]

out = open('var_imp.csv','w')
out.write('var_name,both,pos,neg\n')
for i in range(len(x_var_names)):
	out.write(','.join(map(str,[x_var_names[i],var_imp_comp[i],var_imp_pos[i],var_imp_neg[i]]))+'\n')


out.close()



from matplotlib import pyplot as plt
names = ['all_cooc','positive_cooc','negative_cooc']
var_imps = [var_imp_comp,var_imp_pos,var_imp_neg]
for i in range(3):
	var_imp = var_imps[i]
	plt.rcParams.update({'figure.figsize': (8.0, 6.0)})
	plt.rcParams.update({'font.size': 14})
	sorted_idx = var_imp.argsort()
	plt.barh(array(x_var_names)[sorted_idx][-10:], var_imp[sorted_idx][-10:])
	plt.xlabel("Random Forest Feature Importance")
	plt.savefig(names[i]+'.pdf')
	plt.clf()




