import json,sys
import os.path
import numpy as np
import itertools
import yaml
import shutil


def give_sub_name(dic, namekey):
	name = ''
	for key in namekey:
		name += '%s(%s)'%(key, dic[key])
	return name
def generate_sub_param(name, dic, path):
	with open(path + '/%s.json'%name, 'w') as file:
		json.dump(dic, file, indent=4, sort_keys=True, separators=(',', ': '))
def default_param(dic, key, value):
	if not key in dic:
		dic[key] = value
def train_sub_param(dic):	
	#Shell section:
	default_param(dic, 'integrator', 'langevin')
	default_param(dic, 'dt', 0.01)
	default_param(dic, 'mtheta0', 0.)
	# default_param(dic,'theta0',    2.*np.arcsin(dic['l0']/np.sqrt(12*dic['R0']**2 - 3*dic['l0']**2)))
	return dic

prompt = '> '
path = os.path.abspath(os.getcwd())
paramfile = 'parent_param.yaml'
if len(sys.argv) > 1: paramfile = sys.argv[1]
#read parent parameter file
with open(path + '/' + paramfile, 'r') as f:
	    doc = yaml.load(f, Loader = yaml.FullLoader)
value = []
subkey = []
overwriteflag = False
for key in doc.keys():
	if (key != 'Name') & (key != 'Path'):
		subkey += doc[key].keys()
		value += doc[key].values()

for i, sub_value in enumerate(itertools.product(*value)):
	#generate sub dictionary:
	sub_dict = dict(zip(subkey, sub_value))
	sub_name = give_sub_name(sub_dict, doc['Name']['key'])
	print(sub_name)
	sub_path = os.path.join(path, sub_name)
	print(path, sub_path)
	if not os.path.exists(sub_path):
		os.makedirs(sub_path)
	elif overwriteflag:
		shutil.rmtree(sub_path)
		os.makedirs(sub_path)
	else: 
		print('Sub Directory Exist\nDo you want to cover?')
		cover = input(prompt)
		if cover == 'yes': shutil.rmtree(sub_path); os.makedirs(sub_path)
		elif cover == 'allyes': shutil.rmtree(sub_path); os.makedirs(sub_path); overwriteflag = True
		elif cover == 'no': print('overwrite dir reject')
		elif cover == 'allno': print('overwrite reject for all'); break
		else: print('unknown command')
	sub_dict = train_sub_param(sub_dict)
	generate_sub_param(sub_name, sub_dict, sub_path)
