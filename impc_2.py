import csv
import pandas
import numpy
import sys
import requests
from bs4 import BeautifulSoup
import os
import subprocess
import time
import datetime
import re
import glob
from pprint import pprint

directory = 'C:/Users/svprice/Dropbox/IMPC_Data/'

# get list of urls to protocol parameter pages
def get_protocol_list(pipeline = "IMPC"):
	if pipeline == "IMPC":
		url = "https://www.mousephenotype.org/impress/procedures/7"
	if pipeline == "JAX":
		url = "https://www.mousephenotype.org/impress/procedures/12"

	# Scrape the HTML at the url
	r = requests.get(url)

	# Turn the HTML into a Beautiful Soup object
	soup = BeautifulSoup(r.text, "lxml")

	f = open(directory + 'protocol_list.txt','w+')

	# create a list of links to the protocol pages
	protocol_list = []
	for a in soup.find_all(id=="Display Protocol"):
		if a.has_attr('href'):
			if "protocol" in a['href'] and a['href'] not in protocol_list and "IMPC" in a.text:
				name = a.text.split(" IMPC")[0]
				code = "IMPC" + a.text.split(" IMPC")[1]
				link = a['href']
				protocol_list.append([name, code, link])
				f.write(name + ' ' + code + '\n')
			if "protocol" in a['href'] and a['href'] not in protocol_list and "JAX" in a.text:
				name = a.text.split(" JAX")[0]
				code = "IMPC" + a.text.split(" JAX")[1]
				link = a['href']
				protocol_list.append([name, code, link])
				f.write(name + ' ' + code + '\n')
	return protocol_list

# get procedure group data
def return_procedure_group(procedure_code, number_of_rows, rewrite = False):
		# project directory name
		project_directory = directory + procedure_code

		# create directory if it doesn't exist
		try:
			os.stat(project_directory)
		except:
			os.mkdir(project_directory)

		# output directory name
		output_directory = directory + procedure_code + '/RawData'

		# create directory if it doesn't exist
		try:
			os.stat(output_directory)
		except:
			os.mkdir(output_directory)
		
		# output file name
		output_file = output_directory + '/' + procedure_code + '_data.csv'

		# download raw data
		if rewrite or not os.path.isfile(output_file):
			file = open(output_file,'w+')
			print(procedure_code + ": BEGINNING DOWNLOAD OF RAW DATA")
			curl = subprocess.Popen(["curl", 'http://wwwdev.ebi.ac.uk/mi/impc/dev/solr/experiment/select?q=metadata%3A*&fq=procedure_group:*_' + procedure_code +'&wt=csv&rows=' + str(number_of_rows)], stdout = file)
			curl.wait()
			print(procedure_code + ": DOWNLOAD COMPLETED")
		else:
			print(procedure_code + ": DATA ALREADY DOWNLOADED")

# get data and split by center (AAA input)
# no rewrite parameter - will always run
def get_and_split_by_center(procedure_code, number_of_rows,fileName = ""):
	
	# project directory name
	project_directory = directory + procedure_code

	# create directory if it doesn't exist
	try:
		os.stat(project_directory)
	except:
		os.mkdir(project_directory)


	intermediate_directory = directory + procedure_code + '/IntermediaryFiles'
	# create directory if it doesn't exist
	try:
		os.stat(intermediate_directory)
	except:
		os.mkdir(intermediate_directory)
	

	# output directory name
	output_directory = intermediate_directory + '/ByCenter/'
	
	# create directory if it doesn't exist
	try:
		os.stat(output_directory)
	except:
		os.mkdir(output_directory)
	
	if(fileName == ""):
		fileName = project_directory + '/RawData/' + procedure_code + '_data.csv'

	file_read = csv.reader(open(fileName,'r'))
	
	row_count = sum(1 for row in file_read)
	if row_count > 1:
		# awk code for splitting
		args = ["gawk", "-F','", "NR==1{header=$0}NR>1&&!a[$1]++{print header > (\"" + output_directory + "\"$1\".csv\")}NR>1{print > (\"" + output_directory + "\"$1\".csv\")}", fileName]
		subprocess.Popen(args, shell=True).wait()
		return True
	else:
		return False

# creates a new column to uniquely identify each time series data point
# creates phenotype_production column too
def process_col(procedure_code, phenotyping_center = ''):
	intermediate_directory = directory + procedure_code + '/IntermediaryFiles'
	input_directory = intermediate_directory + '/ByCenter'
	output_directory = intermediate_directory + '/ByCenterProcessed'
	try:
		os.stat(output_directory)
	except:
		os.mkdir(output_directory)

	path = input_directory + '/*.csv'
	if(phenotyping_center == ''):
		for fn in glob.glob(path):
			df = pandas.read_csv(fn, low_memory = False) 
			df["parameter_name"] = df["parameter_name"] + df["discrete_point"].astype(str).str.replace('nan','')
			if "production_phenotype" not in df.columns.values:
				df["production_phenotype"] = df["production_center_id"].astype(str) + "_" + df["phenotyping_center_id"].astype(str)
			print(fn)
			df.to_csv(output_directory + '/' + fn.split('\\')[1].split('.')[0] + "_processed.csv", sep=',', index = False)
			print("COMPLETED PROCESSING COLUMNS: " + fn)
	else:
		fn = input_directory + '/' + phenotyping_center + '.csv'
		df = pandas.read_csv(fn, low_memory = False) 
		df["parameter_name"] = df["parameter_name"] + df["discrete_point"].astype(str).str.replace('nan','')
		if "production_phenotype" not in df.columns.values:
			df["production_phenotype"] = df["production_center_id"].astype(str) + "_" + df["phenotyping_center_id"].astype(str)
		df.to_csv(output_directory + '/' + phenotyping_center + "_processed.csv", sep=',', index = False)
		print("COMPLETED PROCESSING COLUMNS: " + fn)

# pivots data (AAA input)
def get_and_split_by_animal(procedure_code, rewrite = False):
	
	print("STARTED PIVOT")

	intermediate_directory = directory + procedure_code + '/IntermediaryFiles'

	# for printing an error message
	hasMultipleMetadataPerMouse = False

	# project directory name
	project_directory = intermediate_directory + '/ByCenterProcessed'


	# create directory if it doesn't exist
	try:
		os.stat(project_directory)
	except:
		os.mkdir(project_directory)
	# output directory name
	output_directory = intermediate_directory + '/ByCenterProcessedPivot'
	# create directory if it doesn't exist
	try:
		os.stat(output_directory)
	except:
		os.mkdir(output_directory)
	data_empty = False
	for fn in os.listdir(project_directory):
		outputFile = output_directory + '/' + fn.split('.')[0] + '_pivot.csv'
		if not os.path.isfile(outputFile) or rewrite:
			fileName = project_directory + '/' + fn
			df = pandas.read_csv(fileName, low_memory = False, encoding = "ISO-8859-1")
			if(len(df['data_point'].unique()) == 1 and df['data_point'].unique()[0]):
				data_empty = True
			if(data_empty):
				print("ERROR: " + fn.split('.')[0] + ": NO DATA UNABLE TO PIVOT")
			else:
				df.loc[:,'date_of_birth'] = df.loc[:,'date_of_birth'].apply(time_split)
				df.loc[:,'date_of_experiment'] = df.loc[:,'date_of_experiment'].apply(time_split)
				if len(df['weight_date'].unique()) > 1:
					df.loc[:,'weight_date'] = df.loc[:,'weight_date'].apply(time_split)
				# list of all the column names
				col_names = df.columns.values
				# create index for future pivot table
				index_col = ['biological_sample_id']
				for element in col_names:
					if str(element) not in ['data_point', 'parameter_name', 'biological_sample_id', 'Unnamed: 0']:
						index_col.append(element)
				# list of all the unique values in biological_sample_id column
				col_values = df['biological_sample_id'].unique()
				# list of columns without a single value per mouse
				non_pivot_col = []
				# list of non-empty columns
				non_empty_col = []
				subsets = []
				# iterates by mouse (biological_sample_id value)
				index = 1
				for element in col_values:
					# gets specific mouse data frame
					df_subset = df.loc[df['biological_sample_id'] == element]
					df_subset.is_copy = False
					df_subset = df_subset.reset_index(drop = True)
					for column in df_subset:
						# adds columns to non_pivot_col
						if len(df_subset[column].unique())>1:
							if column != 'parameter_name' and column != 'data_point':
								# gets the smallest metadata value - needed for ACS
								if column == 'metadata':
									min_str = df_subset.loc[0,column]
									for i in range(1, len(df_subset.loc[:,column])):
										val = df_subset.loc[i,column]
										if(len(val) < len(min_str)):									
											if(not hasMultipleMetadataPerMouse):
												hasMultipleMetadataPerMouse = True
												print("WARNING: " + procedure_code + ":	" + fn.split('.')[0] + " " + "has multiple metadata values per mouse (in order to pivot will only take the smallest common sample of metadata per mouse)")
											min_str = val
									df_subset.loc[:, column] =	val							
								else:
									df_subset.loc[:, column] = 'Error'
						if str(df_subset[column].unique()[0]) == 'nan':
							df_subset.loc[:, column] = 'NA'
					subsets.append(df_subset)
				df_all = pandas.concat(subsets)

				df_all = df_all.pivot_table(values='data_point', index = index_col, columns = 'parameter_name', aggfunc = lambda x: min(x))

				df_all.to_csv(output_directory + '/' + fn.split('.')[0] + '_pivot.csv', sep = ',')
				print("COMPLETED PIVOT: " + fn.split('.')[0])
				#except:
				#	print("ERROR: " + procedure_code + ": " + fn.split('.')[0] + "	" + "unable to pivot (see impc.py get_and_split_by_animal)")
				#	print(sys.exc_info[0])
		else:
			print("WIDE DATA: " + fn + " already exists")
	if(data_empty):
		print("BECAUSE OF MISSING DATA PROGRAM TERMINATED EARLY")
		return False
	else:
		print("FINISHED PIVOT")
		return True

# used in split_metadata_robust
def time_split(time):
	if isinstance(time, str):
		return time.split('T')[0]

# used in split_metadata_robust
def time_split2(time):
	if isinstance(time, str):
		try:
			output = time.split('T')[1].split('+')[0]
		except:
			output = ''
		if "-" in output:
			output = output.split("-")[0]
		return output

# split up metadata (AAA input)
def split_metadata_robust(procedure_code, rewrite = False):
	
	intermediate_directory = directory + procedure_code + '/IntermediaryFiles'
	# input directory name
	input_directory = intermediate_directory + '/ByCenterProcessedPivot'
	# output directory name
	output_directory = directory + procedure_code + '/WideData'
	# create directory if it doesn't exist
	try:
		os.stat(output_directory)
	except:
		os.mkdir(output_directory)
	# cycle through files in the input directory
	for fn in os.listdir(input_directory):	
		phenotype_center = fn.split('.')[0]
		outputFile = output_directory + '/' + phenotype_center + '_final.csv'
		if not os.path.isfile(outputFile) or rewrite:
			print("STARTED SPLITTING METADATA: " + fn.split('.')[0])
			#get phenotype center from file name
			index = 0
			# create dataframe from csv table
			df = pandas.read_csv(input_directory + '/' + phenotype_center + '.csv', low_memory = False, encoding = 'ISO-8859-1')
			# create matrix to hold new column names and values
			# matrix in form [row num][column num]
			# set width
			# w = numpy.array([len(l.split(',')) for l in df['metadata']]).max()
			# finds the length of the largest number of metadata values to create a numpy array with width of that size
			w = numpy.array([len(re.split("(?<!\\\)\,",l)) for l in df['metadata']]).max()
			# set height
			h = len(df['metadata'])+1
			Matrix = numpy.array([['' for x in range(w)] for y in range(h)],dtype=object)
			# takes only the first of several metadata values
			if re.search("_(?![0-9]{3}|[A-Z]{3})", str(df['metadata'][0])) is not None:
				print("WARNING: " + procedure_code + ": " + phenotype_center + " " + "mouse has metadata from multiple tests (in order to expand will only take the first)")

			pre_split = re.split("_(?![0-9]{3}|[A-Z]{3})", str(df['metadata'][0]))[0]
			# pre_split = str(df['metadata'][0]).split('_')[0]
			# input column names into matrix
			# first_split = pre_split.split(',')
			first_split = re.split("(?<!\\\)\,",pre_split)
			for index in range(len(first_split)):
				# Matrix[0][index] = str(first_split[index].split(" = ")[0])
				Matrix[0][index] = str(re.split("\s=\s|(?<=D):\s",first_split[index])[0])
			# fill matrix with metadata
			for i in range(len(df['metadata'])):
				splited = re.split('(?<!\\\),',str(df['metadata'][i]))
				for j in range(len(splited)):
					if "\\," not in splited[j]:	
						try:
							# val = str(splited[j].split(" = ")[1])
							val = re.split("\s=\s|(?<=D):\s",splited[j])[1]
							# if str(splited[j].split(" = ")[0]) == Matrix[0][j]:
							argument = str(re.split("\s=\s|(?<=D):\s",splited[j])[0])
							if argument == str(Matrix[0][j]):
								Matrix[i+1][j] = val				
							else:
								for k in range(len(Matrix[0,:])):
									# if str(Matrix[0][k]) == splited[j].split(" = ")[0]:	
									if argument == str(Matrix[0][k]):
										Matrix[i+1][k] = val
						except:
							print("ERROR: " + procedure_code + ": " + phenotype_center + "	" + "unable to split metadata (see impc.py split_metadata_robust")
							print("ERRORMESSAGE: " + sys.exc_info[0])
					else:
						try:
							# val = str(splited[j].split(" = ")[1])
							list_of_elements = re.split("\s=\s|(?<=D):\s",splited[j])
							argument = str(list_of_elements[0])
							list_of_elements.pop(0)
							val = ""
							for element in list_of_elements:
								element = element.replace(' =', '_')
								element = element.replace('\\,', '_')
								val += element + "_"
							val = val[:-1]
							if argument == str(Matrix[0][j]):
								Matrix[i+1][j] = val				
							else:
								for k in range(len(Matrix[0,:])):
									# if str(Matrix[0][k]) == splited[j].split(" = ")[0]:	
									if argument == str(Matrix[0][k]):
										Matrix[i+1][k] = val
						except:
							print("ERROR: " + procedure_code + ": " + phenotype_center + "	" + "unable to split metadata (see impc.py split_metadata_robust")
							print("ERRORMESSAGE: " + sys.exc_info[0])
			# get the column names
			newrow = list(df.columns.values)
			# read column names (they are lost when converted to data frame)
			initial_data = numpy.vstack([newrow, df.as_matrix()])
			# combine new metadata columns w/ original
			parsed_data = numpy.concatenate((initial_data,Matrix),axis=1)
			# metadata column index
			metadata_col_index = int(numpy.where(parsed_data[0]=='metadata')[0])
			# deletes metadata column
			parsed_data = numpy.delete(parsed_data, metadata_col_index, axis = 1)
			# creates csv file of the final data set
			parsed_df = pandas.DataFrame(data=parsed_data[1:,:], index= None, columns=parsed_data[0,:])
			if 'Start Time' in parsed_df.columns.values:
				parsed_df['Start_Date_new'] = parsed_df.loc[:,'Start Time'].apply(time_split)
				parsed_df['Start_hms_new'] = parsed_df.loc[:,'Start Time'].apply(time_split2)
			parsed_df.to_csv(output_directory + '/' + phenotype_center + '_final.csv', sep=',', index=False)
			print("COMPLETED SPLITTING METADATA: " + phenotype_center)
		else:
			print("SPLITTED METADATA: " + phenotype_center + " already exists")

# add values to DataDictionary (AAA input)
def add_to_datadict(procedure_code, pipeline = "IMPC"):
	print("STARTED ADDING TO DATA DICTIONARY")
	# output lists
	phenotype_list = []
	metadata_list = []
	add_list = []
	# input directory
	input_directory = directory + procedure_code + '\ScrapedProcedureTable'
	
	# input data
	input_data = directory + procedure_code + '\WideData'

	all_col_names = list()
	# iterate through wide files to get list of the union of the column names
	for file in os.listdir(input_data):
		fileName = input_data + '\\' + file
		data = pandas.read_csv(fileName,low_memory=False,encoding = 'ISO-8859-1')
		current_cols = data.columns.values
		all_col_names = list(set(all_col_names).union(set(current_cols)))

	# iterate through input directory
	for fn in os.listdir(input_directory):
		with open(input_directory + '\\' + fn, 'r') as csvfile:

			# read_string = csvfile.read()
			# print()
			# print(repr(read_string))
			# updated_table = re.sub("(?<!\r)\n",'',read_string)
			# print()
			# print(repr(updated_table))
			# reader = csv.reader(updated_table.split('\r\n'),delimiter='\t')
			

			# print()
			# for row in reader:
			# 	print(row)


			reader = csv.reader(csvfile, delimiter = '\t')

			table_section = 0

			for i in enumerate(reader):
				row=i[1][0].split(',')

				# value of phenotype/metadata name
				val = row[0].split(" " + pipeline)[0]
				# checks to see if in the parameter table
				if table_section ==1 and val not in phenotype_list:
					
					req_level = ""
					if row[3]=='yes':
						if row[4]=='yes':
							req_level = "req_upload_analysis"
						else:
							req_level = "req_upload"
					else:
						if row[4]=='yes':
							req_level = "req_analysis"
						else:
							req_level = "environmental variable"
					
					for j in all_col_names:
						if val in j and val!="":
							phenotype_list.append(j)
							add_list.append([j, req_level])
					
				if table_section ==2 and val not in metadata_list:
					
					req_level = ""
					if row[3]=='yes':
						if row[4]=='yes':
							req_level = "req_upload_analysis"
						else:
							req_level = "req_upload"
					else:
						if row[4]=='yes':
							req_level = "req_analysis"
						else:
							req_level = "environmental variable"

					for j in all_col_names:
						if val in j and val!="":
							metadata_list.append(j)
							add_list.append([j, req_level])
					
				if row[2] == "Type":
					table_section+=1
	# output file
	data_dict = directory + 'IMPC_DataDictionary.csv'


	# creates a DataDictionary if it doesn't exist
	# also sets the column headers of the DataDictionary
	if not os.path.isfile(data_dict):
		with open(data_dict, 'w+',newline='') as dictionary:
			dictionary_writer = csv.writer(dictionary)
			dictionary_writer.writerow(["Protocol", "Parameter", "Field", "Phenotype", "Use", "MetadataCategory"])

	# previous DataDictionary
	df = pandas.read_csv(data_dict, low_memory = False)
	add_parameter = []
	inParameterTable = True
	
	# add production_phenotype factor to DataDictionary
	add_parameter.append([procedure_code, 'production_phenotype', '', 'no', 'yes', 'req_upload_analysis'])


	for element in add_list:
		if element[0]!='':
			if element[0] in phenotype_list:
				add_parameter.append([procedure_code, element[0], 'parameter', 'yes','yes',element[1]])
			if element[0] in metadata_list:
				add_parameter.append([procedure_code, element[0], 'factor', 'no', 'yes',element[1]])

	df2 = pandas.DataFrame(add_parameter, columns = df.columns.values)
	df = df.append(df2)
	# df = df.drop_duplicates(['Protocol','Parameter'])
	df = df.drop_duplicates(['Protocol','Parameter'])

	output_file = directory + 'IMPC_DataDictionary.csv'
	df.to_csv(output_file, sep=',', index = False)
	print("FINISHED ADDING TO DATA DICTIONARY")

# scrape parameter tables (AAA_000 input)
def scrape_tables(procedure_code, rewrite = False, pipeline = "IMPC"):
	# output directory name
	project_directory = directory + procedure_code

	# create directory if it doesn't exist
	try:
		os.stat(project_directory)
	except:
		os.mkdir(project_directory)

	# output directory name
	output_directory = directory + procedure_code + '/ScrapedProcedureTable'

	# create directory if it doesn't exist
	try:
		os.stat(output_directory)
	except:
		os.mkdir(output_directory)
	files = os.listdir(output_directory)

	if rewrite or len(files) == 0: 	
		print("STARTED SCRAPING PROCEDURE TABLES")
		# get list of all protocols
		protocol_list = get_protocol_list(pipeline)
		for element in protocol_list:
			val = element[1].split("_")[1]

			if val == procedure_code:
				url = element[2]
				r = requests.get(url)

				# Turn the HTML into a Beautiful Soup object
				soup = BeautifulSoup(r.text, "lxml")

				# find the table with specified ID
				tables = soup.findAll(lambda tag: tag.name=='table' and tag.has_attr('id') and tag['id']=="paramtable")

				rows = []
				for table in tables:
					# inputs column names
					rows.append([header.text for header in table.find_all('th')])

					# rows holds all the values in the table
					for row in table.find_all('tr'):
						# values in current row
						row_val = []
						for val in row.find_all('td'):
							# if cell contains an image, add "yes" for text
							if val.img:
								row_val.append("yes")
							# if cell does not contain an image, add in the text
							else:
								row_val.append(val.text)
						rows.append(row_val)


				# write table to csv file
				with open(output_directory + '/' + element[1] +'_procedure_table.csv', 'w', newline ='') as f:
					writer = csv.writer(f)
					writer.writerows(row for row in rows if row)
		print("FINISHED SCRAPING PROCEDURE TABLES")
	else:
		print("SCRAPED PROCEDURE TABLES ALREADY EXIST")

# reliant on data dictionary
def get_pivot(procedure_code):
	

	print("STARTED CREATING PIVOT TABLE")

	# # get procedure table
	# procedure_directory = directory + procedure_code + '/ScrapedProcedureTable'

	# for procedure_table in os.listdir(procedure_directory):	
	# 	procedure_data = pandas.read_csv(procedure_directory + '/' + procedure_table, low_memory = False)

	# get DataDictionary
	IMPC_dictionary_name = directory + 'IMPC_DataDictionary.csv'

	IMPC_dictionary = pandas.read_csv(IMPC_dictionary_name, low_memory=False)

	IMPC_dictionary = IMPC_dictionary.loc[IMPC_dictionary['Protocol'] == procedure_code]
	IMPC_dictionary = IMPC_dictionary.loc[IMPC_dictionary['Phenotype'] == 'no']

	IMPC_dictionary = IMPC_dictionary.drop_duplicates(['Protocol','Parameter'])

	metadata_list = list(IMPC_dictionary['Parameter'])
	metadata_category = list(IMPC_dictionary['MetadataCategory'])
	
	# output directory name
	output_directory = directory + procedure_code + '/PivotTable'

	# create directory if it doesn't exist
	try:
		os.stat(output_directory)
	except:
		os.mkdir(output_directory)

	outputFile = output_directory + '/' + procedure_code + '_pivot.csv'
	
	input_directory = directory + procedure_code + '/WideData'

	w = len([name for name in os.listdir(input_directory) if os.path.isfile(os.path.join(input_directory,name))]) + 2 
	h = len(metadata_list) + 1
	Matrix = numpy.array([['' for x in range(w)] for y in range(h)],dtype=object)
	first_col = ["Phenotyping Centers"]
	first_col = first_col + list(metadata_list)
	for x in range(len(first_col)):
		Matrix[x][0] = first_col[x]

	second_col = ["Metadata Category"]
	second_col = second_col + list(metadata_category)

	for x in range(len(second_col)):
		Matrix[x][1] = second_col[x]


	index = 2
	# iterates by center
	for fn in os.listdir(input_directory):
		phenotype_center = fn.split('_')[0]
		# sets column header
		Matrix[0][index] = phenotype_center
		current = pandas.read_csv(input_directory + '/' + fn, low_memory = False, encoding = "ISO-8859-1")
		for metadata in range(1,len(metadata_list)+1):
			if metadata_list[metadata-1] in current.columns.values:
				# replace nan's w/ blank
				current = current.replace(numpy.nan,'', regex=True)
				# number of unique levels
				num = len(current[metadata_list[metadata-1]].unique())
				# subtracts 1 for existence of blanks
				if '' in list(current[metadata_list[metadata-1]].unique()):
					num-=1
				Matrix[metadata][index] = num
			else:
				Matrix[metadata][index] = 0
		index+=1

	df = pandas.DataFrame(data=Matrix[1:,:], index= None, columns=Matrix[0,:])

	df.to_csv(outputFile, sep=',', index=False)

	print("FINISHED CREATING PIVOT TABLE")

# given a procedure code this will download the raw data and convert it to the finalized wide data (AAA input)
def create_wide_data(procedure_code, rewrite = False, number_of_rows = 10000000, pipeline = "IMPC"):	
	
	# create directory if it doesn't exist
	try:
		os.stat(directory)
	except:
		os.mkdir(directory)

	intermediate_directory = directory + procedure_code + '/IntermediaryFiles'
	


	print("START TIME: " + time.strftime("%Y-%m-%d %H:%M:%S"))
	start_time = time.time()
	# add command for table scraping
	scrape_tables(procedure_code, rewrite, pipeline)
	# download to wide outputs
	return_procedure_group(procedure_code, number_of_rows, rewrite = rewrite)
	
	if(rewrite):
		get_and_split_by_center(procedure_code, number_of_rows)
	print("STARTED PROCESSSING COLUMNS")
	if(rewrite):
		process_col(procedure_code)
	else:
		try:
			os.stat(intermediate_directory)
		except:
			os.mkdir(intermediate_directory)

		bycenter_dir = intermediate_directory + '/ByCenter'
		try:
			os.stat(bycenter_dir)
		except:
			os.mkdir(bycenter_dir)

		if len([name for name in os.listdir(bycenter_dir) if os.path.isfile(name)]) == 0:
			get_and_split_by_center(procedure_code, number_of_rows)

		processed_dir = intermediate_directory + '/ByCenterProcessed'
		try:
			os.stat(processed_dir)
		except:
			os.mkdir(processed_dir)
		for fn in os.listdir(bycenter_dir):
			phenotyping_center = fn.split('.')[0]
			filename = phenotyping_center + "_processed.csv"
			if os.path.isfile(processed_dir + '/' + filename):
				print(": " + filename + " already exists")
			else:
				process_col(procedure_code, phenotyping_center)
	print("FINISHED PROCESSING COLUMNS")
	if(get_and_split_by_animal(procedure_code, rewrite)):
		split_metadata_robust(procedure_code)
		add_to_datadict(procedure_code, pipeline)
		get_pivot(procedure_code)
	print("END TIME: " + time.strftime("%Y-%m-%d %H:%M:%S"))
	end_time = time.time()
	elapsed_time = int(end_time - start_time)
	print("ELAPSED TIME: " + str(datetime.timedelta(seconds = elapsed_time)))

# main method
if __name__ == "__main__":
	# create_pivot_data("CBC", rewrite = True, number_of_rows = 1000)
	# split_metadata_robust("CBC", rewrite = True)
	
	# create_wide_data(procedure_code="LDT",rewrite=True,pipeline="JAX")
	# scrape_tables("TLS",rewrite=True,pipeline="JAX")
	# add_to_datadict("TLS","JAX")
	# get_pivot("LDT")
	# process_col("HBD")
	# scrape_tables(procedure_code = "TLS",rewrite = True, pipeline = "JAX")
	# add_to_datadict("URI",  "JAX")
	# create_wide_data(procedure_code="HBD",rewrite=True,pipeline="JAX")
	# create_wide_data(procedure_code="ROT",rewrite=True,pipeline="JAX")
	# # create_wide_data(procedure_code="TLS",rewrite=True,pipeline="JAX")
	# create_wide_data(procedure_code="URI",rewrite=True,pipeline="JAX")
	# add_to_datadict("INS", "IMPC")
	# get_pivot("INS")
	# create_wide_data(procedure_code="ECG",rewrite=True,pipeline="IMPC")
	# create_wide_data(procedure_code="INS",rewrite=True,pipeline="IMPC")

	# split_metadata_robust("INS", rewrite = True)
	# get_and_split_by_animal("INS", rewrite = True)
	# create_wide_data(procedure_code="HWT",rewrite=True,pipeline="IMPC")
	# create_wide_data(procedure_code="CAL",rewrite=True,pipeline="IMPC")
	#create_wide_data(procedure_code="CHL",rewrite=True,pipeline="IMPC")
	# create_wide_data(procedure_code="ABR",rewrite=True,pipeline="IMPC")
	create_wide_data(procedure_code="ALZ",rewrite=True,pipeline="IMPC")
	create_wide_data(procedure_code="BLK",rewrite=True,pipeline="IMPC")

	# get_pivot("ABR")