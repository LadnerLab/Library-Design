These scipts for generating a neural network model for predicting a peptide's z score based on nucleotides, amino acids, and codons

Example metadata formatting:
Library Name	Reads Path	NN Path	AA Path
CT64	IM0182_CT64unq65_raw_3mm_SBonly_CS.tsv	CT64_u65.fna	CT64.faa
CWP IM0156_IM0164_IM0175_CWP_86nt_3mm_i1mm_SblkOnly_CS.tsv	CWP_u40.fna	CWP.faa
FAV FAV/FAV-Uq_3mm_CS.tsv	FAV_coded.fna	FAV_coded.faa

The reads path is optional when using generate_nn_features.py, it is not required when making a matrix for predicting.

# Scripts
generate_nn_features.py - Generate the input matrix for training or predicting using the model
h2o_nn.py - Generate the Neural Network model using the training feature matrix (input matrix includes z scores)
predict_nn_scores.py - Generate z score predictions of each peptide is a prediction feature matrix using the model
model_random_sets.py - Generate multiple models and their performance metrics using random subsets of the input metadata