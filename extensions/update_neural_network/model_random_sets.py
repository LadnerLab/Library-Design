from generate_nn_features import gen_feature_matrix_single_lib, average_sample_counts, calculate_zscores
from h2o_nn import run_h2o
from predict_nn_scores import h2o_predict
import pandas as pd
import os
import fastatools as ft
import argparse
import random
from sklearn.metrics import r2_score, mean_squared_error
import numpy as np

# This scripts generates models for some number of random subsets of libraries and evaluates how well it predicts each indiviual library
# There are some parts that are commented out that were for evaluating the original nn model

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    reqArgs = parser.add_argument_group('Required Arguments')
    reqArgs.add_argument('-m', '--metadata', help='Tab-delimited file with columns for ' \
            'Library Name, Reads Path, NN Path, and AA Path', required=True)
    reqArgs.add_argument('-s', '--set_size', type=int, help='Number of libraries to select to create a model', required=True)
    reqArgs.add_argument('-o', '--output_dir', help='Filename to output models', required=True)

    optArgs = parser.add_argument_group('Optional Arguments')
    optArgs.add_argument('-i', '--iterations', type=int, help='Number of random models to create.', required=False)
    optArgs.add_argument('-t', '--make_test_matrices', action='store_true', help='Generate test matrices for prediction evaluation. ' \
                'Increases time, so this only needs to be done once.', required=False)

    args = parser.parse_args()

    all_libraries = pd.read_csv(args.metadata, sep="\t", index_col=0)

    # make test sets (each library)
    zscore_data = list()
    test_matrices = dict()
    test_out_dir = os.path.join(args.output_dir, "test_matrices")
    os.makedirs(test_out_dir, exist_ok=True)
    for i, row in all_libraries.iterrows():
        counts_file = row["Reads Path"]
        counts_df = pd.read_csv(counts_file, sep="\t", index_col=0)

        # average counts for sampels with the same basename
        counts_df = average_sample_counts(counts_df)
        
        # convert raw_counts to scores
        scores_df = calculate_zscores(counts_df)

        for seq_name in scores_df.index.to_list():
            for sample in scores_df.columns.to_list():
                score = scores_df.loc[seq_name, sample]
                zscore_data.append((seq_name, score))

        '''
        # do the same with the original nn data
        counts_df = pd.read_csv("/Users/scg283/Desktop/NeuralNetworkTest/nn_data/random_model_sets/original_counts_matrix.tsv", sep="\t", index_col=0)
        counts_df = average_sample_counts(counts_df)
        scores_df = calculate_zscores(counts_df)
        for seq_name in scores_df.index.to_list():
            for sample in scores_df.columns.to_list():
                score = scores_df.loc[seq_name, sample]
                zscore_data.append((seq_name, score))
        '''

        
        calculated_score_df = pd.DataFrame(zscore_data, columns=["Sequence name", "Actual Z Score"])

        matrix_out = os.path.join(test_out_dir, f"{i}_test_matrix.tsv")
        
        if args.make_test_matrices:
            counts_file = None
            nt_file = row["NN Path"]
            aa_file = row["AA Path"]

            matrix_df = gen_feature_matrix_single_lib(counts_file, nt_file, aa_file, True)

            matrix_df.to_csv(matrix_out, mode='a', sep="\t", header=False, index=False)
        
        test_matrices[i] = matrix_out
    
    # do the same with the original nn data
    # test_matrices["OriginalData"] = "/Users/scg283/Desktop/NeuralNetworkTest/nn_data/random_model_sets/original_test_matrix.tsv"


    # create a model for each iteration
    for i in range(args.iterations):
        out_i = os.path.join(args.output_dir, f"model_original")
        os.makedirs(out_i, exist_ok=True)
        set_indexes = sorted(random.sample(list(range(0, len(all_libraries))), args.set_size))

        input_set_path = make_metadata_input(all_libraries, set_indexes, out_i)

        feature_matrix_path = os.path.join(out_i, "training_feature_matrix.tsv")

        # parse through each library in metadata and generate matrix
        metadata_df = pd.read_csv(input_set_path, sep="\t", index_col=0)
        for j, row in metadata_df.iterrows():
            counts_file = row["Reads Path"]
            nt_file = row["NN Path"]
            aa_file = row["AA Path"]

            matrix_df = gen_feature_matrix_single_lib(counts_file, nt_file, aa_file, True)
            
            # append to file
            matrix_df.to_csv(feature_matrix_path, mode='a', sep="\t", header=False, index=False)

        # feature_matrix_path = "/Users/scg283/Desktop/NeuralNetworkTest/nn_data/random_model_sets/original_training_matrix.tsv"
        model_out_dir = os.path.join(out_i, f"model")
        os.makedirs(model_out_dir, exist_ok=True)
        run_h2o(feature_matrix_path, model_out_dir, True)

        model_path = find_model_path(model_out_dir)
        print(model_path)

        test_out_dir = os.path.join(out_i, "predictions")
        os.makedirs(test_out_dir, exist_ok=True)
        prediction_metrics = list()
        for name, path in test_matrices.items():
            lib_test_out_path = os.path.join(test_out_dir, f"{name}_predictions.tsv")
            h2o_predict(path, model_path, lib_test_out_path)

            predicted_df = pd.read_csv(lib_test_out_path, sep="\t", index_col=None)
            merged_df = pd.merge(predicted_df, calculated_score_df, how="inner", on="Sequence name")

            x = merged_df["Predicted Z Score"].to_numpy()
            y = merged_df["Actual Z Score"].to_numpy()

            pred = x
            actual = y
            r2 = r2_score(actual, pred)
            mse = mean_squared_error(actual, pred)
            rmse = np.sqrt(mean_squared_error(actual, pred))

            prediction_metrics.append((name, r2, mse, rmse))
        
        pd.DataFrame(
            prediction_metrics, 
            columns=["Sequence name", "R2", "MSE", "RMSE"]
        ).set_index(
            "Sequence name"
        ).to_csv(
            os.path.join(out_i, "prediction_metrics.tsv"),
            sep="\t"
        )



def make_metadata_input(all_libraries, set_indexes, out_dir):
    outfile = os.path.join(out_dir, "input_metadata_set.tsv")
    all_libraries.iloc[set_indexes].to_csv(outfile, sep="\t")
    return outfile

def find_model_path(search_path):
    for root, dirs, files in os.walk(search_path):
        for file in files[::-1]:
            if "DeepLearning_model" in file:
                return os.path.join(root, file)
    return None

if __name__ == "__main__":
    main()




