import h2o
from h2o.estimators.deeplearning import H2ODeepLearningEstimator
import argparse
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from compare_nn_scores import density_scatter 

NUM_COLUMNS = 88

# create h2o Neural Network model for predicted zscores based on peptide encoding features
# NOTE: Start h2o cluster before running this script 

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    reqArgs = parser.add_argument_group('Required Arguments')
    reqArgs.add_argument('-i', '--input_matrix', help='Tab-delimited training matrix generated generate_nn_features.py. ' \
            f'Should be {NUM_COLUMNS} columns with column 0 being the sequence name, column 1 being the zscores, columns 2:{NUM_COLUMNS+1} being ' \
            'the features.', required=True)
    reqArgs.add_argument('-o', '--output_dir', help='Directory to output Neural Network model.', required=True)

    optArgs = parser.add_argument_group('Optional Arguments')
    optArgs.add_argument('-p', '--performance_data', action='store_true', help='Generate performance metrics using test rows. ' \
        'Note: this will significantly increase runtime due to slow generation of visuals.', required=False)

    args = parser.parse_args()

    run_h2o(args.input_matrix, args.output_dir, args.performance_data)


def run_h2o(input_matrix, output_dir, performance_data):
    h2o.init()

    data = h2o.import_file(input_matrix, header=0, sep="\t")

    os.makedirs(output_dir, exist_ok=True)
    
    train, test = data.split_frame(ratios=[0.8], seed=1234)

    print(f"Number of training rows: {train.nrow}")

    print(f"Number of test rows: {test.nrow}")

    x = data.columns[2:]
    y = data.columns[1]

    # define model with the same parameters as old h2o r script
    model = H2ODeepLearningEstimator(
        activation="TanhWithDropout",
        input_dropout_ratio=0.01,
        hidden_dropout_ratios=[0.05],
        hidden=[20],
        epochs=500
    )

    # train model
    model.train(x=x, y=y, training_frame=train)

    # Save model
    h2o.save_model(model, path=output_dir, force=True)

    if performance_data:
        
        performance = model.model_performance(test_data=test)
        print(f"R-squared: {model.r2()}\nMSE: {model.mse()}\nRMSE: {model.rmse()}")

        predictions = model.predict(test[x])
        
        # Combine predictions with actual values
        actual_vs_pred = test[y].cbind(predictions)

        # Show the first few rows
        # print(f"Predictions:\n{actual_vs_pred.head()}")

        # Convert to pandas
        df = actual_vs_pred.as_data_frame()
        df.columns = ['Actual', 'Predicted']

        x = df["Predicted"].to_numpy()
        y = df["Actual"].to_numpy()

        fig , ax = plt.subplots()
        density_scatter( x=x, y=y, fig=fig, ax=ax, bins=(1000,1000) )

        # Line of best fit
        m, b = np.polyfit(x, y, 1)
        ax.plot(x, m * x + b, color='red', label=f'Best fit: y = {m:.2f}x + {b:.2f}')

        # ax.set_aspect('equal', adjustable='box')

        plt.title('Predicted vs. Actual')
        plt.xlabel('Predicted Values')
        plt.ylabel('Actual Values')
        plt.savefig(os.path.join(output_dir, f"predicted_vs_actual.png"))


if __name__ == "__main__":
    main()