# Symbolic Regression GECCO Competition - 2023 - Track 2 - Interpretability

The participants will be free to experiment with this data set until the deadline. 
Analysis on this dataset should include the production of one or more models and a detailed pre and post-analysis regarding the interpretation of that model.
Together with the data set we will also provide the context (what the data is about, how it was extracted, etc.) and a description of each feature.
The interpretability analysis in the extended abstract can contain any information that can be extracted from the symbolic expressions. 
For example, you can try to analyze the behavior of the target value w.r.t. certain features, make a study of how some features interact, measure the uncertainty of the predictions or confidence intervals, and explain whether these results are reasonable given the nature of the data. 
Extra points will be awarded for analysis that is unique to Symbolic Regression models.

At the competition submission deadline the repository of a participating team should contain:

- [**required**] A file called `dataset_best_models` containing all the relevant models created as **sympy-compatible expressions**.
- [**required**] A 4 page extended abstract in PDF format describing the algorithm, pipeline, and the intepretability analysis of the real-world data set (`paper` folder). This PDF must contain the name and contact information of the participants.
- [to be eligible for prize] Reproducibility documentation in the `src` folder.
    - Installation and execution instructions 
    - Scripts (Jupyter Notebook is accepted) with detailed instructions of all steps used to produce models (including hyperparameters search, if any) 
    - Code or instruction to compute any additonal analysis, if applicable.

## Evaluation criteria

Each member of the jury will assign a score to each submission and the final score will be a simple average of the assigned scores. The jury will take into consideration:

- Level of details in the pipeline
- Readability of the model
- Interestingness of the pre and post analysis process
- Analysis of interpretation (with special points for analysis that can only be made using SR models)

Notice that the scores are subjective and these criteria are only a guideline to the jury (e.g., participants who provide a large model with a very good interpretation may score more points than participants that provide a small model with a less good interpretation).

## Data set: Montreal Bike Lane

This data set contains information of the number of bikes that crossed a certain bike lane in Montreal in the year 2015. The goal of the regression model is to predict the number of bikes with a model $f(date, lane, x)$ where date and lane are the date and lane of interest, and $x$ is a vector of any other feature you find relevant.
The file `MontrealBikeLane.csv` contains the data with the information of the count of bikes at every lane in different dates, and the file `WeatherInfo.csv` contains weather information for each of the dates in the first data set.

You are allowed to include more information from external data as long as it is described in the report.

## Repository License

The repositories will be kept private during the whole competition and it will become open to the public **after** GECCO 2023 conference with a BSD-3 license. Please, make sure that you only keep files conforming to such license.

## Deadline

01 June 2023, 00:00 anywhere in the world.

## Question and issues

Any questions and issues can be addressed to folivetti@ufabc.edu.br or at our Discord server (https://discord.gg/Dahqh3Chwy)