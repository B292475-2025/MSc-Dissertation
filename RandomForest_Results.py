#loading packages
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import matplotlib.pyplot as plt

#loading the data:
meta = pd.read_csv("~/Desktop/Project/Results/results_metadata.csv", sep=",")
cal_col = "CAL (ng/g feces)"
meta_indexed = meta.set_index("SampleID")

shannon = pd.read_csv("~/Desktop/Project/Results/exported-shannon/alpha-diversity.tsv", sep="\t")
shannon = shannon.rename(columns={shannon.columns[0]: "SampleID"})
#turning SampleID into the row index for the shannon table so that the dataframes can be matched by sample id
shannon_indexed = shannon.set_index("SampleID")

#alignment step to make sure that
#finds only the sample ids that exist in both tables to prevent mismatches
common_ids = meta_indexed.index.intersection(shannon_indexed.index)
print(f"Samples used: {len(common_ids)}")
#checking the presence of specific columns
predictor_cols = ["CAL (ng/g feces)", "LTF (ng/g feces)", "ALB (ng/g feces)"]
#pulling out only the matched samples and only the three predictor columns that then becomes the model's input data
X = meta_indexed.loc[common_ids, predictor_cols]
y = shannon_indexed.loc[common_ids, "shannon_entropy"]

print("\nPredictor summary:")
print(X.describe())
print("\nShannon diversity summary:")
print(y.describe())

# Check correlation between the three markers (relevant for interpreting
# feature importance, since CAL and LTF are both neutrophil-derived and
# may be correlated with each other)
print("\nCorrelation between predictors:")
print(X.corr())
#training the data (https://www.datacamp.com/tutorial/random-forest-regression)
#splitting the data into training (80%) and testing (20%)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=666
)
print(f"\nTraining set: {X_train.shape[0]} samples")
print(f"Test set: {X_test.shape[0]} samples")

#k-fold cross-validation on the training set
#training the model on the training data
reg = RandomForestRegressor(n_estimators=500, random_state=666)

cv = KFold(n_splits=5, shuffle=True, random_state=666)
cv_r2 = cross_val_score(reg, X_train, y_train, cv=cv, scoring="r2")
cv_mae = -cross_val_score(reg, X_train, y_train, cv=cv, scoring="neg_mean_absolute_error")

print(f"\nCross-validated R2 (5-fold): {cv_r2.mean():.3f} +/- {cv_r2.std():.3f}")
print(f"Individual fold R2 scores: {np.round(cv_r2, 3)}")
print(f"\nCross-validated MAE (5-fold): {cv_mae.mean():.3f} +/- {cv_mae.std():.3f}")

#fitting and training the final model:
reg.fit(X_train, y_train)
#the trained model predicting the target values for the test set and storing in
y_pred = reg.predict(X_test)

#evaluation and assessment of accuracy
#evaluation with R2 between the actual and predicted values - R2 measures the amount of variance in the target the
# model explains:
test_r2 = r2_score(y_test, y_pred)
#evaluation with Root Mean Squared Error (RMSE) between the actual and predicted values - penalises large errors more
# heavily because it squares the residuals
test_rmse = np.sqrt(mean_squared_error(y_test, y_pred))
#evaluation with Mean Absolute Error (MAE) - treats errors equally so it gives a more stable view when the data contains
# outliers
test_mae = mean_absolute_error(y_test, y_pred)

print(f"\nTest set R2: {test_r2:.3f}")
print(f"Test set RMSE: {test_rmse:.3f}")
print(f"Test set MAE: {test_mae:.3f}")
print(f"(for context, Shannon diversity range in this dataset: {y.min():.2f} - {y.max():.2f}, median {y.median():.2f})")

#feature importance - which markers matter most - shows how much each feature contributed to the model's predictions - this number sums to 1 so it represents relative importance among predictors
importances = pd.Series(reg.feature_importances_, index=predictor_cols)
#reorders series from highest importance to lowest
importances_sorted = importances.sort_values(ascending=False)

print("\nFeature importance:")
print(importances_sorted)

importances_sorted.to_csv("rf_shannon_feature_importance.csv")

#creation of plots for data visualisation:
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#creation of the plot (predicted vs actual)
sns.set_theme(style="white", font_scale=0.9)
fig, ax = plt.subplots(figsize=(5, 5))

#plot scatter points
ax.scatter(
    y_test,
    y_pred,
    color="#34495E",
    s=60,
    edgecolor="black",
    linewidth=0.8,
    zorder=3,
)

#plot perfect prediction line (y = x)
min_val = min(min(y_test), min(y_pred)) - 0.1
max_val = max(max(y_test), max(y_pred)) + 0.1
ax.plot(
    [min_val, max_val],
    [min_val, max_val],
    color="grey",
    linestyle="--",
    linewidth=1.5,
    label="Perfect Fit (y = x)",
)

ax.set_xlim(min_val, max_val)
ax.set_ylim(min_val, max_val)
ax.set_xlabel("Actual Shannon Diversity", fontweight="bold")
ax.set_ylabel("Predicted Shannon Diversity", fontweight="bold")
ax.set_title(
    "Random Forest Regression: Predicted vs. Actual",
    fontweight="bold",
    pad=12,
    loc="left",
)

stats_text = f"$R^2$ = {test_r2:.2f}\nMAE = {test_mae:.2f}"
ax.text(
    0.79,
    0.07,
    stats_text,
    transform=ax.transAxes,
    fontsize=9,
    verticalalignment="bottom",
    horizontalalignment="left",
    bbox=dict(boxstyle="square,pad=0.5", facecolor="white", alpha=0.8, edgecolor="#CCCCCC"),
)

sns.despine()
ax.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()

plt.savefig("predicted_vs_actual.png", dpi=300, bbox_inches="tight")
plt.savefig("predicted_vs_actual.pdf", bbox_inches="tight")
plt.close()







#plotting feature importance bar chart (https://stackoverflow.com/questions/44101458/random-forest-feature-importance-chart-using-python)
sns.set_theme(style="white", font_scale=0.9)
#mapping abbreviations to full names
label_mapping = {
    "CAL (ng/g feces)": "Calprotectin (ng/g faeces)",
    "ALB (ng/g feces)": "Albumin (ng/g faeces)",
    "LTF (ng/g feces)": "Lactoferrin (ng/g faeces)",
}
#rename the index before sorting and plotting
sorted_series = importances_sorted.rename(index=label_mapping).sort_values(ascending=True)
#adjust figure sizes based on the number of features
num_features = len(sorted_series)
fig, ax = plt.subplots(figsize=(6, max(3, num_features * 0.35)))
#horizontal bars with a clean aesthetic
bars = ax.barh(
    y=sorted_series.index,
    width=sorted_series.values,
    color="#34495E",  # Modern slate blue/grey hue
    edgecolor="none",
    height=0.65,
)
#cleaning up the axes
sns.despine(top=True, right=True, left=True, bottom=False)
#adding faint gridlines onto the graph
ax.xaxis.grid(True, linestyle="--", alpha=0.5, color="#CCCCCC")
ax.set_axisbelow(True)

#titles and typography
ax.set_xlabel("Mean Decrease in Impurity", fontsize=10, fontweight="bold", labelpad=8)
ax.set_title(
    "Predicting Shannon Diversity: Marker Importance",
    fontsize=11,
    fontweight="bold",
    pad=12,
    loc="left",
)

#expanding teh x-limit so that value labels dont get cut off
ax.set_xlim(0, sorted_series.max() * 1.15)
plt.tight_layout()
#saving the output
plt.savefig("second_shannon_feature_importance.png", dpi=300, bbox_inches="tight")
plt.savefig("second_shannon_feature_importance.pdf", bbox_inches="tight")
plt.close()

#mean decrease in impurity = Impurity is a measure of how mixed up or uncertain the predictions are at a given node in a tree.













