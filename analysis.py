"""
Tip Prediction — Simple Regression Analysis
=============================================
Steps:
 1. Load data
 2. Encode categorical variables (sex, smoker, day, time)
 3. Split into train/test sets
 4. Fit a Linear Regression model to predict `tip`
 5. Report RMSE and R²
 6. Identify the most influential variable
 7. Save a bar chart of feature importance
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------
df = pd.read_csv("tips.csv")
print("Data shape:", df.shape)
print(df.head())

# ---------------------------------------------------------
# 2. Encode categorical variables
#    One-hot encode sex, smoker, day, time (drop_first avoids redundancy)
# ---------------------------------------------------------
df_encoded = pd.get_dummies(
    df, columns=["sex", "smoker", "day", "time"], drop_first=True
)

X = df_encoded.drop(columns=["tip"])
y = df_encoded["tip"]

# ---------------------------------------------------------
# 3. Train/test split (80/20)
# ---------------------------------------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# ---------------------------------------------------------
# 4. Fit Linear Regression
# ---------------------------------------------------------
model = LinearRegression()
model.fit(X_train, y_train)

# ---------------------------------------------------------
# 5. Evaluate performance
# ---------------------------------------------------------
y_pred = model.predict(X_test)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
r2 = r2_score(y_test, y_pred)

print(f"\nModel Performance:")
print(f"RMSE: {rmse:.3f}")
print(f"R^2 : {r2:.3f}")

# ---------------------------------------------------------
# 6. Feature importance
#    Standardize X first so coefficients are comparable
#    (raw coefficients aren't comparable when features have
#     different scales, e.g. total_bill in dollars vs 0/1 dummies)
# ---------------------------------------------------------
X_std = (X - X.mean()) / X.std()
model_std = LinearRegression().fit(X_std, y)

importance = pd.Series(model_std.coef_, index=X.columns).sort_values(
    key=abs, ascending=False
)

print("\nFeature importance (standardized coefficients):")
print(importance)

top_feature = importance.index[0]
top_value = importance.iloc[0]

# ---------------------------------------------------------
# 7. Visualization — bar chart of standardized coefficients
# ---------------------------------------------------------
plt.figure(figsize=(8, 5))
colors = ["#2E86AB" if v > 0 else "#E63946" for v in importance.values]
plt.barh(importance.index[::-1], importance.values[::-1], color=colors[::-1])
plt.axvline(0, color="black", linewidth=0.8)
plt.xlabel("Standardized coefficient (effect on tip)")
plt.title("What drives the tip amount?")
plt.tight_layout()
plt.savefig("feature_importance.png", dpi=150)
print("\nSaved chart to feature_importance.png")

# ---------------------------------------------------------
# Save results for Notes.md
# ---------------------------------------------------------
with open("results.txt", "w") as f:
    f.write(f"RMSE={rmse:.3f}\n")
    f.write(f"R2={r2:.3f}\n")
    f.write(f"TOP_FEATURE={top_feature}\n")
    f.write(f"TOP_VALUE={top_value:.3f}\n")
    f.write("IMPORTANCE:\n")
    f.write(importance.to_string())
