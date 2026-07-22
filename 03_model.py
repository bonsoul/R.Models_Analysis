"""
03_model.py
Train a linear regression model to predict tip amount and report performance.
"""
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

df = pd.read_csv("tips_clean.csv")

X = pd.get_dummies(
    df[["total_bill", "size", "sex", "smoker", "day", "time"]],
    drop_first=True,
)
y = df["tip"]

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

model = LinearRegression().fit(X_train, y_train)
y_pred = model.predict(X_test)

rmse = np.sqrt(mean_squared_error(y_test, y_pred))
mae = mean_absolute_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print("=== Model performance (test set, n=%d) ===" % len(y_test))
print(f"RMSE: {rmse:.3f}")
print(f"MAE:  {mae:.3f}")
print(f"R^2:  {r2:.3f}")

print("\n=== Raw coefficients (unstandardized, $ per unit) ===")
for feat, coef in zip(X.columns, model.coef_):
    print(f"{feat:>12}: {coef:+.4f}")
print(f"{'intercept':>12}: {model.intercept_:+.4f}")

# Baseline comparison: predicting the mean tip every time
baseline_pred = np.full_like(y_test, y_train.mean(), dtype=float)
baseline_rmse = np.sqrt(mean_squared_error(y_test, baseline_pred))
print(f"\nBaseline (predict mean) RMSE: {baseline_rmse:.3f}")
print(f"Model improves RMSE by {100*(1 - rmse/baseline_rmse):.1f}% over baseline")

# Save results for Notes.md reference
with open("model_results.txt", "w") as f:
    f.write(f"RMSE: {rmse:.3f}\nMAE: {mae:.3f}\nR2: {r2:.3f}\n")
    f.write(f"Baseline RMSE: {baseline_rmse:.3f}\n")
